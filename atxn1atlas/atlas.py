"""Build the PAM-altering variant atlas across a locus."""
import numpy as np, pandas as pd, pysam
from . import config as C
from .pam import NUCLEASES, screen_variant


def _genotype_matrix(rec, n_samples, sidx):
    gt = [rec.samples[i]["GT"] for i in range(n_samples)]
    return np.array([[gt[i][0], gt[i][1]] for i in sidx], dtype=np.int8)


def scan_window(vcf_path, chrom, start, end, seq, seq_offset, panel,
                min_het=0.0, require_pam=None, verbose=True):
    """Scan a window for biallelic SNVs; return (DataFrame, haplotype array, sample meta).

    haplotype array is sites x samples x 2 for the retained rows, in row order.
    """
    vf = pysam.VariantFile(str(vcf_path))
    samples = list(vf.header.samples)
    # Phase 3 annotates variant class in VT; the 30x NYGC callset uses
    # VariantType, or neither. Querying an absent INFO key raises, so look once.
    vt_key = next((k for k in ("VT", "VariantType") if k in vf.header.info), None)
    keep = [s for s in samples if s in panel]
    sidx = np.array([samples.index(s) for s in keep])
    sup = np.array([panel[s][1] for s in keep])
    pop = np.array([panel[s][0] for s in keep])

    rows, haps = [], []
    n_records = 0
    for rec in vf.fetch(chrom, start - 1, end):
        n_records += 1
        if len(rec.alts or ()) != 1:
            continue
        ref, alt = rec.ref, rec.alts[0]
        if len(ref) != 1 or len(alt) != 1 or ref not in "ACGT" or alt not in "ACGT":
            continue
        if vt_key is not None:
            vt = rec.info.get(vt_key)
            if vt is not None:
                vals = vt if isinstance(vt, tuple) else (vt,)
                # keep only single-nucleotide classes; ref/alt lengths are
                # already checked above, so this is belt-and-braces
                if not any("SNP" in str(v).upper() or "SNV" in str(v).upper()
                           for v in vals):
                    continue
        i = rec.pos - seq_offset
        if not (0 <= i < len(seq)) or seq[i] != ref:
            continue                      # reference mismatch -> skip (counted below)

        a = _genotype_matrix(rec, len(samples), sidx)
        het = a[:, 0] != a[:, 1]
        if het.mean() < min_het:
            continue

        d = {"pos": rec.pos, "rsid": rec.id or ".", "ref": ref, "alt": alt,
             "het_ALL": int(het.sum()), "n_ALL": len(keep),
             "af_ALL": float(a.sum()) / (2 * len(keep))}
        for k in C.SUPERPOPS:
            m = sup == k
            d[f"het_{k}"] = int(het[m].sum())
            d[f"n_{k}"] = int(m.sum())
            d[f"af_{k}"] = float(a[m].sum()) / (2 * m.sum())

        pam_any = False
        for name, motif, _, _ in NUCLEASES:
            cr, de = screen_variant(seq, i, ref, alt, motif)
            d[f"{name}_created"], d[f"{name}_destroyed"] = cr, de
            d[f"{name}_alter"] = int(bool(cr or de))
            pam_any |= bool(cr or de)
        d["any_pam_alter"] = int(pam_any)

        if require_pam and not d[f"{require_pam}_alter"]:
            continue
        rows.append(d)
        haps.append(a)
        if verbose and len(rows) % 5000 == 0:
            print(f"    {len(rows):,} variants kept / {n_records:,} records", flush=True)

    df = pd.DataFrame(rows)
    for k in ["ALL"] + C.SUPERPOPS:
        df[f"hetf_{k}"] = df[f"het_{k}"] / df[f"n_{k}"]
    meta = dict(samples=keep, superpop=sup, pop=pop, n_records=n_records)
    return df, np.array(haps), meta


def hwe_filter(df):
    """Drop sites whose heterozygosity exceeds the Hardy-Weinberg ceiling."""
    p = df["af_ALL"]
    exp = 2 * p * (1 - p)
    se = np.sqrt(np.maximum(exp * (1 - exp), 1e-9) / df["n_ALL"])
    z = (df["hetf_ALL"] - exp) / np.maximum(se, 1e-9)
    bad = (df["hetf_ALL"] > C.MAX_HET) | (z > C.MAX_HWE_Z)
    return df[~bad].reset_index(drop=True), int(bad.sum()), int((df["hetf_ALL"] > C.MAX_HET).sum())


def pam_profile(df, haps, nuclease="SpCas9"):
    """P[i,h] = 1 where haplotype h carries the PAM-bearing allele at site i.

    The PAM-bearing allele is ALT where the variant creates a PAM, REF where it
    destroys one; a variant may qualify both ways at different offsets.
    """
    created = df[f"{nuclease}_created"].values > 0
    destroyed = df[f"{nuclease}_destroyed"].values > 0
    H = np.concatenate([haps[:, :, 0], haps[:, :, 1]], axis=1)
    P = ((H == 1) & created[:, None]) | ((H == 0) & destroyed[:, None])
    return P, H, created, destroyed
