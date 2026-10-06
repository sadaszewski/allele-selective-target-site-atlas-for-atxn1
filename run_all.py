#!/usr/bin/env python3
"""Reproduce every number in the ATXN1 allele-selective atlas manuscript.

Usage:  python run_all.py [--data-dir data] [--results-dir results] [--skip-htt]

Downloads ~2.1 GB on first run and caches it. Runtime after download is a few
minutes on one core.
"""
import argparse, json, sys, time
from pathlib import Path
import numpy as np, pandas as pd

from atxn1atlas import config as C
from atxn1atlas import data, atlas, eligibility, report
from atxn1atlas.pam import NUCLEASES, pam_states, protospacer, revcomp


def build_locus(locus, build, data_dir, flank, min_het, verbose=True):
    """Scan one gene (+/- flank) on the chosen build; return what follows needs."""
    b = C.BUILDS[build]
    chrom = locus["chrom"]
    k0, k1 = b["coord_keys"]
    start, end = locus[k0] - flank, locus[k1] + flank
    vcf = data.ensure_vcf(chrom, build, data_dir)
    panel = data.load_panel(data.ensure_panel(data_dir))
    # fetch a little extra sequence so PAM windows near the edges are complete
    pad = 50
    seq = data.ensure_reference(b["acc"][chrom], start - pad, end + pad, data_dir)
    seq_offset = start - pad
    contig = b["contig"].format(c=chrom)
    df, haps, meta = atlas.scan_window(vcf, contig, start, end, seq, seq_offset,
                                       panel, min_het=min_het, verbose=verbose)
    if len(meta["samples"]) != 2504:
        raise RuntimeError(f"kept {len(meta['samples'])} samples, expected the 2,504 unrelated")
    return dict(locus=locus, build=build, df=df, haps=haps, meta=meta, seq=seq,
                seq_offset=seq_offset, window=(start, end))


def add_build_positions(df, locus, build):
    """Add explicit pos_grch37 / pos_grch38 columns.

    The analysed build gives real positions; the other is derived by the
    gene's constant inter-build offset (asserted equal at both termini).
    """
    off = C.locus_offset(locus)          # GRCh37 - GRCh38
    if C.BUILDS[build]["assembly"] == "GRCh38":
        df["pos_grch38"] = df["pos"]
        df["pos_grch37"] = df["pos"] + off
    else:
        df["pos_grch37"] = df["pos"]
        df["pos_grch38"] = df["pos"] - off
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default=str(C.DATA_DIR))
    ap.add_argument("--results-dir", default=str(C.RESULTS_DIR))
    ap.add_argument("--skip-htt", action="store_true",
                    help="skip the HTT control (saves a large download)")
    ap.add_argument("--build", default=C.DEFAULT_BUILD, choices=sorted(C.BUILDS),
                    help="which phased callset to analyse")
    a = ap.parse_args()
    dd, rd = Path(a.data_dir), Path(a.results_dir)
    rd.mkdir(parents=True, exist_ok=True)
    R = {}
    t0 = time.time()
    bprof = C.BUILDS[a.build]
    asm = bprof["assembly"]
    R["build"] = a.build
    R["build_label"] = bprof["label"]
    R["assembly"] = asm

    # ---------------------------------------------------------------- 1. atlas
    report.header(f"1. BUILDING THE ATLAS  (ATXN1 +/- 200 kb, {asm})\n   {bprof['label']}")
    A = build_locus(C.ATXN1, a.build, dd, C.FLANK, min_het=0.0)
    raw = A["df"]
    print(f"  records in window          : {A['meta']['n_records']:,}")
    print(f"  biallelic SNVs             : {len(raw):,}")
    clean, n_dropped, n_impossible = atlas.hwe_filter(raw)
    print(f"  dropped (het > {C.MAX_HET} or z > {C.MAX_HWE_Z}) : {n_dropped}"
          f"   (of which het > HWE ceiling: {n_impossible})")
    print(f"  retained                   : {len(clean):,}")
    R["n_records"] = A["meta"]["n_records"]
    R["n_biallelic_snv"] = len(raw)
    R["n_clean"] = len(clean)

    keep_idx = raw.index.isin(clean.index) if False else None
    # realign haplotypes to the filtered rows
    mask = ~raw["pos"].isin(set(raw["pos"]) - set(clean["pos"]))
    haps_clean = A["haps"][mask.values]
    assert len(haps_clean) == len(clean)

    # reference concordance (scan_window already rejects mismatches; prove it)
    seq, off = A["seq"], A["seq_offset"]
    mism = sum(1 for p, r in zip(clean["pos"], clean["ref"]) if seq[p - off] != r)
    print(f"  REF-allele concordance     : {len(clean) - mism:,}/{len(clean):,}")
    assert mism == 0, "reference mismatch"
    R["ref_concordance"] = f"{len(clean) - mism}/{len(clean)}"

    rows = [{"nuclease": n,
             "altering": int(clean[f"{n}_alter"].sum()),
             "created": int((clean[f"{n}_created"] > 0).sum()),
             "destroyed": int((clean[f"{n}_destroyed"] > 0).sum())}
            for n, _, _, _ in NUCLEASES]
    rows.append({"nuclease": "ANY", "altering": int(clean["any_pam_alter"].sum()),
                 "created": None, "destroyed": None})
    report.table(rows, ["nuclease", "altering", "created", "destroyed"],
                 "PAM-altering counts, full window")
    R["pam_counts_window"] = {r["nuclease"]: r["altering"] for r in rows}

    k0, k1 = bprof["coord_keys"]
    g0, g1 = C.ATXN1[k0], C.ATXN1[k1]
    in_gene = (clean["pos"] >= g0) & (clean["pos"] <= g1)
    gene = clean[in_gene]
    R["in_gene_sp_ge20"] = int(((gene["hetf_ALL"] >= .20) & (gene["SpCas9_alter"] == 1)).sum())
    R["in_gene_sp_ge30"] = int(((gene["hetf_ALL"] >= .30) & (gene["SpCas9_alter"] == 1)).sum())
    hmin = gene[[f"hetf_{k}" for k in C.SUPERPOPS]].min(axis=1)
    R["in_gene_sp_allpops_ge20"] = int(((hmin >= .20) & (gene["SpCas9_alter"] == 1)).sum())
    R["tier1_n"] = int(((hmin >= .30) & (gene["SpCas9_alter"] == 1)).sum())
    print(f"\n  in-gene SpCas9 sites, het >= 20%            : {R['in_gene_sp_ge20']:,}")
    print(f"  in-gene SpCas9 sites, het >= 30%            : {R['in_gene_sp_ge30']:,}")
    print(f"  in-gene SpCas9, het >= 20% in ALL 5 superpops: {R['in_gene_sp_allpops_ge20']:,}")
    print(f"  tier 1 (>= 30% in all 5 superpops)          : {R['tier1_n']:,}")

    # repeat tract, located from sequence alone
    hits = data.locate_repeat(seq, off)
    _o = C.locus_offset(C.ATXN1) if asm == "GRCh38" else 0
    R["repeat_tract"] = [int(hits[0][1]), int(hits[0][2])]          # analysed build
    R["repeat_tract_grch37"] = [int(hits[0][1]) + _o, int(hits[0][2]) + _o]
    R["repeat_tract_grch38"] = [int(hits[0][1]) + _o - C.locus_offset(C.ATXN1),
                                int(hits[0][2]) + _o - C.locus_offset(C.ATXN1)]
    print(f"\n  CAG repeat region (CTG on the plus strand, ATXN1 is minus-strand):")
    print(f"    chr6:{hits[0][1]:,}-{hits[0][2]:,}   longest pure run {hits[0][0]} units")
    print(f"    plus strand: {seq[hits[0][1]-off:hits[0][2]-off+1]}")
    print(f"    minus strand (gene orientation): "
          f"{revcomp(seq[hits[0][1]-off:hits[0][2]-off+1])}")

    # write the atlas
    out = add_build_positions(clean.copy(), C.ATXN1, a.build)
    out["region"] = np.where(in_gene, "ATXN1_gene", "flank_200kb")
    out["het_min_superpop"] = clean[[f"hetf_{k}" for k in C.SUPERPOPS]].min(axis=1)
    cols = (["rsid", "pos_grch37", "pos_grch38", "ref", "alt", "region"]
            + [f"hetf_{k}" for k in ["ALL"] + C.SUPERPOPS] + ["het_min_superpop"]
            + [f"{n}_{s}" for n, _, _, _ in NUCLEASES for s in ("created", "destroyed")]
            + ["any_pam_alter"])
    out[cols].sort_values("pos_grch37").to_csv(rd / "ATXN1_allele_selective_atlas_v0.csv.gz", index=False)
    tier1 = out[in_gene & (out["SpCas9_alter"] == 1) & (out["het_min_superpop"] >= .30)]
    tier1[cols].sort_values("het_min_superpop", ascending=False).to_csv(
        rd / "ATXN1_tier1_candidates.csv", index=False)

    # ------------------------------------------------- 2. HTT control + compare
    report.header("2. VALIDATION AGAINST HTT, AND THE HEAD-TO-HEAD COMPARISON")
    comp = []
    for locus, skip in ((C.HTT, a.skip_htt), (C.ATXN1, False)):
        if skip:
            print("  (HTT skipped)")
            continue
        B = build_locus(locus, a.build, dd, 0, min_het=C.BENCH_MIN_HET, verbose=False)
        d = B["df"]
        d = d[d["SpCas9_alter"] == 1].reset_index(drop=True)
        h = B["haps"][(B["df"]["SpCas9_alter"] == 1).values]
        cr = d["SpCas9_created"].values > 0
        de = d["SpCas9_destroyed"].values > 0
        e = eligibility.phase_agnostic(d, h, cr, de, B["meta"]["superpop"])
        comp.append({"gene": locus["name"], "n_sites": len(d),
                     "best_het_EUR": float(d["hetf_EUR"].max()),
                     "het_any_EUR": e["het_any_EUR"], "het_any_ALL": e["het_any_ALL"],
                     "pa_EUR": e["phase_agnostic_EUR"], "pa_ALL": e["phase_agnostic_ALL"]})
        if locus["name"] == "HTT":
            r = C.find_variant(d, "rs363099", asm)
            if r is not None:
                print(f"  rs363099  chr4:{int(r['pos']):,} ({asm})  {r['ref']}>{r['alt']}  "
                      f"created={int(r['SpCas9_created'])} destroyed={int(r['SpCas9_destroyed'])}")
                print(f"    EUR heterozygosity {r['hetf_EUR']*100:.1f}%  x 0.5 phase = "
                      f"{r['hetf_EUR']*50:.1f}%   (published ~20%)")
                R["rs363099"] = {"pos": int(r["pos"]),
                                 "het_EUR": float(r["hetf_EUR"]),
                                 "single_site_elig_EUR": float(r["hetf_EUR"] / 2),
                                 "destroys": int(r["SpCas9_destroyed"])}
    report.table(comp, ["gene", "n_sites", "best_het_EUR", "het_any_EUR",
                        "pa_EUR", "pa_ALL"],
                 "SpCas9, gene body, het >= 2%",
                 pct=("best_het_EUR", "het_any_EUR", "pa_EUR", "pa_ALL"))
    R["comparison"] = comp

    # -------------------------------------------------------- 3. guide panel
    report.header("3. MINIMAL GUIDE PANEL (greedy cover)")
    cand = clean[in_gene & (clean["SpCas9_alter"] == 1)
                 & (clean["hetf_ALL"] >= C.PANEL_MIN_HET)].reset_index(drop=True)
    cand_mask = (in_gene & (clean["SpCas9_alter"] == 1)
                 & (clean["hetf_ALL"] >= C.PANEL_MIN_HET)).values
    cand_haps = haps_clean[cand_mask]
    print(f"  candidate sites: {len(cand):,}")
    steps = eligibility.greedy_panel(cand, cand_haps, A["meta"]["superpop"])
    report.table(steps, ["rsid", "pos", "cumulative"] + [f"cum_{k}" for k in C.SUPERPOPS],
                 None, pct=["cumulative"] + [f"cum_{k}" for k in C.SUPERPOPS])
    R["panel"] = steps
    R["n_panel_candidates"] = len(cand)

    # ---------------------------------------------------------- 4. founder model
    report.header("4. FOUNDER MODEL")
    P, H, cr, de = atlas.pam_profile(cand, cand_haps, "SpCas9")
    sup_h = np.concatenate([A["meta"]["superpop"]] * 2)
    pop_h = np.concatenate([A["meta"]["pop"]] * 2)
    sw = eligibility.founder_sweep(P, sup_h)
    report.table(sw, ["group", "n_hap", "worst", "p1", "median",
                      "distinct_profiles", "commonest_profile_frac"],
                 "Founder and wild-type pool drawn from the same SUPERPOPULATION",
                 pct=("worst", "p1", "median", "commonest_profile_frac"))
    pw = eligibility.founder_sweep(P, pop_h)
    pw_sorted = sorted(pw, key=lambda r: r["worst"])
    report.table(pw_sorted, ["group", "n_hap", "worst", "median", "distinct_profiles"],
                 "Both drawn from a single POPULATION (the bottleneck scenario)",
                 pct=("worst", "median"))
    floor = min(r["worst"] for r in pw)
    print(f"\n  GLOBAL WORST CASE across all {len(pw)} populations: {floor*100:.2f}%"
          f"  ({min(pw, key=lambda r: r['worst'])['group']})")
    R["founder_superpop"] = sw
    R["founder_population"] = pw
    R["founder_floor"] = floor
    pd.DataFrame(pw).to_csv(rd / "founder_worstcase_by_population.csv", index=False)

    # ------------------------------------------------------------ 5. rs2075974
    report.header("5. rs2075974")
    r = C.find_variant(clean, "rs2075974", asm)
    if r is not None:
        i = int(r["pos"]) - off
        created = pam_states(seq, i, r["alt"], "NGG") - pam_states(seq, i, r["ref"], "NGG")
        strand, pstart = sorted(created)[0]
        ps = protospacer(seq, i, r["alt"], pstart, strand)
        p37 = int(r["pos"]) + (C.locus_offset(C.ATXN1) if asm == "GRCh38" else 0)
        p38 = p37 - C.locus_offset(C.ATXN1)
        print(f"  chr6:{p37:,} (GRCh37) / {p38:,} (GRCh38)   {r['ref']}>{r['alt']}")
        print(f"  distance to CAG tract : {abs(int(r['pos']) - hits[0][1])} bp")
        print(f"  SpCas9 created/destroyed: {int(r['SpCas9_created'])}/{int(r['SpCas9_destroyed'])}"
              f"   PAM on the {strand} strand")
        print(f"  protospacer (alt)     : 5'-{ps}-3'")
        print(f"\n  Eligibility if the expansion sits on the {r['alt']} chromosome"
              f" (= 1 - freq({r['alt']})):")
        elig = {}
        for k in C.SUPERPOPS + ["ALL"]:
            f = float(r[f"af_{k}"])
            elig[k] = 1 - f
            print(f"    {k:<4} freq={f*100:5.1f}%   eligible={100*(1-f):5.1f}%"
                  f"   het={r[f'hetf_{k}']*100:5.1f}%")
        R["rs2075974"] = {"pos_grch37": p37, "pos_grch38": p38,
                          "ref": r["ref"], "alt": r["alt"],
                          "created": int(r["SpCas9_created"]),
                          "destroyed": int(r["SpCas9_destroyed"]),
                          "strand": strand, "protospacer": ps,
                          "het_ALL": float(r["hetf_ALL"]),
                          "eligibility": elig}
    else:
        print("  rs2075974 not found in the clean set")

    R["runtime_s"] = round(time.time() - t0, 1)
    (rd / "results.json").write_text(json.dumps(R, indent=2, default=float))
    report.header(f"DONE in {R['runtime_s']}s.  Results written to {rd}/")
    for f in sorted(rd.iterdir()):
        print(f"  {f.name}  ({f.stat().st_size:,} bytes)")
    return R


if __name__ == "__main__":
    main()
