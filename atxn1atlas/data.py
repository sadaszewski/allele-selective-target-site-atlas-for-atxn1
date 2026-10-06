"""Fetch and cache the public inputs. Nothing here needs credentials."""
import subprocess, sys
from pathlib import Path
import numpy as np
from . import config as C


def _curl(url: str, dest: Path, timeout: int = 1800) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    print(f"  downloading {dest.name} ...", flush=True)
    cmd = ["curl", "-sS", "-L", "--retry", "8", "--retry-all-errors",
           "--retry-delay", "5", "-o", str(dest), url]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    if r.returncode != 0 or not dest.exists() or dest.stat().st_size == 0:
        dest.unlink(missing_ok=True)
        raise RuntimeError(f"download failed: {url}\n{r.stderr}")
    return dest


def ensure_panel(data_dir: Path = None) -> Path:
    d = Path(data_dir or C.DATA_DIR)
    return _curl(C.PANEL_URL, d / "integrated_call_samples_v3.20130502.ALL.panel")


def ensure_vcf(chrom: str, build: str = None, data_dir: Path = None) -> Path:
    """Fetch the phased callset for one chromosome, for the chosen build."""
    d = Path(data_dir or C.DATA_DIR)
    b = C.BUILDS[build or C.DEFAULT_BUILD]
    url = b["vcf_url"].format(c=chrom)
    vcf = _curl(url, d / url.rsplit("/", 1)[-1])
    _curl(url + ".tbi", Path(str(vcf) + ".tbi"))
    return vcf


def ensure_reference(acc: str, start: int, end: int, data_dir: Path = None) -> str:
    """Plus-strand reference sequence for acc:start-end, 1-based inclusive."""
    d = Path(data_dir or C.DATA_DIR)
    cache = d / f"ref_{acc}_{start}_{end}.seq"
    if cache.exists():
        seq = cache.read_text().strip()
    else:
        url = (f"{C.EUTILS}/efetch.fcgi?db=nuccore&id={acc}&rettype=fasta"
               f"&retmode=text&seq_start={start}&seq_stop={end}")
        fa = _curl(url, d / f"ref_{acc}_{start}_{end}.fa")
        seq = "".join(l.strip() for l in fa.read_text().splitlines()
                      if not l.startswith(">")).upper()
        cache.write_text(seq)
    if len(seq) != end - start + 1:
        raise RuntimeError(f"reference length {len(seq)} != expected {end - start + 1}")
    return seq


def load_panel(path: Path = None):
    """-> dict sample -> (population, superpopulation) for the 2,504 unrelated samples."""
    p = Path(path) if path else ensure_panel()
    panel = {}
    for i, line in enumerate(p.read_text().splitlines()):
        f = line.split()
        if i == 0 or len(f) < 3:
            continue
        panel[f[0]] = (f[1], f[2])
    if len(panel) != 2504:
        raise RuntimeError(f"panel has {len(panel)} samples, expected 2504")
    return panel


def locate_repeat(seq: str, offset: int, min_units: int = 6, max_gap: int = 12):
    """Locate the CAG repeat region on the plus strand.

    ATXN1 is minus-strand, so its CAG tract reads as CTG here. The wild-type
    tract is interrupted (CAT on the gene strand = ATG here), which splits it
    into several contiguous runs; runs separated by <= max_gap bases (a CAT-CAG-CAT interruption is 9) are merged
    so the returned region spans the whole interrupted repeat.

    Returns a list of (units_in_longest_pure_run, region_start, region_end),
    1-based inclusive, sorted by region length descending.
    """
    import re
    runs = []
    for pat in (r"(?:CTG){%d,}" % min_units, r"(?:CAG){%d,}" % min_units):
        for m in re.finditer(pat, seq):
            runs.append([m.start(), m.end(), len(m.group(0)) // 3])
    if not runs:
        return []
    runs.sort()
    merged = [runs[0]]
    for s, e, u in runs[1:]:
        if s - merged[-1][1] <= max_gap:
            merged[-1][1] = max(merged[-1][1], e)
            merged[-1][2] = max(merged[-1][2], u)
        else:
            merged.append([s, e, u])
    out = [(u, offset + s, offset + e - 1) for s, e, u in merged]
    out.sort(key=lambda r: (r[2] - r[1], r[0]), reverse=True)
    return out
