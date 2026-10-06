#!/usr/bin/env python3
"""Add rsIDs to a 30x atlas by joining against a Phase 3 atlas.

The 30x NYGC callset's ID column holds "chrom:pos:ref:alt" strings, not rsIDs.
Phase 3 carries rsIDs. Both atlases report pos_grch37, so the join is exact on
(pos_grch37, ref, alt) -- no coordinate conversion and no network needed.

Usage:
  python annotate_rsids.py results_30x/ATXN1_allele_selective_atlas_v0.csv.gz \
                           results_phase3/ATXN1_allele_selective_atlas_v0.csv.gz
"""
import sys
from pathlib import Path
import pandas as pd


def main(target, source, out=None):
    t = pd.read_csv(target)
    s = pd.read_csv(source, usecols=["rsid", "pos_grch37", "ref", "alt"])
    s = s[s["rsid"].astype(str).str.startswith("rs")]
    key = ["pos_grch37", "ref", "alt"]
    merged = t.merge(s.rename(columns={"rsid": "rsid_phase3"}), on=key, how="left")
    n = merged["rsid_phase3"].notna().sum()
    print(f"{Path(target).name}: {len(t):,} rows, {n:,} matched to an rsID "
          f"({n/len(t)*100:.1f}%)")
    merged["rsid"] = merged["rsid_phase3"].fillna(merged["rsid"])
    merged = merged.drop(columns=["rsid_phase3"])
    out = out or str(target).replace(".csv", "_rsid.csv")
    merged.to_csv(out, index=False)
    print(f"wrote {out}")
    miss = merged[~merged["rsid"].astype(str).str.startswith("rs")]
    if len(miss):
        print(f"note: {len(miss):,} rows keep their positional ID "
              f"(not present in the Phase 3 callset)")
    return merged


if __name__ == "__main__":
    if not 2 <= len(sys.argv) - 1 <= 3:
        sys.exit(__doc__)
    main(*sys.argv[1:])
