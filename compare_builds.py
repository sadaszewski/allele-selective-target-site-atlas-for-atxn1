#!/usr/bin/env python3
"""Diff two results.json files, so a change of callset is auditable.

Usage: python compare_builds.py results_phase3.json results_30x.json
"""
import json, sys
from pathlib import Path

FIELDS = [
    ("variants scanned",             "n_records",                 "{:,}"),
    ("biallelic SNVs",               "n_biallelic_snv",           "{:,}"),
    ("retained after HWE filter",    "n_clean",                   "{:,}"),
    ("SpCas9 PAM-altering (window)", "pam_counts_window.SpCas9",  "{:,}"),
    ("any-nuclease PAM-altering",    "pam_counts_window.ANY",     "{:,}"),
    ("in-gene SpCas9, het>=20%",     "in_gene_sp_ge20",           "{:,}"),
    ("in-gene SpCas9, het>=30%",     "in_gene_sp_ge30",           "{:,}"),
    ("tier-1 sites",                 "tier1_n",                   "{:,}"),
    ("panel candidate sites",        "n_panel_candidates",        "{:,}"),
    ("HTT PAM-altering sites",       "comparison.0.n_sites",      "{:,}"),
    ("HTT phase-agnostic, EUR",      "comparison.0.pa_EUR",       "{:.2%}"),
    ("HTT phase-agnostic, all",      "comparison.0.pa_ALL",       "{:.2%}"),
    ("ATXN1 PAM-altering sites",     "comparison.1.n_sites",      "{:,}"),
    ("ATXN1 phase-agnostic, EUR",    "comparison.1.pa_EUR",       "{:.2%}"),
    ("ATXN1 phase-agnostic, all",    "comparison.1.pa_ALL",       "{:.2%}"),
    ("rs363099 het, EUR",            "rs363099.het_EUR",          "{:.2%}"),
    ("rs363099 single-site elig",    "rs363099.single_site_elig_EUR", "{:.2%}"),
    ("panel guide 1",                "panel.0.rsid",              "{}"),
    ("panel coverage, 1 guide",      "panel.0.cumulative",        "{:.2%}"),
    ("panel coverage, 6 guides",     "panel.5.cumulative",        "{:.2%}"),
    ("founder floor (26 pops)",      "founder_floor",             "{:.2%}"),
    ("rs2075974 GRCh38 pos",         "rs2075974.pos_grch38",      "{:,}"),
    ("rs2075974 PAM created",        "rs2075974.created",         "{}"),
    ("rs2075974 protospacer",        "rs2075974.protospacer",     "{}"),
    ("rs2075974 het, all",           "rs2075974.het_ALL",         "{:.2%}"),
    ("rs2075974 eligibility, EUR",   "rs2075974.eligibility.EUR", "{:.2%}"),
]


def get(d, path):
    try:
        for k in path.split("."):
            d = d[int(k)] if k.lstrip("-").isdigit() else d[k]
        return d
    except Exception:
        return None


def main(a, b):
    A, B = (json.loads(Path(x).read_text()) for x in (a, b))
    print(f"A = {A.get('build_label', a)}")
    print(f"B = {B.get('build_label', b)}\n")
    print(f"{'quantity':<32}{'A':>18}{'B':>18}   change")
    print("-" * 82)
    for label, path, fmt in FIELDS:
        x, y = get(A, path), get(B, path)
        if x is None or y is None:
            print(f"{label:<32}{'-':>18}{'-':>18}")
            continue
        fx, fy = fmt.format(x), fmt.format(y)
        if isinstance(x, (int, float)) and isinstance(y, (int, float)) and not isinstance(x, bool):
            d = y - x
            ch = "same" if d == 0 else (f"{d:+,}" if fmt == "{:,}" else f"{d:+.2%}"
                                        if fmt.endswith("%}") else f"{d:+}")
        else:
            ch = "same" if x == y else "DIFFERS"
        print(f"{label:<32}{fx:>18}{fy:>18}   {ch}")
    print("-" * 82)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(*sys.argv[1:])
