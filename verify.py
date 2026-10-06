#!/usr/bin/env python3
"""Compare results/results.json against the expected values for its callset.

Usage: python verify.py [results/results.json]
Exits 0 when every check passes.
"""
import json, sys
from pathlib import Path

EXPECTED = json.loads(Path(__file__).with_name("expected")
                     .joinpath("published_values.json").read_text())


def get(d, path):
    for k in path.split("."):
        d = d[int(k)] if k.lstrip("-").isdigit() else d[k]
    return d


def main(results="results/results.json"):
    R = json.loads(Path(results).read_text())
    build = R.get("build", "phase3")
    section = EXPECTED["builds"].get(build)
    print(f"callset : {build}  ({R.get('build_label', '?')})")
    if not section or not section.get("checks"):
        print(f"\nNo expected values recorded for '{build}' yet.")
        print("Nothing to verify against; inspect results.json and record them.")
        return 0
    checks = section["checks"]
    print(f"checks  : {len(checks)}\n")
    print(f"{'check':<44}{'expected':>22} {'got':>22}   status")
    print("-" * 94)
    ok = bad = 0
    for path, exp, tol in checks:
        try:
            got = get(R, path)
        except Exception:
            print(f"{path:<44}{str(exp):>22} {'MISSING':>22}   FAIL")
            bad += 1
            continue
        if isinstance(exp, (int, float)) and isinstance(got, (int, float)) and not isinstance(exp, bool):
            good = abs(got - exp) <= tol
        else:
            good = str(got) == str(exp)
        f = lambda v: f"{v:.4f}" if isinstance(v, float) else f"{v:,}" if isinstance(v, int) else str(v)
        print(f"{path:<44}{f(exp):>22} {f(got):>22}   {'ok' if good else 'FAIL'}")
        ok, bad = (ok + 1, bad) if good else (ok, bad + 1)
    print("-" * 94)
    print(f"{ok} passed, {bad} failed")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:]))
