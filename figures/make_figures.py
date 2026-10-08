#!/usr/bin/env python3
"""Generate print-ready SVG figures from the analysis outputs.

Usage: python make_figures.py [--results-dir results_phase3] [--out figures]

Emits one SVG per figure, sized for a 174 mm double-column journal page, with
literal colours (no CSS variables), real editable text, and no external assets.
"""
import argparse, json, html
from pathlib import Path
import pandas as pd

# --- print palette: blue + greys, colour-blind safe, legible in greyscale -----
INK, QUIET, GRID, AXIS, MUTED = "#1A1A1A", "#595959", "#DCDCDC", "#8C8C8C", "#B0B0B0"
ACCENT, ACCENT_FILL, REF = "#1F6FD0", "#D7E6F7", "#6E6E6E"
FONT = "Helvetica, Arial, sans-serif"
COL_MM = 174.0          # double-column width


def svg_open(h, title, desc):
    w_mm = COL_MM
    h_mm = COL_MM * h / 760.0
    return (f'<?xml version="1.0" encoding="UTF-8"?>\n'
            f'<svg xmlns="http://www.w3.org/2000/svg" version="1.1" '
            f'width="{w_mm:.2f}mm" height="{h_mm:.2f}mm" viewBox="0 0 760 {h}" '
            f'role="img" font-family="{FONT}">\n'
            f'<title>{html.escape(title)}</title>\n<desc>{html.escape(desc)}</desc>\n'
            f'<rect x="0" y="0" width="760" height="{h}" fill="#FFFFFF"/>\n')


def txt(x, y, s, size=13, fill=INK, anchor="start", weight="normal"):
    return (f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" fill="{fill}" '
            f'text-anchor="{anchor}" font-weight="{weight}">{html.escape(s)}</text>\n')


def line(x1, y1, x2, y2, stroke=AXIS, w=1.0, dash=None):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return (f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
            f'stroke="{stroke}" stroke-width="{w}"{d}/>\n')


def rect(x, y, w, h, fill="none", stroke=None, sw=1.0, rx=0):
    s = f' stroke="{stroke}" stroke-width="{sw}"' if stroke else ""
    return f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" rx="{rx}" fill="{fill}"{s}/>\n'


# ----------------------------------------------------------------- figure 1
def figure1(atlas, out):
    g = atlas[(atlas.region == "ATXN1_gene")
              & ((atlas.SpCas9_created + atlas.SpCas9_destroyed) > 0)
              & (atlas.hetf_ALL >= 0.05)].sort_values("pos_grch37")
    g0, g1 = 16_299_343, 16_761_691
    left, right, y0, base, H = 80, 716, 92, 300, 360
    X = lambda p: left + (p - g0) / (g1 - g0) * (right - left)
    Y = lambda h: base - h / 0.5 * (base - y0)
    tract, snp_pos = 16_327_910, 16_327_330
    snp = g[g.pos_grch37 == snp_pos].iloc[0]
    over = int((g.hetf_ALL >= 0.20).sum())

    s = svg_open(H, "Candidate sites span the whole 462 kb of ATXN1",
                 f"{len(g)} in-gene SpCas9 PAM-altering variants by position and "
                 f"observed heterozygosity; {over} are at least 20% heterozygous.")
    s += txt(24, 28, "Candidate sites span the whole 462 kb of ATXN1, not just the repeat", 15, INK, weight="bold")
    s += txt(24, 48, f"Observed heterozygosity, {len(g)} in-gene SpCas9 PAM-altering variants, "
                     "2,504 unrelated samples", 11.5, QUIET)
    for t in [0, .1, .2, .3, .4, .5]:
        s += line(left, Y(t), right, Y(t), GRID, 0.8)
        s += txt(left - 10, Y(t) + 4, "50%" if t == .5 else f"{t*100:.0f}", 11.5, QUIET, "end")
    s += line(left, y0, left, base, AXIS, 1.1)
    for t in range(16_300_000, 16_700_001, 100_000):
        s += txt(X(t), base + 20, f"{t/1e6:.2f}", 11.5, QUIET, "middle")
    s += line(left, base, right, base, AXIS, 1.1)
    s += txt((left + right) / 2, H - 12, "chromosome 6 position (Mb, GRCh37)", 11.5, QUIET, "middle")
    s += line(left, Y(.20), right, Y(.20), REF, 1.1, "4 3")
    s += txt(right, Y(.20) - 8, f"{over} of {len(g)} sites are at least 20% heterozygous", 11.5, QUIET, "end")
    s += line(X(tract), y0, X(tract), base, ACCENT, 1.1, "3 3")
    s += txt(X(tract) + 7, 84, "CAG repeat", 11.5, INK, weight="bold")
    s += '<g>\n'
    for _, r in g.iterrows():
        s += (f'<circle cx="{X(r.pos_grch37):.1f}" cy="{Y(r.hetf_ALL):.1f}" r="2.2" '
              f'fill="{MUTED}" fill-opacity="0.65"/>\n')
    s += '</g>\n'
    s += (f'<circle cx="{X(snp.pos_grch37):.1f}" cy="{Y(snp.hetf_ALL):.1f}" r="4.5" fill="{ACCENT}"/>\n')
    s += txt(X(snp.pos_grch37) + 12, Y(snp.hetf_ALL) + 4, "rs2075974", 11.5, INK, weight="bold")
    (out / "Figure1_ATXN1_locus_candidates.svg").write_text(s + "</svg>\n")
    return len(g), over


# ----------------------------------------------------------------- figure 2
def figure2(out):
    T = {("ATXN1", "expanded"): 45, ("ATXN1", "normal"): 30,
         ("HTT", "expanded"): 45, ("HTT", "normal"): 17}
    ppr, bar0, bar1, outX, H = 3.2, 290, 190, 596, 470
    rowE, rowN, snpX = 124, 202, 228
    s = svg_open(H, "A PAM-altering variant separates the two alleles where repeat length cannot",
                 "Panel A: at rs2075974 the expanded chromosome carries C, completing an NGG PAM, "
                 "and is cut; the normal chromosome carries T, has no PAM and cannot be cut. "
                 "Panel B: the expanded ATXN1 tract is 1.5x the normal one against 2.6x at HTT.")
    s += txt(24, 28, "A PAM-altering variant separates the two alleles where repeat length cannot", 15, INK, weight="bold")
    s += txt(24, 48, "Representative allele sizes; rs2075974 lies 537 bp from the ATXN1 CAG tract", 11.5, QUIET)
    s += txt(24, 82, "A   Cas9 reads a single base at rs2075974", 13, INK, weight="bold")
    for y, base_, on, pam, outw in ((rowE, "C", True, "AGG = PAM", "Cas9 cuts"),
                                    (rowN, "T", False, "AAG = no PAM", "cannot be cut")):
        n = T[("ATXN1", "expanded" if on else "normal")]
        s += rect(160, y - 2, 410, 4, GRID, rx=2)
        s += rect(bar0, y - 9, n * ppr, 18, ACCENT_FILL if on else "#F2F2F0",
                  ACCENT if on else AXIS, 1.1, 4)
        s += rect(snpX - 15, y - 13, 30, 26, ACCENT_FILL if on else "#FFFFFF",
                  ACCENT if on else AXIS, 2 if on else 1.1, 6)
        s += txt(snpX, y + 5, base_, 14, INK, "middle", "bold")
        s += txt(24, y + 5, "Expanded allele" if on else "Normal allele", 13, INK)
        s += txt(snpX, y - 24, pam, 11.5, INK if on else QUIET, "middle", "bold" if on else "normal")
        s += txt(bar0 + n * ppr / 2, y - 15, f"{n} CAG", 11.5, QUIET, "middle")
        s += txt(outX, y + 5, outw, 13, INK if on else QUIET, weight="bold" if on else "normal")
    s += line(24, 258, 736, 258, GRID, 0.8)
    s += txt(24, 288, "B   Repeat length gives too little to work with at ATXN1", 13, INK, weight="bold")
    rows = [("ATXN1 expanded", T[("ATXN1", "expanded")], 322, True),
            ("ATXN1 normal", T[("ATXN1", "normal")], 354, False),
            ("HTT expanded", T[("HTT", "expanded")], 398, False),
            ("HTT normal", T[("HTT", "normal")], 430, False)]
    for lab, n, y, on in rows:
        s += txt(24, y + 5, lab, 11.5, INK)
        s += rect(bar1, y - 8, n * ppr, 16, ACCENT_FILL if on else "#F2F2F0",
                  ACCENT if on else AXIS, 1.1, 4)
        s += txt(bar1 + n * ppr + 10, y + 5, str(n), 11.5, QUIET)
    for y1_, y2_, ratio, strong in ((314, 362, T[("ATXN1", "expanded")] / T[("ATXN1", "normal")], True),
                                    (390, 438, T[("HTT", "expanded")] / T[("HTT", "normal")], False)):
        s += line(520, y1_, 520, y2_, ACCENT if strong else AXIS, 1.1)
        s += txt(532, (y1_ + y2_) / 2 + 4, f"{ratio:.1f}× longer", 11.5,
                 INK if strong else QUIET, weight="bold" if strong else "normal")
    (out / "Figure2_PAM_vs_repeat_length.svg").write_text(s + "</svg>\n")


# ----------------------------------------------------------------- figure 3
def figure3(founder, out):
    f = founder.sort_values("worst").reset_index(drop=True)
    lo, hi, left, right, y0, step = 95, 100, 150, 700, 118, 17
    H = y0 + len(f) * step + 56
    # the axis runs 95-100 in percent; the table stores fractions
    X = lambda frac: left + (frac * 100 - lo) / (hi - lo) * (right - left)
    floor = float(f.worst.min())
    s = svg_open(H, "Worst-case founder eligibility never falls below 96.5%",
                 "Eligibility under an adversarial single-founder model, one point per "
                 "1000 Genomes population; the floor is 96.46% in LWK.")
    s += txt(24, 28, "Worst-case founder eligibility never falls below 96.5%", 15, INK, weight="bold")
    s += txt(24, 50, "Disease chromosome fixed to one haplotype; wild-type chromosome drawn from the same population", 11.5, QUIET)
    s += txt(24, 70, "Eligible individuals (%), each of 26 populations, worst case over all candidate founders", 11.5, QUIET)
    for t in range(lo, hi + 1):
        s += line(X(t / 100), y0 - 10, X(t / 100), y0 + len(f) * step, GRID, 0.8)
        s += txt(X(t / 100), y0 + len(f) * step + 20, str(t), 11.5, QUIET, "middle")
    s += txt((left + right) / 2, y0 + len(f) * step + 42,
             "eligible individuals (%)", 11.5, QUIET, "middle")
    s += line(X(floor), y0 - 16, X(floor), y0 + len(f) * step, REF, 1.1, "4 3")
    s += txt(X(floor), y0 - 24, f"floor {floor*100:.2f}%", 11.5, QUIET, "middle")
    for i, r in f.iterrows():
        yy = y0 + i * step + 4
        is_min = float(r.worst) == floor
        s += txt(24, yy + 4, str(r["super"]), 11.5, QUIET)
        s += txt(74, yy + 4, str(r["pop"]), 11.5, INK)
        s += line(left, yy, X(r.worst), yy, GRID, 0.8)
        s += (f'<circle cx="{X(r.worst):.1f}" cy="{yy:.1f}" r="4.5" '
              f'fill="{ACCENT if is_min else MUTED}"/>\n')
        if is_min:
            s += txt(X(r.worst) + 12, yy + 4, f"{r.worst*100:.2f}%", 11.5, ACCENT, weight="bold")
    (out / "Figure3_founder_floor_by_population.svg").write_text(s + "</svg>\n")


# ----------------------------------------------------------------- figure 4
def figure4(panel, out):
    keys = ["AFR", "AMR", "EAS", "EUR", "SAS"]
    left, right, y0, base, lo, H = 92, 616, 100, 328, 0.40, 384
    n_max = len(panel)
    X = lambda n: left + (n - 1) / (n_max - 1) * (right - left)
    Y = lambda v: base - (v - lo) / (1 - lo) * (base - y0)
    six = panel[5]
    s = svg_open(H, "Six guides cover at least 96.8% of individuals in every ancestry group",
                 "Cumulative share of individuals heterozygous at one or more panel sites, "
                 "greedy cover over 468 candidate sites, by superpopulation.")
    s += txt(24, 28, "Six guides cover at least 96.8% of individuals in every ancestry group", 15, INK, weight="bold")
    s += txt(24, 48, "Cumulative share heterozygous at one or more panel sites, greedy cover over 468 candidates", 11.5, QUIET)
    for t in [.4, .5, .6, .7, .8, .9, 1.0]:
        s += line(left, Y(t), right, Y(t), GRID, 0.8)
        s += txt(left - 10, Y(t) + 4, "100%" if t == 1 else f"{t*100:.0f}", 11.5, QUIET, "end")
    s += line(left, y0, left, base, AXIS, 1.1)
    for r in panel:
        s += txt(X(r["n"]), base + 20, str(r["n"]), 11.5, INK if r["n"] == 6 else QUIET,
                 "middle", "bold" if r["n"] == 6 else "normal")
    s += line(left, base, right, base, AXIS, 1.1)
    s += txt((left + right) / 2, H - 12, "number of guides in the panel", 11.5, QUIET, "middle")
    s += line(X(6), y0, X(6), base, REF, 1.1, "4 3")
    s += txt(X(6), 92, "6 guides", 11.5, INK, "middle", "bold")
    for k in keys:
        if k == "AFR":
            continue
        d = " ".join(f'{"M" if i==0 else "L"}{X(r["n"]):.1f} {Y(r[k]):.1f}' for i, r in enumerate(panel))
        s += f'<path d="{d}" fill="none" stroke="{MUTED}" stroke-width="1.4"/>\n'
    d = " ".join(f'{"M" if i==0 else "L"}{X(r["n"]):.1f} {Y(r["AFR"]):.1f}' for i, r in enumerate(panel))
    s += f'<path d="{d}" fill="none" stroke="{ACCENT}" stroke-width="2.4"/>\n'
    for r in panel:
        s += f'<circle cx="{X(r["n"]):.1f}" cy="{Y(r["AFR"]):.1f}" r="3" fill="{ACCENT}"/>\n'
    s += txt(right + 10, Y(panel[-1]["AFR"]) + 4, "AFR, the lowest", 11.5, INK, weight="bold")
    s += txt(right + 10, Y(panel[-1]["AFR"]) + 22, "AMR, EAS, EUR, SAS", 11.5, QUIET)
    s += (f'<circle cx="{X(6):.1f}" cy="{Y(six["AFR"]):.1f}" r="5.5" fill="none" '
          f'stroke="{ACCENT}" stroke-width="2"/>\n')
    s += txt(X(6) + 14, Y(six["AFR"]) + 28, f'{six["AFR"]*100:.1f}%', 11.5, INK, weight="bold")
    (out / "Figure4_guide_panel_coverage.svg").write_text(s + "</svg>\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", default="results_phase3")
    ap.add_argument("--out", default="figures")
    a = ap.parse_args()
    rd, out = Path(a.results_dir), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    atlas = pd.read_csv(rd / "ATXN1_allele_selective_atlas_v0.csv.gz")
    founder = pd.read_csv(rd / "founder_worstcase_by_population.csv")
    R = json.loads((rd / "results.json").read_text())
    panel = [{"n": i, "AFR": s["cum_AFR"], "AMR": s["cum_AMR"], "EAS": s["cum_EAS"],
              "EUR": s["cum_EUR"], "SAS": s["cum_SAS"]}
             for i, s in enumerate(R["panel"], 1)]
    n, over = figure1(atlas, out)
    figure2(out)
    figure3(founder, out)
    figure4(panel, out)
    print(f"Figure 1: {n} sites, {over} at >=20%")
    for f in sorted(out.glob("*.svg")):
        print(f"  {f.stat().st_size:>8,} B  {f.name}")


if __name__ == "__main__":
    main()
