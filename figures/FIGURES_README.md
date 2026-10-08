# Figures — an allele-selective target-site atlas for ATXN1

Four figures as SVG, plus 600 dpi PNG rasters of each.

| File | Size | Shows |
| --- | --- | --- |
| Figure1_ATXN1_locus_candidates | 174 × 82 mm | The 468 in-gene SpCas9 candidate sites by position and heterozygosity across the locus |
| Figure2_PAM_vs_repeat_length | 174 × 108 mm | The two discrimination axes: PAM presence versus repeat length |
| Figure3_founder_floor_by_population | 174 × 141 mm | Worst-case founder eligibility, 26 populations |
| Figure4_guide_panel_coverage | 174 × 88 mm | Cumulative coverage of the greedy guide panel by superpopulation |

## Specification

- **Width 174 mm**, the usual double-column measure. Each file declares its physical
  size in millimetres with a matching viewBox, so it places at true size without
  rescaling. For a single-column journal, set the width to 85 mm and the type scales
  with it — the smallest label then lands at about 3.6 pt, which is too small, so a
  single-column version needs the type sizes raised rather than the figure shrunk.
- **Text is live text, not outlined paths**, so a copyeditor can correct it. Font is
  declared as `Helvetica, Arial, sans-serif`. If your journal requires embedded or
  outlined fonts at final submission, convert at that point.
- **Type sizes** at 174 mm: title 9.7 pt, labels 8.4 pt, small labels and axis text
  7.4 pt. All above the usual 6–7 pt floor.
- **Colour** is one blue (`#1F6FD0`) against greys. This is colour-blind safe and
  survives greyscale conversion, because the accent is also the only series drawn
  with a heavier stroke or a filled marker.
- **No external assets**, no embedded rasters, no CSS variables, white background
  rectangle included.

## Regenerating

```bash
python make_figures.py --results-dir results_phase3 --out figures
```

Figures 1, 3 and 4 are drawn directly from `results_phase3/`, so they follow the data:
re-run the analysis and the figures change with it. Figure 2 is schematic and carries
representative allele sizes (ATXN1 30/45, HTT 17/45) declared at the top of its
function; the two ratios printed on it are computed from those numbers rather than
typed, so they cannot disagree with the bars.

## A note on Figure 2

The allele sizes are representative, not measured: normal ATXN1 alleles cluster around
29–30 repeats and pathogenic ones begin at 39, while normal HTT is around 17 against 45
or more when expanded. The figure's point is the ratio between the two axes, and the
subtitle says the sizes are representative.
