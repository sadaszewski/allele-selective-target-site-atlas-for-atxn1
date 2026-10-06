# ATXN1 allele-selective target-site atlas

Reproduces every number in *An allele-selective target-site atlas for ATXN1*.

One command, no credentials, no data-use agreements. All inputs are public.

```bash
pip install -r requirements.txt
python run_all.py                      # 30x GRCh38 callset (default)
python run_all.py --build phase3       # the GRCh37 callset the manuscript used
python verify.py                       # check the output against expected values
```

`verify.py` exits 0 when every check passes. For `--build phase3` that is 37 checks
against the values reported in the manuscript: if it passes on your machine, you have
reproduced the paper.

## Two callsets

| `--build` | Callset | Assembly | Download |
| --- | --- | --- | --- |
| `grch38_30x` *(default)* | 1000G 30x NYGC, shapeit2-duohmm phased | GRCh38 | ~4.6 GB |
| `phase3` | 1000G Phase 3, low coverage (4–6x) | GRCh37 | ~2.1 GB |

Both are phased and both are restricted to the same 2,504 unrelated samples, so the
two runs are directly comparable. The 30x callset is higher coverage and
trio-aware-phased; Phase 3 is what the manuscript analysed. The run aborts if the
sample intersection is not exactly 2,504.

```bash
python compare_builds.py results_phase3/results.json results_30x/results.json
```

prints a side-by-side diff of every headline number, so a change of callset is
auditable rather than a matter of trust.

### What the 30x re-run changed

Nothing that matters. The rare-variant tail shrinks — 24,706 biallelic SNVs become
18,819, because the 30x callset is more stringently filtered — but the common-variant
core that the conclusions rest on is unmoved:

| | Phase 3 | 30x |
| --- | --- | --- |
| in-gene SpCas9 sites, het ≥20% | 257 | 258 |
| in-gene SpCas9 sites, het ≥30% | 170 | 170 |
| tier-1 sites | 85 | 86 |
| *ATXN1* phase-agnostic eligibility, EUR | 99.80% | 99.80% |
| *HTT* phase-agnostic eligibility, EUR | 78.13% | 79.13% |
| founder floor, 26 populations | 96.46% | 96.46% |
| rs2075974 eligibility, EUR | 78.53% | 78.53% |

Two incidental confirmations. The HWE filter drops 5 sites from Phase 3 and 0 from
the 30x callset — the mapping artefacts it is designed to catch are already absent
from the higher-coverage data. And the CAG tract, located from sequence alone in each
build, lands at GRCh37 16,327,867–16,327,953 and GRCh38 16,327,636–16,327,722: a
difference of exactly the 231 bp offset.

### rsIDs in the 30x callset

The 30x NYGC VCF's ID column holds `chrom:pos:ref:alt` strings rather than rsIDs.
Named variants are therefore located by position, not by ID (`config.KNOWN_VARIANTS`,
verified against dbSNP in both assemblies and checked against the reference allele at
load). To put rsIDs back on a 30x atlas, join it against a Phase 3 atlas:

```bash
python annotate_rsids.py results_30x/ATXN1_allele_selective_atlas_v0.csv.gz \
                         results_phase3/ATXN1_allele_selective_atlas_v0.csv.gz
```

The join is exact on (pos_grch37, ref, alt) and needs no network. 86% of rows match;
the remainder are variants Phase 3 did not call.

## What it does

| Stage | Output |
| --- | --- |
| 1. Atlas | Scans *ATXN1* ±200 kb for biallelic SNVs, filters mapping artefacts, screens every variant against four nucleases on both strands |
| 2. Validation | Runs the identical code over *HTT*, recovers rs363099, reproduces its published ~20% eligibility |
| 3. Panel | Greedy set cover over candidate sites |
| 4. Founder model | Adversarial single-founder sweep over all 5,008 haplotypes, by superpopulation and by population |
| 5. rs2075974 | PAM direction, protospacer, per-ancestry eligibility |

## Outputs (`results/`)

- `ATXN1_allele_selective_atlas_v0.csv.gz` — all 24,701 variants: both builds, per-ancestry heterozygosity, per-nuclease PAM annotations
- `ATXN1_tier1_candidates.csv` — the 85 sites ≥30% heterozygous in all five superpopulations
- `founder_worstcase_by_population.csv` — founder sweep, 26 populations
- `results.json` — every headline number, machine-readable

## Data sources

Downloaded automatically and cached in `data/`:

- 1000 Genomes Phase 3 chr6 and chr4 phased callsets, and the 2,504-sample panel (AWS S3 open data, ~2.1 GB)
- GRCh37 reference sequence for both windows (NCBI E-utilities, ~1.5 MB)

`--skip-htt` omits the *HTT* control and saves a large download. It also skips the
validation step, so do not use it for a first run.

## Layout

```
atxn1atlas/
  config.py       coordinates, URLs, filter thresholds — the only file to edit
  pam.py          PAM motif logic (the core; ~30 lines, fully unit-tested)
  data.py         download, cache, load; locates the CAG tract from sequence alone
  atlas.py        window scan, HWE filter, PAM profile matrix
  eligibility.py  phase-agnostic eligibility, greedy panel, founder model
  report.py       table printing
run_all.py        pipeline entry point       (--build, --skip-htt, --results-dir)
verify.py         checks results.json against expected/published_values.json
compare_builds.py side-by-side diff of two runs
annotate_rsids.py recover rsIDs for a 30x atlas by joining a Phase 3 atlas
tests/            unit tests for the PAM logic  (python -m pytest tests/ -q)
```

## Notes on the method

**Positions are reported in both builds.** The analysed callset supplies real
positions; the other build is derived from the gene's constant inter-build offset
(ATXN1 −231 bp, HTT −1,727 bp). `config.locus_offset` asserts at import that the
offset is identical at both gene termini, so a mistyped coordinate fails immediately
rather than silently shifting the atlas.

**The reference check is the guard against coordinate error.** Every retained REF
allele is compared against independently fetched reference sequence; a mismatch
aborts the run. Expect 24,701/24,701.

**Phase-agnostic eligibility** counts an individual only when *both* haplotypes are
independently targetable, so a selective guide exists whichever one carries the
expansion. It assumes nothing about phase and therefore cannot overstate eligibility
on that account.

**The founder model is deliberately adversarial.** It fixes the disease chromosome
to one complete haplotype across all 468 sites and draws the wild-type chromosome
from the same single population — strictly more constraining than any real founder
effect, which decays with distance from the repeat.

## A note on mirrors

The 30x phased release is fetched from EBI FTP. The AWS S3 mirror of that same
release contains zero-byte placeholder objects — it will appear to download fine and
yield an empty file. Do not point `BUILDS["grch38_30x"]["vcf_url"]` at S3.

## Runtime

~80 s of compute after the download, single core, ~2 GB RAM. Tested on Python 3.13
with pysam 0.24.1, NumPy 2.5.3, pandas 3.0.5.
