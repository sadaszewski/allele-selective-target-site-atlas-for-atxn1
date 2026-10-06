"""Coordinates, URLs and constants. All coordinates are 1-based inclusive."""
from pathlib import Path

DATA_DIR    = Path("data")
RESULTS_DIR = Path("results")

# --- Loci (GRCh37 / NCBI NC_000006.11, NC_000004.11) -------------------------
# Verified against NCBI Gene: ATXN1 GeneID 6310, HTT GeneID 3064.
ATXN1 = dict(
    name="ATXN1", gene_id=6310, chrom="6", acc37="NC_000006.11",
    start37=16_299_343, end37=16_761_691,
    start38=16_299_112, end38=16_761_460,   # offset = -231 bp, constant at both termini
    strand="-",
)
HTT = dict(
    name="HTT", gene_id=3064, chrom="4", acc37="NC_000004.11",
    start37=3_076_408, end37=3_245_687,
    start38=3_074_681, end38=3_243_960,
    strand="+",
)

FLANK = 200_000                       # atlas window = gene +/- FLANK
# CAG tract located directly in the reference (see locate_repeat); midpoint used
# only for distance reporting.
ATXN1_REPEAT_GRCH37 = (16_327_867, 16_327_953)   # located from sequence; see locate_repeat

SUPERPOPS = ["AFR", "AMR", "EAS", "EUR", "SAS"]

# Variants named in the literature. Anchored by POSITION, not by the VCF ID
# column: Phase 3 carries rsIDs there, the 30x NYGC callset carries
# "chrom:pos:ref:alt" strings. Positions below were confirmed against dbSNP in
# both assemblies and are consistent with each locus's constant offset.
KNOWN_VARIANTS = {
    "rs2075974": dict(locus="ATXN1", chrom="6", pos37=16_327_330, pos38=16_327_099,
                      ref="T", alt="C",
                      note="Mittal 2005: ancestral C allele associated with expanded chromosomes"),
    "rs363099":  dict(locus="HTT", chrom="4", pos37=3_162_056, pos38=3_160_329,
                      ref="C", alt="T",
                      note="Shin 2022: destroys the NGG PAM on the commonest normal HTT haplotype"),
}


def find_variant(df, name, assembly):
    """Locate a named variant in a scanned table by position; verify REF.

    Returns the matching row, or None. Raises if the position is present but
    the reference allele disagrees, which would mean a coordinate error.
    """
    v = KNOWN_VARIANTS[name]
    pos = v["pos38"] if assembly == "GRCh38" else v["pos37"]
    hit = df[df["pos"] == pos]
    if len(hit) == 0:
        return None
    row = hit.iloc[0]
    if row["ref"] != v["ref"]:
        raise RuntimeError(f"{name} at {assembly} pos {pos}: REF is {row['ref']}, "
                           f"expected {v['ref']} - coordinate error")
    return row

# --- Filters -----------------------------------------------------------------
MAX_HET          = 0.55   # above the Hardy-Weinberg ceiling of 0.5 -> mapping artefact
MAX_HWE_Z        = 8.0    # excess-heterozygosity z-score cutoff
PANEL_MIN_HET    = 0.05   # minimum heterozygosity for greedy-panel candidates
BENCH_MIN_HET    = 0.02   # minimum heterozygosity for the HTT/ATXN1 comparison

# --- Data sources (all anonymously downloadable) -----------------------------
S3  = "https://s3.amazonaws.com/1000genomes"
EBI = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

PANEL_URL = f"{S3}/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"

# Two interchangeable callsets. Both are phased and both are restricted to the
# same 2,504 unrelated samples, so results are directly comparable.
BUILDS = {
    "phase3": {
        "label": "1000G Phase 3, low coverage (4-6x), GRCh37",
        "assembly": "GRCh37",
        "vcf_url": S3 + "/release/20130502/ALL.chr{c}.phase3_shapeit2"
                        "_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",
        "contig": "{c}",
        "acc": {"6": "NC_000006.11", "4": "NC_000004.11"},
        "coord_keys": ("start37", "end37"),
        "n_samples_in_vcf": 2504,
    },
    "grch38_30x": {
        "label": "1000G 30x NYGC, shapeit2-duohmm phased, GRCh38",
        "assembly": "GRCh38",
        "vcf_url": EBI + "/data_collections/1000G_2504_high_coverage/working"
                         "/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05"
                         "_chr{c}.filtered.shapeit2-duohmm-phased.vcf.gz",
        "contig": "chr{c}",
        "acc": {"6": "NC_000006.12", "4": "NC_000004.12"},
        "coord_keys": ("start38", "end38"),
        "n_samples_in_vcf": 3202,   # 2,504 unrelated + 698 related; we keep the 2,504
    },
}
DEFAULT_BUILD = "grch38_30x"


def locus_offset(locus):
    """GRCh37 position minus GRCh38 position, constant across the gene.

    Verified at both termini for both loci; asserted at load time below.
    """
    a = locus["start37"] - locus["start38"]
    b = locus["end37"] - locus["end38"]
    assert a == b, f"{locus['name']}: offset differs at the two termini ({a} vs {b})"
    return a


for _loc in (ATXN1, HTT):          # fail fast if a coordinate is ever mistyped
    locus_offset(_loc)
