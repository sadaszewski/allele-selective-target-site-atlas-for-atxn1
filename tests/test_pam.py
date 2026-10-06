"""Unit tests for the PAM logic. Run: python -m pytest tests/ -q"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from atxn1atlas.pam import revcomp, matches, pam_states, screen_variant, protospacer


def test_revcomp():
    assert revcomp("ACGT") == "ACGT"
    assert revcomp("AAG") == "CTT"


def test_iupac():
    assert matches("NGG", "AGG") and matches("NGG", "TGG")
    assert not matches("NGG", "AGA")
    assert matches("NNGRRT", "CCGAGT") and matches("NNGRRT", "CCGGAT")
    assert not matches("NNGRRT", "CCGCCT")
    assert matches("TTTV", "TTTA") and not matches("TTTV", "TTTT")


def test_rs363099_destroys_ngg():
    """HTT rs363099: ref C gives CCG (plus) = CGG on the minus strand; alt T does not."""
    ctx = "TGGAGGGTTTCT" + "C" + "CGCTCAGCCTTG"
    i = 12
    created, destroyed = screen_variant(ctx, i, "C", "T", "NGG")
    assert destroyed == 1 and created == 0


def test_rs2075974_creates_ngg():
    """ATXN1 rs2075974: alt C creates AGG on the minus strand; ref T does not."""
    ctx = "AGAAGGGGAGGC" + "T" + "TCACGATGAGTG"
    i = 12
    created, destroyed = screen_variant(ctx, i, "T", "C", "NGG")
    assert created == 1 and destroyed == 0
    # the rare third allele must NOT create a PAM (fails safe)
    assert screen_variant(ctx, i, "T", "G", "NGG") == (0, 0)


def test_rs2075974_protospacer():
    ctx = "AGAAGGGGAGGC" + "T" + "TCACGATGAGTG"
    i = 12
    created = pam_states(ctx, i, "C", "NGG") - pam_states(ctx, i, "T", "NGG")
    strand, start = sorted(created)[0]
    assert strand == "-"
    ps = protospacer(ctx, i, "C", start, strand, length=8)
    assert len(ps) == 8


def test_both_strands_scanned():
    """A variant creating GG only on the reverse strand must still be detected."""
    ctx = "AAAAAAAAAAAA" + "A" + "CCAAAAAAAAAA"
    i = 12
    created, _ = screen_variant(ctx, i, "A", "C", "NGG")
    assert created >= 1


def test_window_enumeration_count():
    """Every window overlapping the variant is tested: L placements for a motif of length L."""
    ctx = "G" * 25
    hits = pam_states(ctx, 12, "G", "NGG")
    # poly-G matches NGG at all 3 overlapping offsets on the plus strand only
    # (revcomp("GGG") == "CCC", which is not NGG)
    assert sorted(hits) == [("+", 10), ("+", 11), ("+", 12)]


def test_minus_strand_offsets():
    """A poly-C run is an NGG PAM on the minus strand at all 3 offsets, and none on the plus."""
    ctx = "C" * 25
    hits = pam_states(ctx, 12, "C", "NGG")
    assert sorted(hits) == [("-", 10), ("-", 11), ("-", 12)]
