"""PAM motif logic. The whole screen rests on these ~30 lines."""

COMPLEMENT = str.maketrans("ACGT", "TGCA")

IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T",
         "N": "ACGT", "R": "AG", "Y": "CT", "V": "ACG"}

# (name, motif, PAM is 3' of the protospacer, protospacer length)
NUCLEASES = [
    ("SpCas9",    "NGG",    True,  20),
    ("SaCas9",    "NNGRRT", True,  21),
    ("xCas9_NG",  "NG",     True,  20),
    ("Cas12a",    "TTTV",   False, 23),
]


def revcomp(s: str) -> str:
    return s.translate(COMPLEMENT)[::-1]


def matches(motif: str, seq: str) -> bool:
    return len(seq) == len(motif) and all(c in IUPAC[m] for m, c in zip(motif, seq))


def pam_states(seq: str, i: int, base: str, motif: str) -> set:
    """(strand, start) pairs where `motif` matches, with seq[i] replaced by `base`.

    Enumerates every window overlapping position i, on both strands.
    """
    L = len(motif)
    out = set()
    for s in range(i - L + 1, i + 1):
        if s < 0 or s + L > len(seq):
            continue
        w = list(seq[s:s + L])
        w[i - s] = base
        w = "".join(w)
        if matches(motif, w):
            out.add(("+", s))
        if matches(motif, revcomp(w)):
            out.add(("-", s))
    return out


def screen_variant(seq: str, i: int, ref: str, alt: str, motif: str):
    """-> (n_created, n_destroyed) for one variant and one motif."""
    r = pam_states(seq, i, ref, motif)
    a = pam_states(seq, i, alt, motif)
    return len(a - r), len(r - a)


def protospacer(seq: str, i: int, alt: str, pam_start: int, strand: str, length: int = 20) -> str:
    """The protospacer 5' of an NGG PAM, on the given strand, with alt substituted."""
    w = list(seq)
    w[i] = alt
    w = "".join(w)
    if strand == "+":
        return w[pam_start - length:pam_start]
    return revcomp(w[pam_start + 3:pam_start + 3 + length])
