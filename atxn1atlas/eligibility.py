"""Eligibility measures, the guide panel, and the founder model."""
import numpy as np
from . import config as C


def phase_agnostic(df, haps, created, destroyed, superpop):
    """Fraction of individuals for whom a selective guide exists whichever
    haplotype carries the expansion.

    An individual qualifies when BOTH haplotypes are independently targetable:
    some heterozygous site puts the PAM-bearing allele on haplotype 1, and some
    (not necessarily the same) puts it on haplotype 2.
    """
    h1, h2 = haps[:, :, 0], haps[:, :, 1]
    het = h1 != h2
    t1 = het & (((h1 == 1) & created[:, None]) | ((h1 == 0) & destroyed[:, None]))
    t2 = het & (((h2 == 1) & created[:, None]) | ((h2 == 0) & destroyed[:, None]))
    any_het = het.any(axis=0)
    both = t1.any(axis=0) & t2.any(axis=0)
    out = {"het_any_ALL": float(any_het.mean()), "phase_agnostic_ALL": float(both.mean())}
    for k in C.SUPERPOPS:
        m = superpop == k
        out[f"het_any_{k}"] = float(any_het[m].mean())
        out[f"phase_agnostic_{k}"] = float(both[m].mean())
    return out


def greedy_panel(df, haps, superpop, n_guides=12):
    """Greedy set cover: at each step add the site covering the most
    not-yet-covered individuals, where covered == heterozygous there."""
    het = haps[:, :, 0] != haps[:, :, 1]
    remaining = np.ones(het.shape[1], bool)
    steps = []
    for _ in range(n_guides):
        gains = (het & remaining).sum(axis=1)
        j = int(gains.argmax())
        if gains[j] == 0:
            break
        remaining &= ~het[j]
        cov = 1 - remaining.mean()
        by = {k: float(1 - remaining[superpop == k].mean()) for k in C.SUPERPOPS}
        steps.append({"rsid": df.iloc[j]["rsid"], "pos": int(df.iloc[j]["pos"]),
                      "cumulative": float(cov), **{f"cum_{k}": v for k, v in by.items()}})
    return steps


def founder_sweep(P, group_labels):
    """Worst-case eligibility under a single-founder model, per group.

    The disease chromosome is fixed to one complete haplotype across all sites;
    the wild-type chromosome is drawn from the same group. A wild-type haplotype
    fails only if it carries a PAM at EVERY site where the founder does.
    """
    out = []
    for g in sorted(set(group_labels)):
        m = group_labels == g
        Pk = P[:, m]
        n = Pk.shape[1]
        elig = np.zeros(n)
        for j in range(n):
            f = Pk[:, j]
            if not f.any():
                continue
            elig[j] = 1 - Pk[f, :].all(axis=0).mean()
        uniq = len(np.unique(Pk.T, axis=0))
        _, counts = np.unique(Pk.T, axis=0, return_counts=True)
        out.append({"group": g, "n_hap": n,
                    "worst": float(elig.min()), "p1": float(np.percentile(elig, 1)),
                    "median": float(np.percentile(elig, 50)), "best": float(elig.max()),
                    "distinct_profiles": uniq,
                    "commonest_profile_n": int(counts.max()),
                    "commonest_profile_frac": float(counts.max() / n)})
    return out
