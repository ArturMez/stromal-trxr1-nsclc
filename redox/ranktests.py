"""
The rank tests of the patient-level driver analyses.

WHAT: the two-sided Mann-Whitney test, the van Elteren test stratified by histology (or by cohort
      and histology) with a permutation p, and the stratified difference in medians with a
      percentile bootstrap interval.
FROM: plain arrays of per-patient values.
TO:   floats and dictionaries.
WHY:  used by tumor_vs_stroma.py (mw_test: the Mann-Whitney tests of Supplementary Figure 1d and the
      smoking comparison) and by chapter2/07_smoking_atlas.py (van_elteren and strat_med_diff:
      Supplementary Figure 2b). The permutations and the bootstrap draw from a generator that the
      caller passes in and consume it in a fixed order, so a caller that wants a repeatable result
      seeds the generator and keeps the order of its calls.

mw_test         two-sided Mann-Whitney U (scipy, method "auto"), group 1 against group 0
van_elteren     stratified Wilcoxon (van Elteren weights 1/(N_s + 1)); asymptotic Z and p, and a
                permutation p (labels shuffled within stratum, group sizes fixed)
strat_med_diff  weighted mean of the per-stratum differences in medians, weights n1 n0 / N,
                bootstrap within stratum and group
"""
import numpy as np
from scipy import stats


def mw_test(x1, x0):
    """Mann-Whitney U statistic and two-sided p (scipy, method "auto"), or two missing values when a
    group has fewer than two values. Group 1 is x1, group 0 is x0."""
    x1 = np.asarray(x1, float)
    x0 = np.asarray(x0, float)
    if len(x1) < 2 or len(x0) < 2:
        return np.nan, np.nan
    r = stats.mannwhitneyu(x1, x0, alternative="two-sided", method="auto")
    return float(r.statistic), float(r.pvalue)


def _ve_parts(strata):
    """strata: list of (x1, x0). Per stratum: ranks, n1, weight, expectation and variance."""
    parts = []
    for x1, x0 in strata:
        x1 = np.asarray(x1, float)
        x0 = np.asarray(x0, float)
        n1, n0 = len(x1), len(x0)
        if n1 == 0 or n0 == 0:
            continue
        v = np.concatenate([x1, x0])
        N = n1 + n0
        rk = stats.rankdata(v)
        _, cnt = np.unique(v, return_counts=True)
        tie = (cnt ** 3 - cnt).sum()
        var = n1 * n0 / 12.0 * ((N + 1) - tie / (N * (N - 1))) if N > 1 else 0.0
        parts.append(dict(rk=rk, n1=n1, N=N, w=1.0 / (N + 1), E=n1 * (N + 1) / 2.0, V=var))
    return parts


def van_elteren(strata, n_perm, rng):
    """Van Elteren stratified rank test: Z, asymptotic p and permutation p.

    strata is a list of (group 1 values, group 0 values), one pair per stratum. The permutation p
    shuffles the group labels within each stratum with group sizes fixed, n_perm times, drawing from
    the generator rng, and counts shuffles with |Z| at least as large as observed (plus one in
    numerator and denominator).
    """
    parts = _ve_parts(strata)
    if not parts:
        return dict(Z=np.nan, p_asym=np.nan, p_perm=np.nan, n_strata=0)
    T = sum(p["w"] * p["rk"][: p["n1"]].sum() for p in parts)
    ET = sum(p["w"] * p["E"] for p in parts)
    VT = sum(p["w"] ** 2 * p["V"] for p in parts)
    if VT <= 0:
        return dict(Z=np.nan, p_asym=np.nan, p_perm=np.nan, n_strata=len(parts))
    Z = (T - ET) / np.sqrt(VT)
    p_asym = 2 * stats.norm.sf(abs(Z))
    Tp = np.zeros(n_perm)
    for p in parts:
        perm = np.argsort(rng.random((n_perm, p["N"])), axis=1)[:, : p["n1"]]
        Tp += p["w"] * p["rk"][perm].sum(axis=1)
    Zp = (Tp - ET) / np.sqrt(VT)
    p_perm = (1 + np.sum(np.abs(Zp) >= abs(Z) - 1e-12)) / (n_perm + 1)
    return dict(Z=float(Z), p_asym=float(p_asym), p_perm=float(p_perm), n_strata=len(parts))


def strat_med_diff(strata, n_boot, rng, alpha=0.05):
    """Weighted mean of the per-stratum differences in medians, with a bootstrap interval.

    Weights are n1 n0 / N per stratum. The interval is the percentile interval (alpha 0.05, so 2.5th
    and 97.5th percentiles) of n_boot bootstrap resamples, each group resampled with replacement
    within its stratum, drawing from the generator rng. Returns estimate, lower, upper.
    """
    use = [(np.asarray(a, float), np.asarray(b, float)) for a, b in strata if len(a) > 0 and len(b) > 0]
    if not use:
        return np.nan, np.nan, np.nan
    w = np.array([len(a) * len(b) / (len(a) + len(b)) for a, b in use])
    w = w / w.sum()
    est = float(sum(wi * (np.median(a) - np.median(b)) for wi, (a, b) in zip(w, use)))
    boots = np.zeros(n_boot)
    for wi, (a, b) in zip(w, use):
        ba = np.median(rng.choice(a, (n_boot, len(a)), replace=True), axis=1)
        bb = np.median(rng.choice(b, (n_boot, len(b)), replace=True), axis=1)
        boots += wi * (ba - bb)
    lo, hi = np.percentile(boots, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return est, float(lo), float(hi)
