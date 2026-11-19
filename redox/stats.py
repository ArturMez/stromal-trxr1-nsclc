"""
The four statistics the chapters use outside the Cox models.

WHAT: Spearman correlation, the two paired and unpaired rank tests, and a bootstrap confidence
      interval for a median.
FROM: plain arrays of per-patient or per-core values.
TO:   floats.
WHY:  every panel that draws one dot per patient carries a median with a 95 % interval, so the
      bootstrap is written once. The seed is a required argument, never a default, so that a caller
      has to name it and the README can list it. Used by chapter1.py, chapter3.py and chapter4.py.
"""
import numpy as np
from scipy import stats as _sp


def spearman(x, y):
    """Spearman rank correlation of two arrays, pairwise complete. Returns rho and its p value."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    r = _sp.spearmanr(x[ok], y[ok])
    return float(r.statistic), float(r.pvalue)


def wilcoxon_p(v):
    """Two-sided Wilcoxon signed-rank p for a set of paired differences against zero.

    Missing when every value is exactly zero, which the test cannot handle and which means there is
    nothing to test.
    """
    v = np.asarray(v, dtype=float)
    v = v[np.isfinite(v)]
    if not len(v) or not np.any(v != 0):
        return float("nan")
    return float(_sp.wilcoxon(v).pvalue)


def mann_whitney_p(a, b):
    """Two-sided Mann-Whitney p for two unpaired groups."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if not len(a) or not len(b):
        return float("nan")
    return float(_sp.mannwhitneyu(a, b, alternative="two-sided").pvalue)


def bootstrap_median_ci(v, seed, reps=2000):
    """Median of v with a 95 % percentile bootstrap interval. Returns median, lower, upper.

    Percentile bootstrap: reps resamples (default 2000) of the full length with replacement, drawn
    from numpy's default_rng(seed), the median of each, then the 2.5th and 97.5th percentiles.
    """
    v = np.asarray(v, dtype=float)
    v = v[np.isfinite(v)]
    if not len(v):
        return float("nan"), float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    draws = np.median(rng.choice(v, (int(reps), len(v)), replace=True), axis=1)
    lo, hi = np.percentile(draws, [2.5, 97.5])
    return float(np.median(v)), float(lo), float(hi)
