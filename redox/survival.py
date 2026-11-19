"""
The Cox models of the manuscript.

WHAT: the Cox solver (Newton-Raphson with Efron ties), the per-standard-deviation hazard ratio,
      the fixed-effect pooling of two cohorts, and the histology interaction.
FROM: a per-patient DataFrame carrying the exposure, OS_weeks, event and the covariates.
TO:   plain dictionaries of hazard ratio, confidence limits, p value, log hazard ratio and its
      standard error.
WHY:  every hazard ratio of the paper (Figure 1b and 1c, Figure 5g and the supplementary tables)
      comes from this one solver. Users: chapter1.py, chapter4.py and tumor_vs_stroma.py.

The solver also accumulates the log-likelihood in the same loop, so that the interaction can
report a likelihood-ratio check beside its Wald p, and falls back to least squares for a
singular Hessian. Every p value here is a Wald p from the normal distribution; nothing is
permuted or resampled.
"""
import math
import warnings

import numpy as np
from scipy.stats import chi2, norm


def _fit(X, t, d, max_iter=100, tol=1e-9):
    """Newton-Raphson Cox with Efron ties. Returns beta, standard errors and the log-likelihood."""
    X = np.atleast_2d(np.asarray(X, float))
    if X.shape[0] != len(t):
        X = X.T
    t = np.asarray(t, float)
    d = np.asarray(d, int)
    n, p = X.shape
    o = np.argsort(t, kind="stable")
    X, t, d = X[o], t[o], d[o]
    et = np.unique(t[d == 1])
    beta = np.zeros(p)

    def ghl(b):
        """Gradient, Hessian and log partial likelihood (Efron ties) at coefficients b."""
        eta = X @ b
        ex = np.exp(eta)
        g = np.zeros(p)
        H = np.zeros((p, p))
        ll = 0.0
        for T in et:
            risk = t >= T
            die = (t == T) & (d == 1)
            m = int(die.sum())
            if m == 0:
                continue
            rX = X[risk]
            rex = ex[risk]
            S0r = rex.sum()
            S1r = (rex[:, None] * rX).sum(0)
            S2r = (rex[:, None, None] * rX[:, :, None] * rX[:, None, :]).sum(0)
            dX = X[die]
            dex = ex[die]
            S0d = dex.sum()
            S1d = (dex[:, None] * dX).sum(0)
            S2d = (dex[:, None, None] * dX[:, :, None] * dX[:, None, :]).sum(0)
            g += dX.sum(0)
            ll += eta[die].sum()
            for l in range(m):
                c = l / m
                S0 = S0r - c * S0d
                S1 = S1r - c * S1d
                S2 = S2r - c * S2d
                zz = S1 / S0
                g -= zz
                H -= (S2 / S0 - np.outer(zz, zz))
                ll -= math.log(S0)
        return g, H, ll

    for _ in range(max_iter):
        g, H, _ = ghl(beta)
        try:
            step = np.linalg.solve(-H, g)
        except np.linalg.LinAlgError:
            step = np.linalg.lstsq(-H, g, rcond=None)[0]
        beta = beta + step
        if np.max(np.abs(step)) < tol:
            break
    g, H, ll = ghl(beta)
    try:
        se = np.sqrt(np.diag(np.linalg.inv(-H)))
    except np.linalg.LinAlgError:
        warnings.warn("survival._fit: -H is singular at the standard-error step; "
                       "returning missing (NaN) estimates for this stratum")
        beta = np.full(p, np.nan)
        se = np.full(p, np.nan)
    return beta, se, ll


def cox_efron(X, t, d, max_iter=100, tol=1e-9):
    """Cox coefficients and their standard errors. X is n by p, t the time, d the event flag."""
    beta, se, _ = _fit(X, t, d, max_iter=max_iter, tol=tol)
    return beta, se


def _wald(b, se):
    """Two-sided Wald p value of a coefficient b with standard error se."""
    return float(2 * norm.sf(abs(b / se)))


def hr_per_sd(df, col, covars=("stage",), time="OS_weeks", event="event"):
    """Hazard ratio per one standard deviation of col, with the covariates entered as given.

    Rows with a missing exposure, time, event or covariate are dropped, and so are rows with a
    follow-up of zero or less; the exposure is then z-scored over the rows that remain, so that the
    standard deviation is the one of the analyzed set and not of the whole cohort.
    """
    covars = list(covars)
    x = df.dropna(subset=[col, time, event] + covars)
    x = x[x[time] > 0]
    z = (x[col] - x[col].mean()) / x[col].std()
    X = np.column_stack([np.asarray(z, dtype=float)]
                        + [np.asarray(x[c], dtype=float) for c in covars])
    b, se = cox_efron(X, np.asarray(x[time], dtype=float),
                      np.asarray(x[event], dtype=float).astype(int))
    return dict(HR=math.exp(b[0]), lo=math.exp(b[0] - 1.96 * se[0]), hi=math.exp(b[0] + 1.96 * se[0]),
                p=_wald(b[0], se[0]), logHR=float(b[0]), se=float(se[0]),
                n=int(len(x)), events=int(x[event].sum()))


def meta_fixed(results):
    """Inverse-variance fixed-effect pooling of results on the log hazard scale.

    Each entry needs a standard error and an effect: logHR as written by hr_per_sd, or beta as
    written by histology_interaction, so that the same function pools hazard ratios and interaction
    terms. Patient and event counts are summed when the entries carry them.
    """
    eff = np.array([float(r["logHR"] if "logHR" in r else r["beta"]) for r in results])
    w = np.array([1 / float(r["se"]) ** 2 for r in results])
    lg = float(np.sum(w * eff) / w.sum())
    se = math.sqrt(1 / w.sum())
    return dict(HR=math.exp(lg), lo=math.exp(lg - 1.96 * se), hi=math.exp(lg + 1.96 * se),
                p=_wald(lg, se), logHR=lg, se=se,
                n=int(sum(r.get("n", 0) for r in results)),
                events=int(sum(r.get("events", 0) for r in results)))


def histology_interaction(df, col, covars=("stage",), time="OS_weeks", event="event"):
    """Does the effect of col on survival differ between squamous and adenocarcinoma disease.

    Model, on LUSC and LUAD together: z + LUAD indicator + z x LUAD + covariates, where z is col
    standardized once over both histologies so that the interaction term is the per-standard-
    deviation difference in log hazard between the two. The reported p is the Wald p on that term;
    a likelihood-ratio p is returned beside it as a check. Pool two cohorts with meta_fixed.
    """
    covars = list(covars)
    x = df.dropna(subset=[col, time, event, "histology"] + covars)
    x = x[(x[time] > 0) & (x.histology.isin(["LUSC", "LUAD"]))].copy()
    x["_z"] = (x[col] - x[col].mean()) / x[col].std()
    x["_luad"] = (x.histology == "LUAD").astype(float)
    x["_zx"] = x._z * x._luad
    full = ["_z", "_luad", "_zx"] + covars
    reduced = ["_z", "_luad"] + covars
    Xf = np.column_stack([np.asarray(x[c], dtype=float) for c in full])
    Xr = np.column_stack([np.asarray(x[c], dtype=float) for c in reduced])
    t = np.asarray(x[time], dtype=float)
    d = np.asarray(x[event], dtype=float).astype(int)
    bf, sef, llf = _fit(Xf, t, d)
    _, _, llr = _fit(Xr, t, d)
    i = full.index("_zx")
    lr = 2 * (llf - llr)
    return dict(p_wald=_wald(bf[i], sef[i]), beta=float(bf[i]), se=float(sef[i]),
                HR_ratio=math.exp(bf[i]), HR_LUSC=math.exp(bf[0]), HR_LUAD=math.exp(bf[0] + bf[i]),
                chi2=float(lr), p_lr=float(chi2.sf(lr, 1)),
                n=int(len(x)), n_LUSC=int((x._luad == 0).sum()), n_LUAD=int((x._luad == 1).sum()),
                events=int(d.sum()))
