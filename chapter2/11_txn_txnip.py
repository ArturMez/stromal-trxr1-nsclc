#!/usr/bin/env python3
"""
Chapter 2, step 11: TXN and TXNIP within each cell set, the direct test that thioredoxin and TXNIP regulate each
other (Maimaiti et al., Antioxid Redox Signal 2026, PMID 41736403: thioredoxin suppresses TXNIP).

ITEMS: Supplementary Figure 2d (file panel_txn_txnip) and the sentence on TXN and TXNIP in the thioredoxin
      section; the sheet TXN and TXNIP of Supplementary Data 1 (through 08_supptable2.py s3 and
      10_supplementary_data1.py).

WHAT: per cell set (fib, mac, tnk, and can when its working set exists):
        - the split: TXN is detected far more often than TXNRD1 (in 34 to 87 % of the cells of a cell set), so the
          cells are split by TXN level instead of detection.
          TXN-high is a TXN count fraction above the median of the cell's batch x histology stratum (a cell at the
          median counts as low; where fewer than half of the cells of a stratum have TXN counts, this is detection);
        - effect: the mean paired difference in the TXNIP score (the within-batch z of the TXNIP count fraction, the
          single-gene score of 06), TXN-high minus a TXN-low cell of the same stratum matched on genes detected (the
          matcher of 03 and 06, SEED = 6, the same cap of 2,000 cells per batch x histology);
        - z, normal p and empirical p from 1,000 permutations of the TXN-high label within patient, the matcher re-run
          in each (the permutation stream of 06);
        - calibration, the analogue of the comparison-gene rule of 06: every tested gene whose mean count fraction lies
          within 20 % of TXN's (widened in steps of 5 % until 200 genes) is split at its own stratum medians and
          paired the same way; the percentile of the TXN effect among their TXNIP effects (TXN, TXNIP and TXNRD1 are
          not comparison genes).
FROM: work_root/intermediate/chapter2/work_<cell set>/ (02_working_sets.py).
TO:   work_root/output/chapter2/txn_txnip.csv, txn_txnip_numbers.csv (the counts the text quotes) and
      txn_txnip_sets.gmt (each with a sibling .md), the panel
      panel_txn_txnip.png/.pdf (one row per cell set: the TXN effect on TXNIP against the range of the comparison
      genes), and work_root/logs/chapter2/11_<cell set>.log.
WHY:  06 splits by TXNRD1 detection, and TXNIP was lower in TXNRD1-positive cancer cells (1st percentile of the
      comparison genes). This asks the question of the paper directly, with the same pairing, permutations and
      calibration as 06.
RUN:  python3 chapter2/11_txn_txnip.py [--cells fib mac tnk can] [--nperm 1000] [--only-panel]
      --only-panel redraws the panel and rewrites the counts from txn_txnip.csv without any computation.
      Seeds: SEED = 6 (cap, observed pairing, permutation stream, comparison-gene splits); 1,000 permutations.
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ch2lib as L  # noqa: E402
from redox import paths  # noqa: E402
from redox import plotting as RP  # noqa: E402  (chooses the Agg backend)
import matplotlib.pyplot as plt  # noqa: E402

SEED = L.SEED
CAP = 2000
MINCTRL = 200
TOL0 = 0.20
LOGS = paths.WORK_ROOT / "logs" / "chapter2"
EXCLUDE = {"TXN", "TXNIP", L.SPLITTER}
COLUMNS = [
    ("cell", "fib, mac, tnk or can"), ("cell_type", "cell set as printed"),
    ("cells", "cells after the cap"), ("patients", "patients"), ("batches", "sequencing batches"),
    ("txn_detected", "fraction of cells with TXN detected"), ("txnip_detected", "fraction of cells with TXNIP detected"),
    ("matched_pairs", "pairs of the observed pairing, TXN-high with TXN-low"),
    ("effect", "mean paired difference in the TXNIP score, TXN-high minus TXN-low"),
    ("z", "(effect minus the permutation mean) / permutation SD"), ("p_gauss", "two-sided normal p of z"),
    ("p_emp", "empirical permutation p (centered)"), ("n_replicates", "permutations used"),
    ("n_controls", "control genes (mean count fraction within the tolerance of TXN's)"),
    ("tolerance", "the window around TXN's mean count fraction"),
    ("percentile", "percentile of the TXN effect among the control effects"),
    ("ctrl_p2_5", "2.5th percentile of the control effects"), ("ctrl_median", "their median"),
    ("ctrl_p97_5", "their 97.5th percentile"),
    ("reading", "unusual, lower (<= 2.5th), unusual, higher (>= 97.5th), ordinary, or not testable")]


def write_gmt():
    """Write the one-set gmt (TXNIP alone) that is scored; returns its path."""
    path = L.out("txn_txnip_sets.gmt")
    path.write_text("G_TXNIP\tTXNIP\tTXNIP\n")
    L.sidecar(path, "The one set 11_txn_txnip.py scores: TXNIP alone (fields: set id, label, genes).", [])
    return path


def median_split(f, groups):
    """True where f lies above the median of its stratum."""
    high = np.zeros(len(f), bool)
    for ix in groups:
        high[ix] = f[ix] > np.median(f[ix])
    return high


def run(cell, a, gmt):
    """The TXN split and its calibration for one cell set (see the header); returns one row of txn_txnip.csv."""
    spec = L.cell_spec(cell)
    log = L.Log()
    S = L.load_set(cell, CAP)
    N, X, tot, GI, genes = S["N"], S["X"], S["tot"], S["GI"], S["genes"]
    depth, strat, pidx, keepg = S["depth"], S["strat"], S["pidx"], S["keepg"]
    log("=" * 100)
    log(f"TXN AND TXNIP.  {spec['label']}.  nperm {a.nperm}, cap {CAP}, seed {SEED}")
    rec = dict(cell=cell, cell_type=spec["label"], cells=N, patients=len(np.unique(S["pat"])),
               batches=len(S["studies"]))
    if "TXN" not in GI or "TXNIP" not in S["UNI"]:
        rec["reading"] = "not testable (TXN or TXNIP fails the detection rule)"
        log("  " + rec["reading"])
        log.write(LOGS / f"11_{cell}.log")
        return rec
    Xn = X.multiply(1.0 / np.maximum(tot, 1).reshape(-1, 1)).tocsc()
    col = lambda g: np.asarray(Xn[:, g].todense()).ravel()      # count fraction of gene column g in every cell
    rec.update(txn_detected=float((col(GI["TXN"]) > 0).mean()), txnip_detected=float((col(GI["TXNIP"]) > 0).mean()))
    log(f"  {N:,} cells, {rec['patients']} patients, {rec['batches']} batches; TXN detected in "
        f"{rec['txn_detected']:.3f}, TXNIP in {rec['txnip_detected']:.3f}")
    P = L.score_programs(S, gmt, minset=1)
    k = P["names"].index("G_TXNIP")
    scorer = L.Scorer.build(P, S)
    strata = L.Strata(strat)
    groups = [np.where(strat == s)[0] for s in np.unique(strat)]

    E = median_split(col(GI["TXN"]), groups)
    L.check_matcher(depth, E, strat)
    aa, bb = L.match_fast(depth, E, strata, np.random.default_rng(SEED))
    obs = float(scorer.base(aa, bb)[0][k])
    r = np.random.default_rng(SEED)
    nul = []
    for i in range(a.nperm):
        pa, pb = L.match_fast(depth, L.perm_expo(E, pidx, r), strata, r)
        if len(pa) < L.MINPAIR:
            continue
        nul.append(float(scorer.base(pa, pb)[0][k]))
        if (i + 1) % 250 == 0:
            log(f"    replicate {i + 1} / {a.nperm}   [{log.elapsed()}]")
    nul = np.array(nul)
    mu, sd = nul.mean(), nul.std()
    z = (obs - mu) / sd if sd > 0 else np.nan
    rec.update(matched_pairs=len(aa), effect=obs, z=z, p_gauss=2 * stats.norm.sf(abs(z)),
               p_emp=(np.sum(np.abs(nul - mu) >= abs(obs - mu)) + 1) / (len(nul) + 1), n_replicates=len(nul))
    log(f"  pairs {len(aa):,}; TXNIP effect {obs:+.4f}, z {z:+.2f}   [{log.elapsed()}]")

    # comparison genes: mean count fraction within TOL0 (20 %) of TXN's, the window widened by 5 % until at least MINCTRL genes
    mfrac = np.asarray(Xn.mean(0)).ravel()
    t = mfrac[GI["TXN"]]
    tol = TOL0
    while True:
        pool = [g for g in keepg if abs(mfrac[g] - t) <= tol * t and genes[g] not in EXCLUDE]
        if len(pool) >= MINCTRL or tol >= 1.0:
            break
        tol += 0.05
    ref = []
    for g in pool:
        pa, pb = L.match_fast(depth, median_split(col(g), groups), strata, np.random.default_rng(SEED))
        if len(pa) >= L.MINPAIR:
            ref.append(float(scorer.base(pa, pb)[0][k]))
    ref = np.array([x for x in ref if np.isfinite(x)])
    pct = 100 * float((ref < obs).mean())
    rec.update(n_controls=len(ref), tolerance=tol, percentile=pct, ctrl_p2_5=np.percentile(ref, 2.5),
               ctrl_median=np.median(ref), ctrl_p97_5=np.percentile(ref, 97.5),
               reading="unusual, lower" if pct <= 2.5 else "unusual, higher" if pct >= 97.5 else "ordinary")
    log(f"  controls {len(ref):,} genes within {tol:.0%} of TXN's mean count fraction; percentile {pct:.1f} "
        f"({rec['reading']})   [{log.elapsed()}]")
    log.write(LOGS / f"11_{cell}.log")
    return rec


BLUE, RED, AXIS, BAND = "#2a78d6", "#e34948", "#c3c2b7", "#efeee9"     # the colors of 06_thioredoxin.py


def panel(T):
    """One row per cell set: the TXNIP difference of the TXN split (circle) against the range of 95 % of the control
    genes (gray bar, tick at their median), drawn as the panels of 06_thioredoxin.py."""
    RP.style()
    plt.rcParams["axes.edgecolor"] = AXIS
    T = T.dropna(subset=["effect"]).reset_index(drop=True)
    n = len(T)
    fig = plt.figure(figsize=(7.2, 0.95 + 0.145 * n))
    ax = fig.add_axes([0.355, 0.42, 0.30, 0.42])
    ax.set_xlim(-0.25, 0.25)
    ax.set_ylim(n - 0.5, -0.5)
    ax.axvline(0, color=AXIS, lw=0.6, zorder=1)
    ax.set_yticks(range(n))
    ax.set_yticklabels([f"{r.cell_type}, {int(r.matched_pairs):,} pairs" for r in T.itertuples()], fontsize=6.6)
    ax.tick_params(axis="y", length=0)
    ax.set_xticks([-0.2, -0.1, 0, 0.1, 0.2])
    ax.set_xticklabels(["−0.2", "−0.1", "0", "0.1", "0.2"])
    ax.tick_params(axis="x", labelsize=6, length=2)
    ax.spines[["top", "right", "left"]].set_visible(False)
    for j, r in enumerate(T.itertuples()):
        ax.barh(j, r.ctrl_p97_5 - r.ctrl_p2_5, left=r.ctrl_p2_5, height=0.5, color=BAND, edgecolor=AXIS, lw=0.5, zorder=2)
        ax.plot([r.ctrl_median] * 2, [j - 0.25, j + 0.25], color=RP.INK2, lw=1.0, zorder=3, solid_capstyle="butt")
        edge = BLUE if r.reading.endswith("higher") else RED if r.reading.endswith("lower") else RP.MUTED
        ax.scatter(r.effect, j, facecolor=edge if edge != RP.MUTED else "white", edgecolor=edge, lw=0.9, s=24, zorder=4,
                   clip_on=False)
    ax.set_xlabel("difference in TXNIP score, cells with high TXN minus partner cells", fontsize=6.4, labelpad=2)
    fig.text(0.03, 0.97, "TXNIP in cells with high TXN, against genes expressed as strongly as TXN", fontsize=7.4,
             fontweight="bold", va="top")
    RP.save(fig, 2, "panel_txn_txnip")
    print("wrote panel_txn_txnip.png and .pdf")


def numbers(T):
    """The counts the text quotes (txn_txnip_numbers.csv)."""
    t = T.dropna(subset=["effect"])
    items = [("cell_sets_tested", len(t), "cell sets with a testable TXN split"),
             ("txnip_lower_cells", int((t.z <= -L.Z95).sum()),
              "cell sets in which TXNIP is lower in the cells with high TXN at z <= -1.96 ('in all four cell types')"),
             ("txnip_beyond_range_cells", int(t.reading.str.startswith("unusual").sum()),
              "cell sets in which this difference lies beyond the 2.5th to 97.5th range of the control genes (0; "
              "'within the range obtained with other genes expressed as strongly as TXN')"),
             ("tolerance_max_percent", int(round(100 * t.tolerance.max())),
              "widest window around TXN's mean count fraction that held 200 control genes, percent"),
             ("controls_min", int(t.n_controls.min()), "control genes, smallest over the cell sets")]
    N = pd.DataFrame(items, columns=["item", "value", "meaning"])
    o = L.out("txn_txnip_numbers.csv")
    N.to_csv(o, index=False)
    L.sidecar(o, "Counts quoted with the TXN split (Supplementary Figure 2d).",
              [("item", "key"), ("value", "value"), ("meaning", "what it counts")])
    print(N.to_string(index=False))


def main():
    """Run the TXN split for the cell sets named with --cells, write the table, the counts and the panel."""
    ap = argparse.ArgumentParser(description=__doc__.split("RUN:")[0].strip(),
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    have = [c for c in L.EXTRA_CELL_SETS if (L.work_dir(c) / "genes.npy").exists()]
    ap.add_argument("--cells", nargs="*", default=list(L.CELL_SETS) + have,
                    choices=list(L.CELL_SETS) + list(L.EXTRA_CELL_SETS))
    ap.add_argument("--nperm", type=int, default=1000)
    ap.add_argument("--only-panel", action="store_true", help="redraw the panel from txn_txnip.csv, no computation")
    a = ap.parse_args()
    LOGS.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    o = L.out("txn_txnip.csv")
    if a.only_panel:
        T = pd.read_csv(o, float_precision="round_trip")
        numbers(T)
        panel(T)
        return
    gmt = write_gmt()
    T = pd.DataFrame([run(c, a, gmt) for c in a.cells])
    T.to_csv(o, index=False)
    L.sidecar(o, "TXNIP under the TXN split (TXN-high against matched TXN-low cells of the same stratum), per cell set, "
              "calibrated against genes expressed as strongly as TXN (11_txn_txnip.py)", COLUMNS)
    with pd.option_context("display.width", 220):
        print(T[[c for c in ("cell", "cells", "patients", "matched_pairs", "effect", "z", "p_emp", "percentile",
                             "n_controls", "reading") if c in T]].to_string(index=False, float_format=lambda v: f"{v:.4g}"))
    numbers(T)
    panel(T)
    print(f"11_txn_txnip.py done in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
