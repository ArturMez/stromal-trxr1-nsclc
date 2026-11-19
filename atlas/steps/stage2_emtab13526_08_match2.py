#!/usr/bin/env python
"""
Atlas stage 02, step match: resolve each author sample group to a raw library.

Paper items: part of the atlas behind Figure 2, Figure 3, Supplementary Figures 2 to 4 and Supplementary Data 1.
Run: python atlas/steps/stage2_emtab13526_08_match2.py   (normally through atlas/02_arm_emtab13526.py).
No random numbers are drawn.

Exact equality of author total_counts to the raw UMI total is
the wrong statistic: the authors computed total_counts on their 25,718-gene
filtered matrix, so author total_counts <= raw UMI, with equality only when the
dropped genes carried nothing. Containment of barcodes is necessary but not
sufficient, because 10x barcode sets overlap.

The discriminating statistic used here is the ratio author_total_counts /
raw_umi over the barcodes a library contains. For the true library the ratio
sits just below 1 for essentially every cell. For a library that merely shares
barcodes it is scattered, and frequently above 1, which is outright impossible
for the true library.

Score = share of the group's cells with 0.5 <= ratio <= 1.0 in that library.
A global one-to-one assignment is then solved with the Hungarian algorithm so
no library can be claimed twice, and the result is compared against the greedy
per-group argmax.

Writes: output/RNA/meta_dataset/stage2/fingerprint_assignment.csv
"""
import glob
import os

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

import os as _os  # noqa: E402
import sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # atlas/
from common import ROOT, SCRNA, SCRATCH, DATA, RSCRIPT  # noqa: E402,F401  paths from config.json
S2 = os.path.join(ROOT, "output/RNA/meta_dataset/stage2")
FP = os.path.join(S2, "fingerprint")

pd.set_option("display.width", 260)
pd.set_option("display.max_rows", 400)


def main():
    """Score every group against every library, solve the one-to-one assignment (Hungarian) and write fingerprint_scores.csv and fingerprint_assignment.csv."""
    obs = pd.read_parquet(os.path.join(S2, "author_obs.parquet"),
                          columns=["patient_sample", "author_barcode",
                                   "total_counts", "n_genes"])
    obs["barcode"] = obs["author_barcode"].str.extract(r"^([ACGT]+-1)", expand=False)
    obs["total_counts"] = obs["total_counts"].astype(np.int64)
    obs["n_genes"] = obs["n_genes"].astype(np.int64)
    obs = obs.rename(columns={"patient_sample": "grp"})[
        ["grp", "barcode", "total_counts", "n_genes"]]

    fp = pd.concat([pd.read_parquet(f).assign(library_id=os.path.basename(f)[:-8])
                    for f in sorted(glob.glob(os.path.join(FP, "*.parquet")))],
                   ignore_index=True)

    size = obs.groupby("grp", observed=True).size().rename("n_cells")
    m = obs.merge(fp, on="barcode", how="inner")
    m["ratio"] = m["total_counts"] / m["raw_umi"]
    m["gene_ratio"] = m["n_genes"] / m["raw_ngenes"]
    m["plausible"] = (m["ratio"] >= 0.5) & (m["ratio"] <= 1.0) & \
                     (m["gene_ratio"] >= 0.5) & (m["gene_ratio"] <= 1.0)

    agg = (m.groupby(["grp", "library_id"], observed=True)
             .agg(n_hit=("ratio", "size"),
                  n_plausible=("plausible", "sum"),
                  median_ratio=("ratio", "median"),
                  frac_ratio_gt1=("ratio", lambda r: float((r > 1).mean())))
             .reset_index().merge(size.reset_index(), on="grp"))
    agg["containment"] = agg["n_hit"] / agg["n_cells"]
    agg["score"] = agg["n_plausible"] / agg["n_cells"]
    agg.to_csv(os.path.join(S2, "fingerprint_scores.csv"), index=False)

    grps = sorted(agg["grp"].unique())
    libs = sorted(fp["library_id"].unique())
    gi = {g: i for i, g in enumerate(grps)}
    li = {l: i for i, l in enumerate(libs)}
    M = np.zeros((len(grps), len(libs)))
    for _, r in agg.iterrows():
        M[gi[r["grp"]], li[r["library_id"]]] = r["score"]

    # greedy per-group argmax
    greedy = {g: libs[int(np.argmax(M[gi[g]]))] for g in grps}
    # global one-to-one
    ri, ci = linear_sum_assignment(-M)
    hung = {grps[a]: libs[b] for a, b in zip(ri, ci)}

    rows = []
    for g in grps:
        row = M[gi[g]]
        order = np.argsort(-row)
        rows.append({
            "grp": g,
            "n_cells": int(size[g]),
            "greedy_library": greedy[g],
            "greedy_score": round(float(row[order[0]]), 4),
            "runner_up": libs[order[1]],
            "runner_up_score": round(float(row[order[1]]), 4),
            "hungarian_library": hung[g],
            "hungarian_score": round(float(M[gi[g], li[hung[g]]]), 4),
            "agree": greedy[g] == hung[g],
        })
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(S2, "fingerprint_assignment.csv"), index=False)

    print(out.to_string(index=False))
    print("\n=== summary ===")
    print("groups:", len(out))
    print("greedy == hungarian:", int(out["agree"].sum()), "of", len(out))
    print("greedy assignment injective:", out["greedy_library"].nunique() == len(out))
    print("min winning score : %.4f" % out["greedy_score"].min())
    print("max runner-up     : %.4f" % out["runner_up_score"].max())
    print("groups with score >= 0.99 :", int((out["greedy_score"] >= 0.99).sum()))
    print("\nwrote", os.path.join(S2, "fingerprint_assignment.csv"))


if __name__ == "__main__":
    main()
