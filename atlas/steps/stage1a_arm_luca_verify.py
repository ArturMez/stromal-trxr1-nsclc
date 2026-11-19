#!/usr/bin/env python
"""
Atlas stage 01, step arm_luca_verify: verify the LuCA arm written by stage1a_arm_luca.py.

Paper items: none directly; it checks the LuCA part of the atlas behind Figure 2, Figure 3,
Supplementary Figures 2 to 4 and Supplementary Data 1. No random numbers are drawn.

Reads stage1a/arm_luca.h5ad; every check re-derives its number from that file, not from the build
run (shape 1,061,555 x 2,058, 42 LUSC and 155 LUAD donors, unique cell ids, obs schema, compartments).
Prints PASS or FAIL per check, and writes the donor table (arm_luca_donor_table.csv) and the label map
(arm_luca_label_map.csv) to stage1a/. Run: python atlas/steps/stage1a_arm_luca_verify.py
(normally through atlas/01_arms_luca_lambrechts.py).
"""
from __future__ import annotations

import sys

import anndata as ad
import numpy as np
import pandas as pd

import os as _os  # noqa: E402
import sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # atlas/
from common import ROOT, SCRNA, SCRATCH, DATA, RSCRIPT  # noqa: E402,F401  paths from config.json
META = f"{ROOT}/output/RNA/meta_dataset"
OUT = f"{META}/stage1a"
H5AD = f"{OUT}/arm_luca.h5ad"

SCHEMA_COLUMNS = [
    "cell_id", "source_key", "study", "donor_id", "library_id", "histology",
    "histology_source", "histology_confidence", "tissue_class", "sampling_site",
    "platform", "modality", "enrichment", "composition_safe", "treatment",
    "block_id", "level1_compartment", "level2_celltype", "level3_celltype_native",
    "celltype_source", "celltype_confidence", "is_malignant", "malignancy_score",
    "n_counts", "n_genes", "counts_are_estimates", "doublet_score", "sex", "age",
    "stage", "smoking",
]


def main() -> int:
    """Run every check on arm_luca.h5ad and print the results; returns 0 when all pass."""
    pd.set_option("display.width", 250)
    pd.set_option("display.max_rows", 400)
    a = ad.read_h5ad(H5AD, backed="r")
    o = a.obs
    X = a.X

    ok, bad = [], []

    def check(name, cond, detail=""):
        """Record and print one PASS or FAIL check."""
        (ok if cond else bad).append(f"{name}: {detail}")
        print(f"{'PASS' if cond else 'FAIL'}  {name}  {detail}")

    check("shape", a.shape == (1_061_555, 2_058), f"{a.shape[0]:,} x {a.shape[1]:,}")
    nd = o.groupby("histology", observed=True)["donor_id"].nunique()
    check("donors LUSC 42", int(nd.get("LUSC", 0)) == 42, str(int(nd.get("LUSC", 0))))
    check("donors LUAD 155", int(nd.get("LUAD", 0)) == 155, str(int(nd.get("LUAD", 0))))
    check("cell_id unique", o.index.is_unique and o["cell_id"].is_unique,
          f"{o.index.nunique():,} distinct")
    check("cell_id form source|library|barcode",
          bool(o["cell_id"].str.count(r"\|").eq(2).all()) and
          bool(o["cell_id"].str.startswith("LuCA_census|").all()), "")
    check("donor_id never empty",
          o["donor_id"].notna().all() and (o["donor_id"].astype(str) != "").all(), "")
    check("every schema column present",
          all(c in o.columns for c in SCHEMA_COLUMNS),
          f"missing {[c for c in SCHEMA_COLUMNS if c not in o.columns]}")
    check("platform on every cell", o["platform"].notna().all(), "")
    check("modality on every cell",
          o["modality"].notna().all() and set(o["modality"].unique()) == {"single_cell"}, "")
    cs = o.groupby("enrichment", observed=True)["composition_safe"].agg(["mean", "size"])
    check("composition_safe True only for whole_tissue",
          bool(cs.loc["whole_tissue", "mean"] == 1.0) and
          bool((cs.drop(index="whole_tissue")["mean"] == 0).all()), "")
    check("block_id is study_platform",
          bool((o["study"].astype(str) + "_" + o["platform"].astype(str)
                == o["block_id"].astype(str)).all()), "")
    check("var index is ensembl", bool(a.var.index.str.match(r"^ENSG\d{11}$").all()),
          f"{a.n_vars} genes")
    check("var index unique", a.var.index.is_unique, "")
    check("RGS5 present via ensembl join",
          "ENSG00000143248" in set(a.var.index),
          "stored in LuCA as RGS5_ENSG00000143248, a symbol join drops it")

    # the matrix must be raw counts, not the log-normalized X of the source
    sub = X[:20000].tocsr() if hasattr(X, "tocsr") else X[:20000]
    d = np.asarray(sub.data)
    check("X is counts, not lognorm", bool((d >= 0).all() and d.max() > 50),
          f"max {d.max():.1f} over the first 20,000 cells")
    frac = float((d != np.floor(d)).mean())
    check("first 20,000 cells are integer counts", frac == 0.0, f"{100 * frac:.2f} pct fractional")
    est = o.groupby("study", observed=True)["counts_are_estimates"].max()
    check("counts_are_estimates set on the two Smart-seq2 studies",
          set(est[est].index) == {"Maynard_Bivona_2020", "Guo_Zhang_2018"},
          ", ".join(sorted(est[est].index)))
    check("n_counts positive everywhere", bool((o["n_counts"] > 0).all()),
          f"min {int(o['n_counts'].min())}")
    check("no cell lost every panel gene", bool((o["n_genes_panel"] > 0).all()),
          f"min {int(o['n_genes_panel'].min())} panel genes per cell")

    cmap = pd.read_csv(f"{META}/build-plan_compartment_map.csv")
    m = dict(zip(cmap["census_cell_type"], cmap["level1_compartment"]))
    agree = o["census_cell_type"].astype(str).map(m).values == o["level1_compartment"].astype(str).values
    check("level1 matches the explicit map on every cell", bool(agree.all()),
          f"{int((~agree).sum())} disagreements")
    mal = o.loc[o["census_cell_type"] == "malignant cell", "level1_compartment"].value_counts()
    check("malignant cells are not in T_NK",
          set(mal[mal > 0].index) == {"epithelial_malignant"},
          f"{int(mal.get('epithelial_malignant', 0)):,} malignant cells in epithelial_malignant")
    mast = o.loc[o["census_cell_type"] == "mast cell", "level1_compartment"].value_counts()
    check("mast cells are not in T_NK", set(mast[mast > 0].index) == {"other"},
          f"{int(mast.get('other', 0)):,} mast cells in other")

    print("\n=== per-block donor counts ===")
    b = o.groupby("block_id", observed=True).apply(
        lambda g: pd.Series({
            "cells": len(g), "donors": g.donor_id.nunique(),
            "LUSC": g.loc[g.histology == "LUSC", "donor_id"].nunique(),
            "LUAD": g.loc[g.histology == "LUAD", "donor_id"].nunique(),
            "NSCLC_unspecified": g.loc[g.histology == "NSCLC_unspecified", "donor_id"].nunique(),
            "normal": g.loc[g.histology == "normal", "donor_id"].nunique(),
            "COPD": g.loc[g.histology == "COPD", "donor_id"].nunique(),
            "libraries": g.library_id.nunique(),
            "composition_safe_cells": int(g.composition_safe.sum()),
        }), include_groups=False).sort_values("cells", ascending=False)
    print(b.to_string())
    b.to_csv(f"{OUT}/arm_luca_block_table.csv")

    don = o.groupby(["study", "block_id", "donor_id", "histology"], observed=True).agg(
        cells=("donor_id", "size"), libraries=("library_id", "nunique"),
        patient_id=("patient_id", "first"),
        enrichment=("enrichment", lambda s: ",".join(sorted(set(map(str, s))))),
        tissue_class=("tissue_class", lambda s: ",".join(sorted(set(map(str, s))))),
        n_counts_median=("n_counts", "median"))
    comp = pd.crosstab(o["donor_id"], o["level1_compartment"])
    comp.columns = [str(c) for c in comp.columns]
    don = don.join(comp, on="donor_id")
    don.reset_index().to_csv(f"{OUT}/arm_luca_donor_table.csv", index=False)
    print(f"\ndonor table: {len(don)} donor x block rows -> {OUT}/arm_luca_donor_table.csv")

    lab = o.groupby(["level3_celltype_native", "census_cell_type", "level1_compartment",
                     "level2_celltype"], observed=True).size().rename("cells").reset_index()
    lab = lab[lab["cells"] > 0]
    lab.to_csv(f"{OUT}/arm_luca_label_map.csv", index=False)
    print(f"label map: {len(lab)} rows -> {OUT}/arm_luca_label_map.csv")

    print("\n=== donors, patients, libraries, cells by histology ===")
    print(o.groupby("histology", observed=True).agg(
        donors=("donor_id", "nunique"), patients=("patient_id", "nunique"),
        libraries=("library_id", "nunique"), cells=("histology", "size")).to_string())
    print("\n=== compartment x histology ===")
    print(pd.crosstab(o["level1_compartment"], o["histology"],
                      margins=True, margins_name="all").to_string())
    print("\n=== compartment x histology, composition_safe cells only ===")
    s = o[o["composition_safe"]]
    print(pd.crosstab(s["level1_compartment"], s["histology"],
                      margins=True, margins_name="all").to_string())

    print(f"\n{len(ok)} checks passed, {len(bad)} failed")
    for f in bad:
        print(f"  FAILED {f}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
