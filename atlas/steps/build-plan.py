#!/usr/bin/env python
"""
Atlas stage 01, step build_plan: the explicit compartment map of the meta-dataset build.

Writes the map from CELLxGENE Census cell-type labels to the level-1 compartments, and, when the
cached LuCA obs parquet exists, the donor and block tables with assertions on the numbers quoted
in the build plan. Touches no matrix and downloads nothing. Paper items: none directly; the map
is read by later steps of the atlas that feeds Figure 2, Figure 3, Supplementary Figures 2 to 4
and Supplementary Data 1. No random numbers are drawn.

Run: python atlas/steps/build-plan.py   (normally through atlas/01_arms_luca_lambrechts.py)

The Census obs cache is optional here. Without it the script writes build-plan_compartment_map.csv,
the only one of its outputs any later step reads, and skips the donor and block tables and their
assertions. The cache (CELLxGENE Census 2023-12-15) is not part of the package.

Outputs (all in output/RNA/meta_dataset/):
  build-plan_donor_table.csv   one row per Census donor, compartment counts
  build-plan_block_generated.csv  study x platform blocks holding both histologies
  build-plan_compartment_map.csv  the explicit 33-label map, no substring matching
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

import os as _os  # noqa: E402
import sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # atlas/
from common import ROOT as _ROOT  # noqa: E402  paths from config.json, see atlas/common.py
ROOT = Path(_ROOT)
OUT = ROOT / "output/RNA" / "meta_dataset"
CACHE = OUT / "_luca_obs_cache.parquet"
LUCA_DATASET_ID = "1e6a6ef9-7ec9-4c90-bbfb-2ad3c3165fd1"

# ---------------------------------------------------------------------------
# The compartment map.  EXPLICIT, exhaustive, one entry per Census label.
# Substring matching is forbidden here: "malignant cell" and "mast cell" both
# contain the string "t cell" and a substring mapper silently sent 33,234 LUSC
# malignant cells into T/NK in an early version of the acquisition.
# ---------------------------------------------------------------------------
COMPARTMENT = {
    "CD4-positive, alpha-beta T cell": "T_NK",
    "CD8-positive, alpha-beta T cell": "T_NK",
    "regulatory T cell": "T_NK",
    "natural killer cell": "T_NK",
    "B cell": "B_plasma",
    "plasma cell": "B_plasma",
    "alveolar macrophage": "myeloid",
    "macrophage": "myeloid",
    "classical monocyte": "myeloid",
    "non-classical monocyte": "myeloid",
    "myeloid cell": "myeloid",
    "dendritic cell": "myeloid",
    "conventional dendritic cell": "myeloid",
    "CD1c-positive myeloid dendritic cell": "myeloid",
    "plasmacytoid dendritic cell": "myeloid",
    "neutrophil": "myeloid",
    "malignant cell": "epithelial_malignant",
    "epithelial cell of lung": "epithelial_malignant",
    "club cell": "epithelial_malignant",
    "multi-ciliated epithelial cell": "epithelial_malignant",
    "type I pneumocyte": "epithelial_malignant",
    "type II pneumocyte": "epithelial_malignant",
    "fibroblast of lung": "fibroblast",
    "bronchus fibroblast of lung": "fibroblast",
    "stromal cell": "fibroblast",
    "pericyte": "mural",
    "smooth muscle cell": "mural",
    "capillary endothelial cell": "endothelial",
    "endothelial cell of lymphatic vessel": "endothelial",
    "pulmonary artery endothelial cell": "endothelial",
    "vein endothelial cell": "endothelial",
    "mast cell": "other",
    "mesothelial cell": "other",
}
MACROPHAGE_STRICT = {
    "alveolar macrophage",
    "macrophage",
    "classical monocyte",
    "non-classical monocyte",
}
HISTOLOGY = {
    "lung adenocarcinoma": "LUAD",
    "squamous cell lung carcinoma": "LUSC",
    "non-small cell lung carcinoma": "NSCLC_unspecified",
    "normal": "normal",
    "chronic obstructive pulmonary disease": "COPD",
}

# Lambrechts donors are excluded from the Census pull and taken from the local
# E-MTAB-6149 / E-MTAB-6653 files instead: the local route gives +28.3 % cells,
# +101 % fibroblasts in 6149, the core/middle/edge/normal sampling site the
# Census structurally cannot carry, and the author 52-subtype taxonomy.  Keeping
# both routes would double-count 8 donors.
EXCLUDE_DONOR_PREFIX = "Lambrechts"


def write_compartment_map():
    """Write build-plan_compartment_map.csv (Census label, level-1 compartment, strict macrophage flag)."""
    pd.DataFrame(
        {"census_cell_type": list(COMPARTMENT),
         "level1_compartment": [COMPARTMENT[k] for k in COMPARTMENT],
         "macrophage_strict": [k in MACROPHAGE_STRICT for k in COMPARTMENT]}
    ).to_csv(OUT / "build-plan_compartment_map.csv", index=False)


def main() -> int:
    """Write the compartment map and, if the Census obs cache exists, the donor and block tables."""
    if not CACHE.exists():
        print(f"no Census obs cache at {CACHE}: writing the compartment map only")
        write_compartment_map()
        return 0

    obs = pd.read_parquet(CACHE)
    obs = obs[obs.dataset_id == LUCA_DATASET_ID]
    obs = obs[obs.is_primary_data]

    unmapped = set(obs.cell_type.astype(str)) - set(COMPARTMENT)
    assert not unmapped, f"unmapped Census cell types: {sorted(unmapped)}"

    obs = obs.assign(
        compartment=obs.cell_type.astype(str).map(COMPARTMENT),
        histology=obs.disease.astype(str).map(HISTOLOGY),
        # study = the donor_id prefix up to and including the publication year.
        # A naive rsplit on "_" splits one study into several blocks whenever a
        # donor label itself contains an underscore (Goveia..._patient_1,
        # Maynard..._TH103_2), which silently shrinks the block it matters most
        # for.  Donor ids without a year (UKIM-V, UKIM-V-2) keep their prefix.
        study=obs.donor_id.astype(str).str.extract(r"^(.*?_\d{4})", expand=False)
        .fillna(obs.donor_id.astype(str).str.split("_").str[0]),
        macrophage_strict=obs.cell_type.astype(str).isin(MACROPHAGE_STRICT),
    )

    # ---- assertions against the numbers the build plan quotes ----------------
    full = obs
    n_lusc = full.loc[full.histology == "LUSC", "donor_id"].nunique()
    n_luad = full.loc[full.histology == "LUAD", "donor_id"].nunique()
    assert (n_lusc, n_luad) == (45, 158), f"LuCA donors drifted: {n_lusc}/{n_luad}"
    fib = full[full.compartment == "fibroblast"]
    assert len(fib[fib.histology == "LUSC"]) == 3423
    assert len(fib[fib.histology == "LUAD"]) == 6383

    kept = obs[~obs.donor_id.astype(str).str.startswith(EXCLUDE_DONOR_PREFIX)]
    dropped = len(obs) - len(kept)
    assert dropped == 66296, f"Lambrechts cell count drifted: {dropped}"

    # ---- donor table --------------------------------------------------------
    donor = (
        kept.pivot_table(
            index=["study", "donor_id", "histology", "assay"],
            columns="compartment",
            values="cell_type",
            aggfunc="size",
            fill_value=0,
        )
        .reset_index()
    )
    donor["n_cells"] = donor.select_dtypes("number").sum(axis=1)
    donor["macrophage_strict"] = (
        kept[kept.macrophage_strict].groupby("donor_id").size().reindex(donor.donor_id).fillna(0).astype(int).values
    )
    donor["block_id"] = donor.study + " | " + donor.assay.astype(str)
    donor.to_csv(OUT / "build-plan_donor_table.csv", index=False)

    # ---- blocks holding both histologies ------------------------------------
    tumour = donor[donor.histology.isin(["LUAD", "LUSC"])]
    rows = []
    comps = ["epithelial_malignant", "fibroblast", "mural", "endothelial",
             "myeloid", "T_NK", "B_plasma"]
    for block, grp in tumour.groupby("block_id"):
        n = grp.groupby("histology").donor_id.nunique()
        if n.get("LUSC", 0) == 0 or n.get("LUAD", 0) == 0:
            continue
        row = {"block_id": block,
               "LUSC_donors": int(n.get("LUSC", 0)),
               "LUAD_donors": int(n.get("LUAD", 0))}
        for c in comps:
            for h in ("LUSC", "LUAD"):
                sub = grp[grp.histology == h]
                row[f"{c}_{h}_cells"] = int(sub[c].sum())
                row[f"{c}_{h}_donors_ge20"] = int((sub[c] >= 20).sum())
        rows.append(row)
    blocks = pd.DataFrame(rows).sort_values("LUSC_donors", ascending=False)
    blocks.to_csv(OUT / "build-plan_block_generated.csv", index=False)

    write_compartment_map()

    print(f"Census donors kept: LUSC {tumour[tumour.histology=='LUSC'].donor_id.nunique()}, "
          f"LUAD {tumour[tumour.histology=='LUAD'].donor_id.nunique()}")
    print(f"both-histology blocks: {len(blocks)}")
    print(blocks[["block_id", "LUSC_donors", "LUAD_donors",
                  "fibroblast_LUSC_donors_ge20", "fibroblast_LUAD_donors_ge20",
                  "epithelial_malignant_LUSC_donors_ge20",
                  "epithelial_malignant_LUAD_donors_ge20"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
