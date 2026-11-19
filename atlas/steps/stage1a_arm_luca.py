#!/usr/bin/env python
"""
Atlas stage 01, step arm_luca: the LuCA arm at the panel gene axis.
Builds output/RNA/meta_dataset/stage1a/arm_luca.h5ad.

Paper items: part of the atlas behind Figure 2, Figure 3, Supplementary Figures 2 to 4 and
Supplementary Data 1 (the largest of the 13 public studies; the Lambrechts copies are excluded here
and included once through E-MTAB-6149 and E-MTAB-6653). No random numbers are drawn.

Reads LuCA/luca_extended.h5ad (00_download.sh) through h5py only, one row block at a time, so the
17.36 GB file is never materialised, plus stage1a/panel_genes.csv and
build-plan_compartment_map.csv. Writes arm_luca.h5ad and the tables arm_luca_enrichment.csv,
arm_luca_counts_integrality.csv, arm_luca_block_table.csv, arm_luca_donor_table.csv and
arm_luca_label_map.csv in stage1a/. Peak RSS is reported at the end.

What it does, in order
  1.  is_primary_data == True AND donor_id does not start with "Lambrechts".
      Asserts the exclusion removes exactly 66,296 cells, leaves 1,061,555, and
      that 42 LUSC / 155 LUAD donors remain.  Stops on any mismatch.
  2.  Subsets to the 1,906-gene panel, joined on ENSEMBL id, never on symbol:
      LuCA stores 22 of 17,811 features as "<SYMBOL>_<ENSG>" and RGS5 is one of
      them, so a symbol join silently drops the marker that level1 mural, the
      vCAF signature and the pericyte disambiguation all rest on.
  3.  Maps the 33 Census cell_type labels to level1_compartment with the
      EXPLICIT table in build-plan_compartment_map.csv.  No substring matching:
      "malignant cell" and "mast cell" both contain "t cell".
  4.  Populates the obs schema of build-plan_obs_schema.csv.
  5.  Writes raw counts (layers/count) into .X.

Run: python atlas/steps/stage1a_arm_luca.py   (normally through atlas/01_arms_luca_lambrechts.py)
"""
from __future__ import annotations

import gc
import os
import re
import resource
import sys
import time

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

import os as _os  # noqa: E402
import sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # atlas/
from common import ROOT, SCRNA, SCRATCH, DATA, RSCRIPT  # noqa: E402,F401  paths from config.json
LUCA = f"{SCRNA}/LuCA/luca_extended.h5ad"
META = f"{ROOT}/output/RNA/meta_dataset"
OUT = f"{META}/stage1a"
H5AD = f"{OUT}/arm_luca.h5ad"

SOURCE_KEY = "LuCA_census"
ROW_BLOCK = 20_000
EXCLUDE_DONOR_PREFIX = "Lambrechts"

EXPECT_DROPPED = 66_296
EXPECT_CELLS = 1_061_555
EXPECT_LUSC_DONORS = 42
EXPECT_LUAD_DONORS = 155

# ---------------------------------------------------------------------------
# fixed vocabularies
# ---------------------------------------------------------------------------

HISTOLOGY = {
    "lung adenocarcinoma": "LUAD",
    "squamous cell lung carcinoma": "LUSC",
    "non-small cell lung carcinoma": "NSCLC_unspecified",
    "normal": "normal",
    "chronic obstructive pulmonary disease": "COPD",
}

# Census assay -> the platform vocabulary of the obs schema.  The coarse LuCA
# "platform" column is NOT usable here: it collapses Leader_Merad's four
# chemistries into one value "10x", and in the obs schema two
# Leader_Merad chemistries are two blocks of one study.
PLATFORM = {
    "10x 3' v1": "10x 3' v1",
    "10x 3' v2": "10x 3' v2",
    "10x 3' v3": "10x 3' v3",
    "10x 5' v1": "10x 5' v1",
    "BD Rhapsody Whole Transcriptome Analysis": "BD Rhapsody",
    "GEXSCOPE technology": "GEXSCOPE",
    "Smart-seq2": "Smart-seq2",
    "inDrop": "inDrop",
    "Drop-seq": "Drop-seq",
}

# (origin, tissue) -> tissue_class.  EXPLICIT pairs.  origin, not origin_fine:
# origin_fine relabels 187,728 cells from libraries literally named
# "<patient>_normal_adjacent" as tumor_primary (UKIM-V-2, six Leader_Merad
# datasets) and relabels Mayr_Schiller normal_adjacent as normal.  origin agrees
# with the library names everywhere except the five Kim_Lee EFFUSION libraries,
# where origin says tumor_metastasis and both origin_fine and the tissue
# ontology say pleural effusion; the (origin, tissue) pair resolves that.
TISSUE_CLASS = {
    ("tumor_primary", "lung"): "primary_tumour",
    ("normal_adjacent", "lung"): "adjacent_normal",
    ("normal", "lung"): "normal_lung",
    ("normal", "lymph node"): "normal_LN",
    ("tumor_metastasis", "lymph node"): "LN_met",
    ("tumor_metastasis", "brain"): "brain_met",
    ("tumor_metastasis", "pleural effusion"): "effusion",
    ("tumor_metastasis", "liver"): "liver_met",      
    ("tumor_metastasis", "adrenal tissue"): "adrenal_met",
    # the only cells LuCA leaves without an origin are the two Maynard "_NAT"
    # donors, whose donor label spells out normal adjacent tissue
    ("nan", "lung"): "adjacent_normal",
}

# ann_fine (LuCA's finest label) -> level2_celltype.  EXPLICIT.  Verified 1:1
# nested inside the Census cell_type, so level2 can never contradict level1.
# No CAF subtype is assigned: LuCA carries no myCAF/iCAF/apCAF annotation and
# inventing one from a coarse label is exactly the trap this schema exists to
# stop.  All fibroblasts get "fibroblast" and keep their native label.
LEVEL2 = {
    "T cell CD4": "CD4_T",
    "T cell CD4 dividing": "CD4_T",
    "T cell CD8 activated": "CD8_T",
    "T cell CD8 dividing": "CD8_T",
    "T cell CD8 effector memory": "CD8_T",
    "T cell CD8 naive": "CD8_T",
    "T cell CD8 terminally exhausted": "CD8_T",
    "T cell NK-like": "CD8_T",
    "T cell regulatory": "Treg",
    "NK cell": "NK",
    "NK cell dividing": "NK",
    "B cell": "B_cell",
    "B cell dividing": "B_cell",
    "Plasma cell": "plasma_cell",
    "Plasma cell dividing": "plasma_cell",
    "Macrophage": "macrophage",
    "Macrophage alveolar": "macrophage",
    "Monocyte classical": "monocyte",
    "Monocyte non-classical": "monocyte",
    "cDC1": "DC",
    "cDC2": "DC",
    "DC mature": "DC",
    "pDC": "pDC",
    "Myeloid dividing": "myeloid_other",
    "Neutrophils": "neutrophil",
    "Mast cell": "mast",
    "Cancer cells": "malignant",
    "Alveolar cell type 1": "AT1",
    "Alveolar cell type 2": "AT2",
    "Transitional Club/AT2": "club",
    "Club": "club",
    "Ciliated": "ciliated",
    "ROS1+ healthy epithelial": "epithelial_other",
    "Fibroblast adventitial": "fibroblast",
    "Fibroblast alveolar": "fibroblast",
    "Fibroblast peribronchial": "fibroblast",
    "Stromal dividing": "fibroblast",
    "Pericyte": "pericyte",
    "Smooth muscle cell": "SMC",
    "Endothelial cell venous": "endothelial",
    "Endothelial cell capillary": "endothelial",
    "Endothelial cell arterial": "endothelial",
    "Endothelial cell lymphatic": "endothelial_lymphatic",
    "Mesothelial": "mesothelial",
}

# enrichment.  The Census records no sort field, so this is set from the two
# things that are actually recorded: the LuCA sub-dataset name, and the observed
# compartment composition of that sub-dataset.  A dataset whose non-immune
# fraction (epithelial + fibroblast + mural + endothelial) is under 3 pct cannot
# be a whole-tissue digest; the observed values fall into two clean groups,
# 0.0-1.2 pct and 7.4-67.8 pct, so the cut is not near any dataset.
NON_IMMUNE_SORT_CUT = 3.0
ENDOTHELIAL_ENRICH_CUT = 20.0
ENRICHMENT_NAMED = {
    # T-cell sorted study: 99.8 pct T/NK
    "Guo_Zhang_2018": ("T_sort", "brief; 99.8 pct T/NK"),
    # LuCA spells the selection out in the sub-dataset name
    "Leader_Merad_2021_10x_5p_v1_CD2": ("T_sort", "dataset name CD2; 94.3 pct T/NK"),
    "Leader_Merad_2021_10x_3p_v1_sort": ("CD45pos", "dataset name sort"),
    "Leader_Merad_2021_10x_3p_v2_sort": ("CD45pos", "dataset name sort"),
}


def clean_symbol(name: str) -> str:
    """Strip the Ensembl id suffix from a LuCA feature name of the form SYMBOL_ENSG..."""
    m = re.match(r"^(?P<sym>.+)_(?P<ens>ENSG\d{11})$", name)
    return m.group("sym") if m else name


def rss_gb() -> float:
    """Peak resident memory of this process in GB."""
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 ** 3


def log(msg: str) -> None:
    """Print a message with a time stamp and the peak memory."""
    print(f"[{time.strftime('%H:%M:%S')}] [{rss_gb():5.2f} GB] {msg}", flush=True)


# ---------------------------------------------------------------------------
# 1. obs
# ---------------------------------------------------------------------------

def read_obs(f):
    """Read the LuCA obs columns, keep primary non-Lambrechts cells (asserts the 66,296 / 1,061,555 counts); returns the table and kept row numbers."""
    o = f["obs"]

    def cat(k):
        """Decode one categorical obs column to an array of strings."""
        g = o[k]
        return (np.asarray(g["categories"]).astype(str))[np.asarray(g["codes"])]

    prim = np.asarray(o["is_primary_data"])
    donor = cat("donor_id")
    keep = prim & ~pd.Series(donor).str.startswith(EXCLUDE_DONOR_PREFIX).values

    n_prim, n_keep = int(prim.sum()), int(keep.sum())
    dropped = n_prim - n_keep
    log(f"cells {prim.size:,} | is_primary_data True {n_prim:,} | "
        f"Lambrechts removed {dropped:,} | kept {n_keep:,}")
    assert dropped == EXPECT_DROPPED, f"Lambrechts exclusion removed {dropped}, expected {EXPECT_DROPPED}"
    assert n_keep == EXPECT_CELLS, f"{n_keep} cells kept, expected {EXPECT_CELLS}"

    fields = ["study", "assay", "platform", "dataset", "sample", "disease", "origin",
              "origin_fine", "tissue", "suspension_type", "cell_type", "ann_fine",
              "cell_type_tumor", "sex", "ever_smoker", "uicc_stage", "doublet_status"]
    d = pd.DataFrame({k: cat(k)[keep] for k in fields})
    d["donor_id"] = donor[keep]
    d["barcode"] = np.asarray(o["_index"]).astype(str)[keep]
    d["age"] = np.asarray(o["age"])[keep]
    d["total_counts_atlas"] = np.asarray(o["total_counts"])[keep]
    return d, np.where(keep)[0]


def build_obs(d: pd.DataFrame) -> pd.DataFrame:
    """Map LuCA/Census labels to the harmonized obs schema (compartment, level 2, histology, tissue class, platform, enrichment)."""
    cmap_csv = pd.read_csv(f"{META}/build-plan_compartment_map.csv")
    cmap = dict(zip(cmap_csv["census_cell_type"], cmap_csv["level1_compartment"]))
    mstrict = dict(zip(cmap_csv["census_cell_type"], cmap_csv["macrophage_strict"]))

    present = sorted(set(d["cell_type"]))
    unmapped = [c for c in present if c not in cmap]
    assert not unmapped, f"cell_type labels with no compartment: {unmapped}"
    dupe = cmap_csv["census_cell_type"].duplicated()
    assert not dupe.any(), f"a label maps to two compartments: {cmap_csv[dupe]}"
    log(f"compartment map: {len(present)} labels present, all mapped explicitly")

    unmapped2 = sorted(set(d["ann_fine"]) - set(LEVEL2))
    assert not unmapped2, f"ann_fine labels with no level2: {unmapped2}"
    unmapped_h = sorted(set(d["disease"]) - set(HISTOLOGY))
    assert not unmapped_h, f"disease values with no histology: {unmapped_h}"
    unmapped_p = sorted(set(d["assay"]) - set(PLATFORM))
    assert not unmapped_p, f"assay values with no platform: {unmapped_p}"
    pairs = sorted(set(zip(d["origin"], d["tissue"])))
    unmapped_t = [p for p in pairs if p not in TISSUE_CLASS]
    assert not unmapped_t, f"(origin, tissue) pairs with no tissue_class: {unmapped_t}"

    assert set(d["suspension_type"]) == {"cell"}, set(d["suspension_type"])

    comp = d["cell_type"].map(cmap)
    dataset = d["dataset"]
    ct = pd.crosstab(dataset, comp, normalize="index") * 100
    for c in ["epithelial_malignant", "fibroblast", "mural", "endothelial"]:
        if c not in ct:
            ct[c] = 0.0
    non_immune = ct[["epithelial_malignant", "fibroblast", "mural", "endothelial"]].sum(axis=1)
    enrich, basis = {}, {}
    for ds in ct.index:
        if ds in ENRICHMENT_NAMED:
            enrich[ds], basis[ds] = ENRICHMENT_NAMED[ds]
        elif non_immune[ds] < NON_IMMUNE_SORT_CUT:
            enrich[ds] = "CD45pos"
            basis[ds] = f"non-immune {non_immune[ds]:.1f} pct"
        elif ct.loc[ds, "endothelial"] > ENDOTHELIAL_ENRICH_CUT:
            enrich[ds] = "stromal_sort"
            basis[ds] = f"endothelial {ct.loc[ds, 'endothelial']:.1f} pct"
        else:
            enrich[ds] = "whole_tissue"
            basis[ds] = f"non-immune {non_immune[ds]:.1f} pct"
    enr = dataset.map(enrich)

    hist = d["disease"].map(HISTOLOGY)
    platform = d["assay"].map(PLATFORM)

    obs = pd.DataFrame({
        "cell_id": SOURCE_KEY + "|" + d["sample"] + "|" + d["barcode"],
        "source_key": SOURCE_KEY,
        "study": d["study"],
        "donor_id": d["donor_id"],
        "library_id": d["sample"],
        "histology": hist,
        "histology_source": "census_disease",
        "histology_confidence": np.where(hist == "NSCLC_unspecified", "unknown", "high"),
        "tissue_class": [TISSUE_CLASS[p] for p in zip(d["origin"], d["tissue"])],
        "sampling_site": "NA",
        "platform": platform,
        "modality": "single_cell",
        "enrichment": enr,
        "composition_safe": (enr == "whole_tissue").values,
        "treatment": "unknown",
        "block_id": d["study"] + "_" + platform,
        "level1_compartment": comp,
        "level2_celltype": d["ann_fine"].map(LEVEL2),
        "level3_celltype_native": d["ann_fine"],
        "celltype_source": "census",
        "celltype_confidence": np.nan,
        "is_malignant": pd.array([pd.NA] * len(d), dtype="boolean"),
        "malignancy_score": np.nan,
        "n_counts": 0,          # filled from layers/count while streaming
        "n_genes": 0,
        "counts_are_estimates": False,   # filled from layers/count while streaming
        "doublet_score": np.nan,
        "sex": d["sex"],
        "age": d["age"],
        "stage": d["uicc_stage"],
        "smoking": d["ever_smoker"],
        # ---- extras, clearly named, not part of the 30-column schema ----------
        "dataset_luca": d["dataset"],
        "patient_id": d["donor_id"].str.replace("_NAT$", "", regex=True),
        "census_cell_type": d["cell_type"],
        "celltype_tumor_native": d["cell_type_tumor"],
        "macrophage_strict": d["cell_type"].map(mstrict).astype(bool),
        "tissue_native": d["tissue"],
        "origin_native": d["origin"],
        "n_counts_panel": 0,
        "n_genes_panel": 0,
        "total_counts_atlas": d["total_counts_atlas"],
    })
    obs.index = pd.Index(obs["cell_id"].values, name=None)
    assert obs.index.is_unique, "cell_id is not unique"
    assert obs["donor_id"].notna().all() and (obs["donor_id"] != "").all()
    assert not obs["composition_safe"][obs["enrichment"] != "whole_tissue"].any()

    side = pd.DataFrame({
        "dataset_luca": ct.index,
        "n_cells": dataset.value_counts().reindex(ct.index).values,
        "enrichment": [enrich[i] for i in ct.index],
        "composition_safe": [enrich[i] == "whole_tissue" for i in ct.index],
        "basis": [basis[i] for i in ct.index],
        "non_immune_pct": non_immune.round(2).values,
    })
    for c in sorted(ct.columns):
        side[f"pct_{c}"] = ct[c].round(2).values
    side.sort_values("n_cells", ascending=False).to_csv(
        f"{OUT}/arm_luca_enrichment.csv", index=False)
    return obs


# ---------------------------------------------------------------------------
# 2. var / panel
# ---------------------------------------------------------------------------

def build_var(f):
    """Resolve the panel genes in LuCA var by Ensembl id; returns var, column positions, total genes and the missing genes."""
    v = f["var"]
    gene_id = np.asarray(v["gene_id"]).astype(str)
    cats = np.asarray(v["feature_name"]["categories"]).astype(str)
    feature_name = cats[np.asarray(v["feature_name"]["codes"])]
    var_all = pd.DataFrame({"gene_id": gene_id, "feature_name": feature_name})
    var_all["gene_symbol"] = [clean_symbol(x) for x in feature_name]

    panel = pd.read_csv(f"{OUT}/panel_genes.csv")
    want = panel[panel["ensembl_id"].astype(str) != ""].set_index("ensembl_id")
    hit = var_all["gene_id"].isin(want.index).values
    var = var_all[hit].copy()
    var["source_set"] = want.loc[var["gene_id"], "source_set"].values
    var["is_hvg"] = want.loc[var["gene_id"], "is_hvg"].values
    var["panel_symbol"] = want.loc[var["gene_id"], "gene_symbol"].values
    var.index = pd.Index(var["gene_id"].values, name=None)

    missing = panel[~panel["ensembl_id"].isin(set(var["gene_id"]))]
    log(f"panel {len(panel)} genes | resolved in LuCA var {len(var)} | "
        f"missing {len(missing)}")
    for _, r in missing.iterrows():
        log(f"  missing: {r.gene_symbol} ({r.ensembl_id or 'no ensembl id'}) "
            f"[{r.source_set}]")
    mism = var[var["gene_symbol"] != var["panel_symbol"]]
    assert mism.empty, f"symbol disagreement on an ensembl join: {mism}"
    return var, np.where(hit)[0], len(var_all), missing


# ---------------------------------------------------------------------------
# 3. matrix
# ---------------------------------------------------------------------------

def stream_counts(f, rows, gene_pos, n_genes_all, study_codes, study_cats):
    """Read layers/count block by block; keep `rows` and the panel columns."""
    keep_gene = np.zeros(n_genes_all, dtype=bool)
    keep_gene[gene_pos] = True
    gmap = np.full(n_genes_all, -1, dtype=np.int32)
    gmap[gene_pos] = np.arange(gene_pos.size, dtype=np.int32)

    g = f["layers"]["count"]
    indptr = np.asarray(g["indptr"], dtype=np.int64)
    D, I = g["data"], g["indices"]
    n_rows_all = indptr.size - 1

    row_keep = np.zeros(n_rows_all, dtype=bool)
    row_keep[rows] = True

    data_parts, ind_parts, nnz_parts = [], [], []
    n_counts, n_genes_c, n_counts_p, n_genes_p = [], [], [], []
    nonint_by_study = np.zeros(len(study_cats), dtype=np.int64)
    nnz_by_study = np.zeros(len(study_cats), dtype=np.int64)

    t0 = time.time()
    for r0 in range(0, n_rows_all, ROW_BLOCK):
        r1 = min(r0 + ROW_BLOCK, n_rows_all)
        m = row_keep[r0:r1]
        if not m.any():
            continue
        loc = indptr[r0:r1 + 1] - indptr[r0]
        p0, p1 = int(indptr[r0]), int(indptr[r1])
        dat = np.asarray(D[p0:p1], dtype=np.float32)
        ind = np.asarray(I[p0:p1])

        def rowsum(x, dtype=np.float64):
            """Sum of x within each row of the current row block (via cumulative sums)."""
            cs = np.empty(x.size + 1, dtype=dtype)
            cs[0] = 0
            np.cumsum(x, dtype=dtype, out=cs[1:])
            return cs[loc[1:]] - cs[loc[:-1]]

        nonint = dat != np.floor(dat)
        ni_row = rowsum(nonint.astype(np.int64), np.int64)
        rn = np.diff(loc)
        s = study_codes[r0:r1][m]
        nonint_by_study += np.bincount(s, weights=ni_row[m],
                                       minlength=len(study_cats)).astype(np.int64)
        nnz_by_study += np.bincount(s, weights=rn[m],
                                    minlength=len(study_cats)).astype(np.int64)

        n_counts.append(rowsum(dat)[m])
        n_genes_c.append(rn[m])

        keep = keep_gene[ind]
        rownnz = rowsum(keep.astype(np.int64), np.int64)
        n_genes_p.append(rownnz[m])
        n_counts_p.append(rowsum(np.where(keep, dat, np.float32(0)))[m])

        ptr = np.empty(r1 - r0 + 1, dtype=np.int64)
        ptr[0] = 0
        np.cumsum(rownnz, out=ptr[1:])
        blk = sp.csr_matrix(
            (dat[keep], gmap[ind[keep]], ptr), shape=(r1 - r0, gene_pos.size))
        blk = blk[m]
        data_parts.append(blk.data.astype(np.float32, copy=False))
        ind_parts.append(blk.indices.astype(np.int32, copy=False))
        nnz_parts.append(np.diff(blk.indptr).astype(np.int64))
        del dat, ind, keep, blk
        if (r0 // ROW_BLOCK) % 10 == 0:
            log(f"  rows {r1:,}/{n_rows_all:,}  ({time.time() - t0:.0f}s)")

    gc.collect()
    nnz_row = np.concatenate(nnz_parts)
    del nnz_parts
    total_nnz = int(nnz_row.sum())
    assert total_nnz < 2 ** 31 - 1, f"nnz {total_nnz} needs int64 indptr"
    ptr = np.empty(nnz_row.size + 1, dtype=np.int32)
    ptr[0] = 0
    np.cumsum(nnz_row, out=ptr[1:])
    X = sp.csr_matrix(
        (np.concatenate(data_parts), np.concatenate(ind_parts), ptr),
        shape=(rows.size, gene_pos.size))
    del data_parts, ind_parts
    gc.collect()
    stats = pd.DataFrame({
        "n_counts": np.concatenate(n_counts),
        "n_genes": np.concatenate(n_genes_c),
        "n_counts_panel": np.concatenate(n_counts_p),
        "n_genes_panel": np.concatenate(n_genes_p),
    })
    est = pd.DataFrame({"study": list(study_cats),
                        "nonint": nonint_by_study,
                        "nnz": nnz_by_study})
    est = est[est["nnz"] > 0].copy()
    est["pct"] = 100 * est["nonint"] / est["nnz"]
    return X, stats, est


# ---------------------------------------------------------------------------

def main() -> int:
    """Build the LuCA arm, write arm_luca.h5ad and its tables, and print the report."""
    os.makedirs(OUT, exist_ok=True)
    import anndata as ad

    f = h5py.File(LUCA, "r")
    d, rows = read_obs(f)

    nd = d.groupby("disease", observed=True)["donor_id"].nunique()
    n_lusc = int(nd.get("squamous cell lung carcinoma", 0))
    n_luad = int(nd.get("lung adenocarcinoma", 0))
    log(f"donors: LUSC {n_lusc}, LUAD {n_luad} "
        f"(NSCLC_unspecified {int(nd.get('non-small cell lung carcinoma', 0))}, "
        f"normal {int(nd.get('normal', 0))}, "
        f"COPD {int(nd.get('chronic obstructive pulmonary disease', 0))})")
    assert (n_lusc, n_luad) == (EXPECT_LUSC_DONORS, EXPECT_LUAD_DONORS), \
        f"donor assertion failed: LUSC {n_lusc} LUAD {n_luad}"

    obs = build_obs(d)
    var, gene_pos, n_genes_all, missing = build_var(f)

    sg = f["obs"]["study"]
    study_cats = np.asarray(sg["categories"]).astype(str)
    study_codes = np.asarray(sg["codes"]).astype(np.int64)

    log("streaming layers/count ...")
    X, stats, est = stream_counts(f, rows, gene_pos, n_genes_all,
                                  study_codes, study_cats)
    f.close()
    log(f"X {X.shape} nnz {X.nnz:,} density {X.nnz / (X.shape[0] * X.shape[1]):.4f}")

    obs["n_counts"] = stats["n_counts"].values.astype(np.int64)
    obs["n_genes"] = stats["n_genes"].values.astype(np.int32)
    obs["n_counts_panel"] = stats["n_counts_panel"].values.astype(np.float32)
    obs["n_genes_panel"] = stats["n_genes_panel"].values.astype(np.int32)
    est_studies = set(est.loc[est["nonint"] > 0, "study"])
    obs["counts_are_estimates"] = obs["study"].isin(est_studies).values
    est.sort_values("pct", ascending=False).to_csv(
        f"{OUT}/arm_luca_counts_integrality.csv", index=False)
    log("non-integer entries in layers/count, by study:")
    for _, r in est[est["nonint"] > 0].sort_values("pct", ascending=False).iterrows():
        log(f"  {r.study:24s} {r.nonint:12,d} / {r.nnz:12,d}  ({r.pct:.2f} pct)")

    for c in ["source_key", "study", "histology", "histology_source",
              "histology_confidence", "tissue_class", "sampling_site", "platform",
              "modality", "enrichment", "treatment", "block_id",
              "level1_compartment", "level2_celltype", "celltype_source", "sex",
              "dataset_luca", "census_cell_type", "celltype_tumor_native",
              "level3_celltype_native", "tissue_native", "origin_native",
              "library_id"]:
        obs[c] = pd.Categorical(obs[c])

    a = ad.AnnData(X=X, obs=obs, var=var)
    a.uns["build"] = {
        "stage": "1a",
        "arm": "luca",
        "source": LUCA,
        "source_title": "The single-cell lung cancer atlas (LuCA) -- extended atlas",
        "census_dataset_id": "1e6a6ef9-7ec9-4c90-bbfb-2ad3c3165fd1",
        "X_layer": "layers/count of the source object, raw counts, not normalised",
        "filter": "is_primary_data == True and donor_id does not start with Lambrechts",
        "panel": f"{OUT}/panel_genes.csv",
        "compartment_map": f"{META}/build-plan_compartment_map.csv",
        "var_index": "ENSEMBL gene id",
        "built": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    log(f"writing {H5AD} ...")
    a.write_h5ad(H5AD)
    log(f"written, {os.path.getsize(H5AD) / 1024 ** 3:.2f} GB")

    # ---- report -----------------------------------------------------------
    pd.set_option("display.width", 250)
    pd.set_option("display.max_rows", 400)
    o = a.obs
    print("\n=== shape ===")
    print(f"{a.shape[0]:,} cells x {a.shape[1]:,} genes")

    print("\n=== donors and cells by histology ===")
    t = o.groupby("histology", observed=True).agg(
        donors=("donor_id", "nunique"), patients=("patient_id", "nunique"),
        libraries=("library_id", "nunique"), cells=("histology", "size"))
    print(t.to_string())

    print("\n=== cells by compartment ===")
    print(pd.crosstab(o["level1_compartment"], o["histology"],
                      margins=True, margins_name="all").to_string())

    print("\n=== per-block donor counts ===")
    b = o.groupby("block_id", observed=True).apply(
        lambda g: pd.Series({
            "cells": len(g),
            "donors": g.donor_id.nunique(),
            "LUSC": g.loc[g.histology == "LUSC", "donor_id"].nunique(),
            "LUAD": g.loc[g.histology == "LUAD", "donor_id"].nunique(),
            "NSCLC_unspecified": g.loc[g.histology == "NSCLC_unspecified", "donor_id"].nunique(),
            "normal": g.loc[g.histology == "normal", "donor_id"].nunique(),
            "COPD": g.loc[g.histology == "COPD", "donor_id"].nunique(),
            "composition_safe_cells": int(g.composition_safe.sum()),
        }), include_groups=False).sort_values("cells", ascending=False)
    print(b.to_string())
    b.to_csv(f"{OUT}/arm_luca_block_table.csv")

    print("\n=== tissue_class ===")
    print(o["tissue_class"].value_counts().to_string())
    print("\n=== enrichment ===")
    print(o.groupby("enrichment", observed=True).agg(
        cells=("enrichment", "size"), donors=("donor_id", "nunique")).to_string())
    print("\n=== level2_celltype ===")
    print(o["level2_celltype"].value_counts().to_string())
    print("\n=== depth ===")
    print(o.groupby("platform", observed=True).agg(
        n_counts_median=("n_counts", "median"), n_genes_median=("n_genes", "median"),
        panel_counts_median=("n_counts_panel", "median"),
        panel_genes_median=("n_genes_panel", "median"), cells=("n_counts", "size")).to_string())

    don = o.groupby(["study", "block_id", "donor_id", "histology"], observed=True).agg(
        cells=("donor_id", "size"), libraries=("library_id", "nunique"),
        enrichment=("enrichment", lambda s: ",".join(sorted(set(map(str, s))))),
        tissue_class=("tissue_class", lambda s: ",".join(sorted(set(map(str, s))))))
    comp = pd.crosstab(o["donor_id"], o["level1_compartment"])
    comp.columns = [str(c) for c in comp.columns]      # a CategoricalIndex cannot be joined on
    don = don.join(comp, on="donor_id")
    don.reset_index().to_csv(f"{OUT}/arm_luca_donor_table.csv", index=False)
    lab = o.groupby(["level3_celltype_native", "census_cell_type", "level1_compartment",
                     "level2_celltype"], observed=True).size().rename("cells").reset_index()
    lab.to_csv(f"{OUT}/arm_luca_label_map.csv", index=False)
    print(f"\npeak RSS {rss_gb():.2f} GB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
