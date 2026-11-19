#!/usr/bin/env python
"""
Atlas stage 01, step arm_emtab6653: the E-MTAB-6653 arm (Lambrechts / Thienpont 2018, validation
cohort, patients 6/7/8) at the panel gene axis.

Paper items: part of the atlas behind Figure 2, Figure 3, Supplementary Figures 2 to 4 and
Supplementary Data 1 (the Lambrechts copies are included once, here and in E-MTAB-6149).
No random numbers are drawn.

Builds output/RNA/meta_dataset/stage1a/arm_emtab6653.h5ad to the 30-column obs schema in
output/RNA/meta_dataset/build-plan_obs_schema.csv.

Source: EBI Single Cell Expression Atlas aggregated Matrix Market export of E-MTAB-6653.
  23,979 genes (Ensembl ENSG, unversioned) x 32,341 cells, values are Alevin EM estimates and
  are NOT integer UMI counts -> counts_are_estimates is True for every cell.

Everything that is not derivable from the deposited files is either taken from a spec file in
output/RNA/meta_dataset/ or from Supplementary table 1 of the paper, and the origin is recorded
in .uns['provenance'].  No cell type labels are invented: the EBI deposit ships k-means cluster
integers and unnamed marker tables only, so level1/level2/level3 are left empty and
celltype_source is 'transferred'.

Read-only inputs.  Writes only into output/RNA/meta_dataset/stage1a/.
Run: python atlas/steps/stage1a_arm_emtab6653.py   (normally through atlas/01_arms_luca_lambrechts.py)
"""

import json
import os
from datetime import datetime, timezone

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

import os as _os  # noqa: E402
import sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # atlas/
from common import ROOT, SCRNA, SCRATCH, DATA, RSCRIPT  # noqa: E402,F401  paths from config.json
SRC = f"{SCRNA}/E-MTAB-6653"
MM = f"{SRC}/raw_counts/matrix_market"
META = f"{SRC}/metadata"
OUTDIR = f"{ROOT}/output/RNA/meta_dataset/stage1a"
PANEL = f"{OUTDIR}/panel_genes.csv"
INVENTORY = f"{ROOT}/output/RNA/meta_dataset/lambrechts.csv"
OUT_H5AD = f"{OUTDIR}/arm_emtab6653.h5ad"

SOURCE_KEY = "EMTAB6653"
STUDY = "Lambrechts_Thienpont_2018"
BLOCK_ID = "EMTAB6653_10xv2"          # verbatim block_id of this arm's row in build-plan_blocks.csv
PLATFORM = "10x 3' v2"                # SDRF Comment[library construction] = 10xV2, polyA RNA
MODALITY = "single_cell"
ENRICHMENT = "whole_tissue"           # no sort anywhere in the SDRF/IDF protocols
TREATMENT = "naive"                   # Suppl. table 1: "primary, previously-untreated lesions"

os.makedirs(OUTDIR, exist_ok=True)
log = []


def say(msg):
    """Print a message and keep it in the run log."""
    print(msg, flush=True)
    log.append(msg)


def check(name, ok, detail=""):
    """Print a PASS or FAIL line and stop the script when the check fails."""
    say(f"[{'PASS' if ok else 'FAIL'}] {name}{(' :: ' + detail) if detail else ''}")
    if not ok:
        raise SystemExit(f"ASSERTION FAILED: {name} :: {detail}")


# ---------------------------------------------------------------- 1. metadata: SDRF
sdrf = pd.read_csv(f"{META}/E-MTAB-6653.sdrf.txt", sep="\t", dtype=str)
check("SDRF has 12 rows", len(sdrf) == 12, f"n={len(sdrf)}")

# EXPLICIT maps only.  Never substring: 'tumor middle (in between core and edge sample)'
# contains the substring 'core', and 'normal tissue adjacent to tumor' contains 'tumor'.
SITE_MAP = {
    "tumour core": "core",
    "tumour middle (in between core and edge sample)": "middle",
    "tumour edge": "edge",
    "normal tissue adjacent to tumour": "normal",
}
TISSUE_CLASS_MAP = {
    "core": "primary_tumour",
    "middle": "primary_tumour",
    "edge": "primary_tumour",
    "normal": "adjacent_normal",
}
# Histology is NOT in the SDRF (disease == 'lung carcinoma' for all 12 samples) nor in the IDF
# nor in the EBI experiment-design file.  It comes from Supplementary table 1, page 21 of
# E-MTAB-6149/source_support/Lambrechts_2018_supplement.pdf:  6 Adenomatous, 7 Squamous, 8 Pleiomorphic.
HISTOLOGY_MAP = {"6": "LUAD", "7": "LUSC", "8": "pleomorphic"}
HISTOLOGY_PAPER_WORD = {"6": "Adenomatous", "7": "Squamous", "8": "Pleiomorphic"}
STAGE_MAP = {"6": "IIIA", "7": "IB", "8": "IIB"}
SMOKING_MAP = {"6": "former", "7": "former", "8": "active"}
TNM_MAP = {"6": "pT4N1M0", "7": "pT2aN0M0", "8": "pT3N0M0"}
COPD_MAP = {"6": "yes_GOLD2", "7": "no", "8": "no"}

unmapped = sorted(set(sdrf["Characteristics[sampling site]"]) - set(SITE_MAP))
check("every SDRF sampling site is in the explicit map", not unmapped, str(unmapped))
check("SDRF disease is uninformative for histology",
      set(sdrf["Characteristics[disease]"]) == {"lung carcinoma"},
      str(sorted(set(sdrf["Characteristics[disease]"]))))
check("SDRF individuals are 6/7/8", sorted(set(sdrf["Characteristics[individual]"])) == ["6", "7", "8"])
check("SDRF chemistry is 10xV2 throughout", set(sdrf["Comment[library construction]"]) == {"10xV2"})

lib = pd.DataFrame({
    "library_id": sdrf["Comment[ENA_RUN]"].values,             # ERR27440xx, the key the matrix uses
    "library_name_native": sdrf["Source Name"].values,         # 'Sample 1a' .. 'Sample 3d'
    "patient_native": sdrf["Characteristics[individual]"].values,
    "sampling_site": [SITE_MAP[s] for s in sdrf["Characteristics[sampling site]"]],
    "sex": sdrf["Characteristics[sex]"].values,
    "age": sdrf["Characteristics[age]"].astype(float).values,
}).set_index("library_id")
lib["tissue_class"] = lib["sampling_site"].map(TISSUE_CLASS_MAP)
lib["donor_id"] = [f"{STUDY}_{p}" for p in lib["patient_native"]]
lib["histology"] = lib["patient_native"].map(HISTOLOGY_MAP)
check("12 unique libraries", lib.index.is_unique and len(lib) == 12)
check("3 donors x 4 sites", (lib.groupby(["donor_id", "sampling_site"]).size() == 1).all())

# sex/age cross-check between the SDRF and the paper supplement
PAPER_SEX_AGE = {"6": ("male", 65.0), "7": ("male", 60.0), "8": ("female", 55.0)}
for p, (sx, ag) in PAPER_SEX_AGE.items():
    sub = lib[lib["patient_native"] == p]
    check(f"patient {p} sex/age agree SDRF vs supplement",
          set(sub["sex"]) == {sx} and set(sub["age"]) == {ag},
          f"SDRF={sorted(set(sub['sex']))}/{sorted(set(sub['age']))} supplement={sx}/{ag}")

# ---------------------------------------------------------------- 2. matrix
say("\nloading Matrix Market ...")
coo = scipy.io.mmread(f"{MM}/E-MTAB-6653.aggregated_filtered_counts.mtx")
genes_raw = pd.read_csv(f"{MM}/E-MTAB-6653.aggregated_filtered_counts.mtx_rows",
                        sep="\t", header=None, dtype=str)
cells_raw = pd.read_csv(f"{MM}/E-MTAB-6653.aggregated_filtered_counts.mtx_cols",
                        sep="\t", header=None, dtype=str)
gene_ids = genes_raw[0].values
cell_ids_native = cells_raw[0].values
say(f"  mtx {coo.shape[0]} genes x {coo.shape[1]} cells, {coo.nnz:,} nonzeros, dtype {coo.data.dtype}")
check("shape is 23,979 genes x 32,341 cells", coo.shape == (23979, 32341), str(coo.shape))
check("rows file matches gene axis", len(gene_ids) == coo.shape[0])
check("cols file matches cell axis", len(cell_ids_native) == coo.shape[1])
check("gene ids unique", pd.Index(gene_ids).is_unique)
check("cell ids unique", pd.Index(cell_ids_native).is_unique)
check("gene ids are unversioned ENSG", bool(pd.Series(gene_ids).str.match(r"^ENSG\d{11}$").all()))
check("the two columns of mtx_rows are identical", bool((genes_raw[0] == genes_raw[1]).all()))

# --------------------------------------------- 2a. THE FRACTIONAL-ENTRY CHECK (task item 2)
data = coo.data
n_nz = data.size
frac_exact = int(np.sum(data != np.rint(data)))
resid = np.abs(data - np.rint(data))
frac_1e9 = int(np.sum(resid > 1e-9))
frac_1e6 = int(np.sum(resid > 1e-6))
# the figure quoted for the deposit was computed with np.isclose(d, round(d)), whose default
# rtol=1e-5 is RELATIVE, so an entry of 1000.005 is scored as integer. Reproduced here for
# comparison, but the exact test is the one that governs counts_are_estimates.
frac_isclose = int(np.sum(~np.isclose(data, np.rint(data))))
below_one = int(np.sum(data < 1.0))
say(f"\n  nonzero entries                    {n_nz:,}")
say(f"  fractional (exact != rint)         {frac_exact:,}  = {100*frac_exact/n_nz:.4f} pct")
say(f"  fractional (|resid| > 1e-9)        {frac_1e9:,}  = {100*frac_1e9/n_nz:.4f} pct")
say(f"  fractional (|resid| > 1e-6)        {frac_1e6:,}  = {100*frac_1e6/n_nz:.4f} pct")
say(f"  fractional (np.isclose, rel 1e-5)  {frac_isclose:,}  = {100*frac_isclose/n_nz:.4f} pct")
say(f"  entries < 1.0                      {below_one:,}  = {100*below_one/n_nz:.4f} pct")
say(f"  min {data.min():.3g}  max {data.max():.6g}  total signal {data.sum():,.3f}")
PCT_FRAC = 100.0 * frac_exact / n_nz
PCT_FRAC_ISCLOSE = 100.0 * frac_isclose / n_nz
check("matrix is NOT integer, i.e. counts_are_estimates is justified", frac_exact > 0,
      f"{frac_exact:,} fractional entries")
check("exact fractional share is in the reported 15-16 pct band",
      15.0 < PCT_FRAC < 16.5, f"measured {PCT_FRAC:.4f} pct")
check("np.isclose definition reproduces the 15.91 pct inventory figure",
      abs(PCT_FRAC_ISCLOSE - 15.91) < 0.01, f"measured {PCT_FRAC_ISCLOSE:.4f} pct")

# --------------------------------------------- 2b. per-cell depth over the FULL gene space
csc = coo.tocsc()
del coo
n_counts_exact = np.asarray(csc.sum(axis=0)).ravel()
n_genes = np.diff(csc.indptr).astype(np.int32)
cs_resid = np.abs(n_counts_exact - np.rint(n_counts_exact))
cellsum_frac = int(np.sum(cs_resid > 0))
cellsum_frac_1e3 = int(np.sum(cs_resid > 1e-3))
say(f"  per-cell totals not exactly integer   {cellsum_frac:,} of {csc.shape[1]:,} cells "
    f"({100*cellsum_frac/csc.shape[1]:.2f} pct), median residual {np.median(cs_resid):.2g}")
say(f"  per-cell totals off by more than 1e-3 {cellsum_frac_1e3:,} of {csc.shape[1]:,} cells "
    f"({100*cellsum_frac_1e3/csc.shape[1]:.2f} pct)  <- the fractions mostly redistribute "
    f"within a cell")
check("n_genes positive for every cell", bool((n_genes > 0).all()))
check("n_counts positive for every cell", bool((n_counts_exact > 0).all()))

# ---------------------------------------------------------------- 3. obs
run = np.array([c.split("-", 1)[0] for c in cell_ids_native])
barcode = np.array([c.split("-", 1)[1] for c in cell_ids_native])
check("every run prefix is an SDRF ENA run", set(run) <= set(lib.index), str(sorted(set(run) - set(lib.index))))
check("all 12 libraries carry cells", len(set(run)) == 12, f"{len(set(run))} present")
check("barcodes are 14nt", bool((pd.Series(barcode).str.len() == 14).all()))

meta = lib.loc[run].reset_index()
cell_id = np.array([f"{SOURCE_KEY}|{r}|{b}" for r, b in zip(run, barcode)])
check("cell_id globally unique", pd.Index(cell_id).is_unique)
check("cell_id form <source_key>|<library_id>|<barcode>",
      bool(pd.Series(cell_id).str.match(rf"^{SOURCE_KEY}\|ERR\d+\|[ACGT]{{14}}$").all()))

n = len(cell_id)
EMPTY = pd.Categorical([""] * n, categories=[""])
obs = pd.DataFrame(index=pd.Index(cell_id, name="cell_id"))
obs["source_key"] = pd.Categorical([SOURCE_KEY] * n)
obs["study"] = pd.Categorical([STUDY] * n)
obs["donor_id"] = meta["donor_id"].values                      # <study>_<native donor label>
obs["library_id"] = meta["library_id"].values                  # ENA run, the native matrix key
obs["histology"] = pd.Categorical(meta["histology"].values,
                                  categories=["LUAD", "LUSC", "pleomorphic"])
obs["histology_source"] = pd.Categorical(["paper_supplement"] * n)
obs["histology_confidence"] = pd.Categorical(["high"] * n)
obs["tissue_class"] = pd.Categorical(meta["tissue_class"].values,
                                     categories=["primary_tumour", "adjacent_normal"])
obs["sampling_site"] = pd.Categorical(meta["sampling_site"].values,
                                      categories=["core", "middle", "edge", "normal"])
obs["platform"] = pd.Categorical([PLATFORM] * n)
obs["modality"] = pd.Categorical([MODALITY] * n)
obs["enrichment"] = pd.Categorical([ENRICHMENT] * n)
obs["composition_safe"] = (obs["enrichment"].astype(str) == "whole_tissue").values
obs["treatment"] = pd.Categorical([TREATMENT] * n)
obs["block_id"] = pd.Categorical([BLOCK_ID] * n)
# No cell type labels are deposited.  See uns['celltype_decision'].
obs["level1_compartment"] = EMPTY
obs["level2_celltype"] = EMPTY
obs["level3_celltype_native"] = EMPTY
obs["celltype_source"] = pd.Categorical(["transferred"] * n)
obs["celltype_confidence"] = np.full(n, np.nan, dtype=np.float32)
obs["is_malignant"] = pd.array([pd.NA] * n, dtype="boolean")
obs["malignancy_score"] = np.full(n, np.nan, dtype=np.float32)
obs["n_counts"] = np.rint(n_counts_exact).astype(np.int64)
obs["n_genes"] = n_genes
obs["counts_are_estimates"] = np.ones(n, dtype=bool)
obs["doublet_score"] = np.full(n, np.nan, dtype=np.float32)
obs["sex"] = pd.Categorical(meta["sex"].values, categories=["male", "female", "unknown"])
obs["age"] = meta["age"].values.astype(np.float32)
obs["stage"] = pd.Categorical(meta["patient_native"].map(STAGE_MAP).values)
obs["smoking"] = pd.Categorical(meta["patient_native"].map(SMOKING_MAP).values)

# extra provenance columns, outside the 30-column schema but lossless
obs["barcode_native"] = barcode
obs["cell_id_native"] = cell_ids_native
obs["library_name_native"] = pd.Categorical(meta["library_name_native"].values)
obs["patient_native"] = pd.Categorical(meta["patient_native"].values)
obs["accession"] = pd.Categorical(["E-MTAB-6653"] * n)
obs["histology_paper_word"] = pd.Categorical(meta["patient_native"].map(HISTOLOGY_PAPER_WORD).values)
obs["tnm"] = pd.Categorical(meta["patient_native"].map(TNM_MAP).values)
obs["copd"] = pd.Categorical(meta["patient_native"].map(COPD_MAP).values)
obs["n_counts_exact"] = n_counts_exact.astype(np.float64)

# ---------------------------------------------------------------- 4. schema conformance
schema = pd.read_csv(f"{DATA}/build-plan_obs_schema.csv")  # hand-made table, shipped in atlas/data/
schema = schema.dropna(subset=["column"])
required = schema.loc[schema["required"] == "yes", "column"].tolist()
missing = [c for c in schema["column"] if c not in obs.columns and c != "cell_id"]
check("every schema column present", not missing, str(missing))
check("donor_id never empty", bool((obs["donor_id"].astype(str).str.len() > 0).all()))
check("platform never empty", bool((obs["platform"].astype(str).str.len() > 0).all()))
check("modality never empty", bool((obs["modality"].astype(str).str.len() > 0).all()))
check("counts_are_estimates True everywhere", bool(obs["counts_are_estimates"].all()))
check("composition_safe True iff whole_tissue",
      bool((obs["composition_safe"] == (obs["enrichment"].astype(str) == "whole_tissue")).all()))
for c in required:
    if c in ("level1_compartment", "level2_celltype", "level3_celltype_native", "cell_id"):
        continue
    check(f"required column {c} has no missing value", bool(obs[c].notna().all()))

# per-site-per-patient counts must reproduce the independent inventory table
inv = pd.read_csv(INVENTORY)
inv = inv[inv["accession"] == "E-MTAB-6653"]
got = obs.groupby(["patient_native", "sampling_site"], observed=True).size()
for _, r in inv.iterrows():
    check(f"patient {r['patient']} {r['site']} = {r['cells_local']} cells",
          int(got.loc[(str(r["patient"]), r["site"])]) == int(r["cells_local"]),
          f"got {int(got.loc[(str(r['patient']), r['site'])])}")
for _, r in inv.iterrows():
    check(f"patient {r['patient']} histology matches inventory",
          HISTOLOGY_MAP[str(r["patient"])] == r["histology_class"])

# ---------------------------------------------------------------- 5. panel subset
panel = pd.read_csv(PANEL)
check("panel has 2,059 genes", len(panel) == 2059, f"n={len(panel)}")
check("panel ensembl ids unique", panel["ensembl_id"].is_unique)
gene_pos = pd.Series(np.arange(len(gene_ids)), index=gene_ids)
panel["row"] = panel["ensembl_id"].map(gene_pos)
present = panel[panel["row"].notna()].copy()
present["row"] = present["row"].astype(int)
say(f"\n  panel genes found in this source: {len(present)} of {len(panel)}")

presence = pd.read_csv(f"{OUTDIR}/panel_genes_source_presence.csv")
expected = int(presence["in_EMTAB6653"].sum())
check("panel-in-source count matches panel_genes_source_presence.csv",
      len(present) == expected, f"got {len(present)}, expected {expected}")
ids_here = set(present["ensembl_id"])
ids_flagged = set(presence.loc[presence["in_EMTAB6653"], "ensembl_id"])
check("the same gene ids, not just the same count", ids_here == ids_flagged,
      f"sym-diff {len(ids_here ^ ids_flagged)}")

X = csc[present["row"].values, :].T.tocsr().astype(np.float32)
check("panel subset shape", X.shape == (n, len(present)), str(X.shape))
say(f"  panel matrix {X.shape[0]:,} cells x {X.shape[1]:,} genes, {X.nnz:,} nonzeros")
frac_panel = int(np.sum(X.data != np.rint(X.data)))
say(f"  fractional entries inside the panel subset: {frac_panel:,} = "
    f"{100*frac_panel/X.nnz:.4f} pct of {X.nnz:,}")

var = pd.DataFrame({
    "gene_symbol": present["gene_symbol"].values,
    "ensembl_id": present["ensembl_id"].values,
    "source_set": present["source_set"].values,
    "is_hvg": present["is_hvg"].values,
    "row_in_source_matrix": present["row"].values,
}).set_index(pd.Index(present["gene_symbol"].values, name="gene_symbol"))
check("var index unique", var.index.is_unique)

obs["n_counts_panel"] = np.asarray(X.sum(axis=1)).ravel().astype(np.float64)
obs["n_genes_panel"] = np.diff(X.indptr).astype(np.int32)

# ---------------------------------------------------------------- 6. write
missing_panel = panel[panel["row"].isna()]
adata = ad.AnnData(X=X, obs=obs, var=var)
adata.uns["arm"] = SOURCE_KEY
adata.uns["accession"] = "E-MTAB-6653"
adata.uns["built_utc"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
adata.uns["script"] = os.path.abspath(__file__)
adata.uns["X_is"] = ("EBI Single Cell Expression Atlas aggregated_filtered_counts, Alevin EM "
                     "abundance estimates, NOT integer UMI counts, unmodified, subset to the "
                     "Stage 1a panel genes present in this source")
adata.uns["fractional_entries"] = {
    "nonzeros_full_matrix": int(n_nz),
    "fractional_full_matrix": int(frac_exact),
    "pct_fractional_full_matrix": float(PCT_FRAC),
    "test": "data != np.rint(data), exact",
    "pct_fractional_np_isclose_rel1e5": float(PCT_FRAC_ISCLOSE),
    "note_on_15.91": ("lambrechts.md quotes 15.91 pct, computed with np.isclose whose default "
                      "rtol=1e-5 is relative, so entries like 1000.005 scored as integer. The "
                      "exact share is 15.98 pct. Both are reproduced by this script."),
    "pct_entries_below_one": float(100 * below_one / n_nz),
    "cells_with_non_integer_total_exact": int(cellsum_frac),
    "cells_with_non_integer_total_tol_1e-3": int(cellsum_frac_1e3),
    "nonzeros_panel_subset": int(X.nnz),
    "pct_fractional_panel_subset": float(100 * frac_panel / X.nnz),
    "rounding_applied": "none",
}
adata.uns["celltype_decision"] = (
    "No cell type labels are deposited for E-MTAB-6653. E-MTAB-6653_clusters.tsv holds 9 k-means "
    "partitions (K=13,22,26,32,40,53,60,75,82; sel.K=TRUE on K=40) of unnamed integer clusters, and "
    "E-MTAB-6653_marker_genes.zip holds one marker table per K in Ensembl ids with no cell type "
    "names anywhere. The deposited partitions are also not donor corrected: at K=13 all 1,219 cells "
    "of the malignant cluster come from patient 7 and patients 6 and 8 contribute none. Therefore "
    "level1_compartment, level2_celltype and level3_celltype_native are empty, celltype_source is "
    "'transferred' and celltype_confidence is NaN. Labels must come from a donor-aware transfer in "
    "a later stage. The deposited cluster integers are carried in obsm['ebi_kmeans_clusters'] as "
    "provenance only and must not be read as cell types.")
adata.uns["histology_provenance"] = (
    "NOT from the SDRF. E-MTAB-6653.sdrf.txt records Characteristics[disease]='lung carcinoma' for "
    "all 12 samples and the words squamous/adenocarcinoma/pleomorphic appear in neither the SDRF, "
    "the IDF nor the EBI experiment-design file. Histology is Supplementary table 1, page 21 of "
    "original_data/external_scRNA/E-MTAB-6149/source_support/Lambrechts_2018_supplement.pdf: "
    "patient 6 Adenomatous, patient 7 Squamous, patient 8 Pleiomorphic. Sex and age in that table "
    "agree with the SDRF for all three patients. Note the Census/LuCA call patient 8 "
    "'non-small cell lung carcinoma', a coarser label, not a contradiction.")
adata.uns["treatment_provenance"] = (
    "Suppl. table 1 legend: 'All tumours were primary, previously-untreated lesions of the "
    "non-small-cell type, resected for curative intent.' Also 'naive' in build-plan_blocks.csv.")
adata.uns["block_id_note"] = (
    "block_id is the verbatim block_id of this arm's row in build-plan_blocks.csv "
    "(EMTAB6653_10xv2). build-plan_donor_table.csv writes its block ids as '<study> | <platform>' "
    "instead, so the merge step must harmonise the two forms once, globally. E-MTAB-6653 must stay "
    "a separate block from E-MTAB-6149 even though both are 10x 3' v2 and one paper: median depth "
    "differs 2.2-fold between the two accessions.")
adata.uns["composition_note"] = (
    "enrichment is whole_tissue and composition_safe is True, which refers to the dissociation "
    "protocol: nothing was sorted. It does not license the deposited k-means partitions for "
    "composition, which are donor confounded.")
adata.uns["panel_genes_absent"] = sorted(missing_panel["gene_symbol"].tolist())
adata.uns["source_files"] = [
    f"{MM}/E-MTAB-6653.aggregated_filtered_counts.mtx",
    f"{MM}/E-MTAB-6653.aggregated_filtered_counts.mtx_rows",
    f"{MM}/E-MTAB-6653.aggregated_filtered_counts.mtx_cols",
    f"{META}/E-MTAB-6653.sdrf.txt",
    f"{META}/E-MTAB-6653.idf.txt",
    f"{META}/E-MTAB-6653_clusters.tsv",
    f"{SCRNA}/E-MTAB-6149/source_support/Lambrechts_2018_supplement.pdf",
]
adata.uns["full_source_shape"] = {"genes": 23979, "cells": 32341, "nonzeros": int(n_nz)}

# deposited k-means partitions, kept as provenance, never as labels
clu = pd.read_csv(f"{META}/E-MTAB-6653_clusters.tsv", sep="\t")
sel_k = int(clu.loc[clu["sel.K"].astype(str).str.upper() == "TRUE", "K"].iloc[0])
kmat = clu.set_index("K").drop(columns=["sel.K"])
kmat = kmat[cell_ids_native]                       # align to matrix column order
adata.obsm["ebi_kmeans_clusters"] = kmat.T.to_numpy(dtype=np.int16)
adata.uns["ebi_kmeans_K"] = [int(k) for k in kmat.index]
adata.uns["ebi_kmeans_sel_K"] = sel_k
check("kmeans block aligned to all cells", adata.obsm["ebi_kmeans_clusters"].shape == (n, len(kmat)))

adata.write_h5ad(OUT_H5AD, compression="gzip")
say(f"\nwrote {OUT_H5AD}  ({os.path.getsize(OUT_H5AD)/1e6:.1f} MB)")

# ---------------------------------------------------------------- 7. read back and report
rb = ad.read_h5ad(OUT_H5AD)
check("round trip shape", rb.shape == adata.shape, str(rb.shape))
check("round trip cell_id unique", rb.obs_names.is_unique)
say(f"\nFINAL shape: {rb.shape[0]:,} cells x {rb.shape[1]:,} genes")

say("\ndonors by histology:")
dh = rb.obs.groupby("donor_id", observed=True)["histology"].agg(lambda s: s.iloc[0])
for h, grp in dh.groupby(dh, observed=True):
    cells = int(rb.obs["histology"].astype(str).eq(h).sum())
    say(f"  {h:<12} {len(grp)} donor(s): {', '.join(sorted(grp.index))}   {cells:,} cells")

say("\ncells per site per patient:")
tab = pd.crosstab(rb.obs["donor_id"], rb.obs["sampling_site"])
tab = tab[["core", "middle", "edge", "normal"]]
tab["total"] = tab.sum(axis=1)
say(tab.to_string())
say("\ncells per tissue_class: " + str(rb.obs["tissue_class"].value_counts().to_dict()))
say(f"\nmedian n_counts {np.median(rb.obs['n_counts']):.0f}, median n_genes {np.median(rb.obs['n_genes']):.0f}")
say(f"measured fractional nonzero entries: {PCT_FRAC:.4f} pct ({frac_exact:,} of {n_nz:,})")

with open(f"{OUTDIR}/arm_emtab6653_buildlog.txt", "w") as fh:
    fh.write("\n".join(log) + "\n")
print("\nDONE")
