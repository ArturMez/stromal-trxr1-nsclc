#!/usr/bin/env python
"""
Atlas stage 01, step merge_panel_object: concatenate arm_luca + arm_emtab6149 + arm_emtab6653 into
panel_object.h5ad (1,146,594 cells x 1,906 genes, 19 blocks).

Paper items: the three-arm panel object of the atlas behind Figure 2, Figure 3, Supplementary Figures
2 to 4 and Supplementary Data 1. No random numbers are drawn.
Reads stage1a/arm_*.h5ad and panel_genes.csv; writes stage1a/panel_object.h5ad,
panel_object_summary.csv and panel_object_gene_intersection.csv.
Run: python atlas/steps/stage1a_merge_panel_object.py   (normally through atlas/01_arms_luca_lambrechts.py)

Design decisions, all of them explicit and all of them reported:
  - Gene join key is gene_symbol, NOT the Ensembl id.  arm_emtab6149 documented 26 panel
    genes whose Ensembl id drifted between the 2018 featuredata build and the LuCA
    reference (TXNIP and MRC1 among them).  An Ensembl join would silently drop them.
  - Intersection is over genes REALLY MEASURED in each arm.  The 65 structural-zero
    columns written into arm_emtab6149 (present_in_source == False) are treated as
    absent, because an all-zero column is a fabricated non-detection.
  - Every categorical remap in this script is an explicit dict with a completeness
    assertion.  No substring matching anywhere.
"""
import sys, json, time
import numpy as np
import pandas as pd
import scipy.sparse as sp
import h5py
import anndata as ad
from anndata.io import read_elem

T0 = time.time()
import os as _os  # noqa: E402
import sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # atlas/
from common import MD  # noqa: E402  paths from config.json, see atlas/common.py
BASE = MD + "/"
S1A = BASE + "stage1a/"
OUT_H5 = S1A + "panel_object.h5ad"
OUT_SUM = S1A + "panel_object_summary.csv"
OUT_GENES = S1A + "panel_object_gene_intersection.csv"

def log(*a):
    """Print a message with the elapsed time."""
    print(f"[{time.time()-T0:7.1f}s]", *a, flush=True)

def die(msg):
    """Print a stop message and exit with status 1."""
    print("\n*** STOP: " + msg, flush=True)
    sys.exit(1)

SCHEMA_COLS = [
    "source_key", "study", "donor_id", "library_id", "histology", "histology_source",
    "histology_confidence", "tissue_class", "sampling_site", "platform", "modality",
    "enrichment", "composition_safe", "treatment", "block_id", "level1_compartment",
    "level2_celltype", "level3_celltype_native", "celltype_source", "celltype_confidence",
    "is_malignant", "malignancy_score", "n_counts", "n_genes", "counts_are_estimates",
    "doublet_score", "sex", "age", "stage", "smoking",
]
REQUIRED = [
    "source_key", "study", "donor_id", "library_id", "histology", "histology_source",
    "histology_confidence", "tissue_class", "platform", "modality", "enrichment",
    "composition_safe", "treatment", "block_id", "level1_compartment", "level2_celltype",
    "level3_celltype_native", "celltype_source", "n_counts", "n_genes",
    "counts_are_estimates",
]
EXTRAS = [
    "patient_id", "donor_id_native", "block_id_native", "level2_celltype_native",
    "technical_replicate_group", "macrophage_strict", "tnm",
    "n_counts_intersect", "n_genes_intersect",
]

# ---------------------------------------------------------------- 1. gene axis
log("reading var of the three arms")
with h5py.File(S1A + "arm_luca.h5ad", "r") as h:
    var_l = read_elem(h["var"])
with h5py.File(S1A + "arm_emtab6149.h5ad", "r") as h:
    var_a = read_elem(h["var"])
with h5py.File(S1A + "arm_emtab6653.h5ad", "r") as h:
    var_b = read_elem(h["var"])

sym_l = pd.Index(var_l["panel_symbol"].astype(str))
present_a = var_a["present_in_source"].astype(bool).values
sym_a = pd.Index(var_a.index.astype(str))
sym_a_present = sym_a[present_a]
sym_b = pd.Index(var_b.index.astype(str))

for nm, ix in [("luca", sym_l), ("6149", sym_a_present), ("6653", sym_b)]:
    if ix.duplicated().any():
        die(f"duplicate gene symbols in {nm}")

panel = pd.read_csv(S1A + "panel_genes.csv")
inter = [g for g in panel["gene_symbol"] if g in set(sym_l) and g in set(sym_a_present) and g in set(sym_b)]
inter = pd.Index(inter)
log(f"panel {len(panel)}  luca {len(sym_l)}  6149 present {len(sym_a_present)} (of {len(sym_a)} written)  6653 {len(sym_b)}")
log(f"INTERSECTION = {len(inter)} genes")

lost = panel[~panel["gene_symbol"].isin(inter)].copy()
lost["in_luca"] = lost["gene_symbol"].isin(set(sym_l))
lost["in_emtab6149"] = lost["gene_symbol"].isin(set(sym_a_present))
lost["in_emtab6149_written_as_structural_zero"] = lost["gene_symbol"].isin(set(sym_a)) & ~lost["gene_symbol"].isin(set(sym_a_present))
lost["in_emtab6653"] = lost["gene_symbol"].isin(set(sym_b))
lost["curated_not_hvg_only"] = lost["source_set"].str.replace("hvg_luca_2000", "", regex=False).str.strip(",").str.len() > 0
lost = lost.sort_values(["curated_not_hvg_only", "gene_symbol"], ascending=[False, True])
lost.to_csv(OUT_GENES, index=False)
log(f"lost {len(lost)} panel genes, {int(lost.curated_not_hvg_only.sum())} of them curated (non HVG-only)")

# ---------------------------------------------------------------- 2. remap tables
LEVEL2_MAP = {
    # LuCA vocabulary
    "AT1": "AT1", "AT2": "AT2", "B_cell": "B_cell", "CD4_T": "CD4_T", "CD8_T": "CD8_T",
    "DC": "DC", "NK": "NK", "SMC": "SMC", "Treg": "Treg", "ciliated": "ciliated",
    "club": "club", "endothelial": "endothelial",
    "endothelial_lymphatic": "endothelial_lymphatic", "epithelial_other": "epithelial_other",
    "fibroblast": "fibroblast", "macrophage": "macrophage", "malignant": "malignant",
    "mast": "mast", "mesothelial": "mesothelial", "monocyte": "monocyte",
    "myeloid_other": "myeloid_other", "neutrophil": "neutrophil", "pDC": "pDC",
    "pericyte": "pericyte", "plasma_cell": "plasma_cell",
    # E-MTAB-6149 vocabulary that differs
    "B": "B_cell", "plasma": "plasma_cell", "lymphatic_EC": "endothelial_lymphatic",
    "epithelial": "epithelial_other", "basal": "basal", "erythroid": "erythroid",
    # unlabelled arm
    "": "",
}
COMPARTMENTS = ["epithelial_malignant", "T_NK", "myeloid", "B_plasma", "fibroblast",
                "mural", "endothelial", "other"]

# ---------------------------------------------------------------- 3. per-arm load
def load_arm(fname, key):
    """Read one arm, restrict it to the intersection genes and return its CSR matrix and obs."""
    log(f"reading {fname}")
    A = ad.read_h5ad(S1A + fname)
    log(f"  {A.shape}  X nnz {A.X.nnz:,}  dtype {A.X.dtype}")
    if key == "LuCA_census":
        A.var_names = pd.Index(A.var["panel_symbol"].astype(str))
    A = A[:, inter].copy()
    X = A.X.tocsr()
    X.data = X.data.astype(np.float32, copy=False)
    X.eliminate_zeros()
    obs = A.obs.copy()
    del A
    return X, obs

X_l, obs_l = load_arm("arm_luca.h5ad", "LuCA_census")
X_a, obs_a = load_arm("arm_emtab6149.h5ad", "EMTAB6149")
X_b, obs_b = load_arm("arm_emtab6653.h5ad", "EMTAB6653")

# ---------------------------------------------------------------- 4. harmonize obs
def build_obs(obs, arm):
    """Harmonize the obs of one arm to the schema columns and the canonical block, donor and level-2 vocabularies."""
    o = pd.DataFrame(index=obs.index.astype(str))
    for c in SCHEMA_COLS:
        o[c] = obs[c].values
    o["donor_id_native"] = obs["donor_id"].astype(str).values
    o["block_id_native"] = obs["block_id"].astype(str).values
    o["level2_celltype_native"] = obs["level2_celltype"].astype(str).values
    if arm == "LuCA_census":
        o["patient_id"] = obs["patient_id"].astype(str).values
        o["macrophage_strict"] = pd.array(obs["macrophage_strict"].values, dtype="boolean")
        o["technical_replicate_group"] = ""
        o["tnm"] = ""
    elif arm == "EMTAB6149":
        o["patient_id"] = ("Lambrechts_Thienpont_2018_P" + obs["patient"].astype(str)).values
        o["macrophage_strict"] = pd.array([pd.NA] * len(o), dtype="boolean")
        o["technical_replicate_group"] = obs["technical_replicate_group"].astype(str).values
        o["tnm"] = obs["tnm"].astype(str).values
        # donor_id already Lambrechts_Thienpont_2018_P1..P5
    else:
        o["patient_id"] = ("Lambrechts_Thienpont_2018_P" + obs["patient_native"].astype(str)).values
        o["macrophage_strict"] = pd.array([pd.NA] * len(o), dtype="boolean")
        o["technical_replicate_group"] = ""
        o["tnm"] = obs["tnm"].astype(str).values
        # HARMONISE donor_id: the 6653 arm wrote Lambrechts_Thienpont_2018_6/7/8, the
        # 6149 arm wrote ..._P1.._P5.  Same study, same patient series (paper patients
        # 1 to 8).  Relabel to _P6/_P7/_P8 so the two accessions read as one series.
        # This is a rename only: the donor partition is unchanged.
        o["donor_id"] = ("Lambrechts_Thienpont_2018_P" + obs["patient_native"].astype(str)).values
    return o

o_l = build_obs(obs_l, "LuCA_census")
o_a = build_obs(obs_a, "EMTAB6149")
o_b = build_obs(obs_b, "EMTAB6653")
obs = pd.concat([o_l, o_a, o_b], axis=0)
del o_l, o_a, o_b, obs_l, obs_a, obs_b
log(f"obs concatenated {obs.shape}")

# level2 harmonisation, explicit and total
seen = set(obs["level2_celltype"].astype(str).unique())
missing = seen - set(LEVEL2_MAP)
if missing:
    die(f"level2_celltype values with no explicit map entry: {sorted(missing)}")
obs["level2_celltype"] = obs["level2_celltype"].astype(str).map(LEVEL2_MAP)

# level1 sanity: only the 8 allowed compartments, or empty for the unlabelled arm
bad = set(obs["level1_compartment"].astype(str).unique()) - set(COMPARTMENTS) - {""}
if bad:
    die(f"level1_compartment outside the allowed vocabulary: {sorted(bad)}")

# block_id harmonisation.  Canonical form is <study>_<platform>, per the obs schema.
# Where one study spans more than one source_key the accession is appended, because
# E-MTAB-6149 and E-MTAB-6653 are both Lambrechts 10x 3' v2 and must NOT be one block
# (median depth differs 2.2-fold, NQO1 detection moves 5 pct -> 31 pct between them).
study_srcs = obs.groupby("study", observed=True)["source_key"].nunique()
split_studies = set(study_srcs[study_srcs > 1].index)
log(f"studies spanning more than one source_key: {sorted(split_studies)}")
bid = obs["study"].astype(str) + "_" + obs["platform"].astype(str)
need = obs["study"].astype(str).isin(split_studies).values
bid = bid.where(~need, bid + "_" + obs["source_key"].astype(str))
obs["block_id"] = bid.values

# no block may span two source_keys after harmonisation
chk = obs.groupby("block_id", observed=True)["source_key"].nunique()
if (chk > 1).any():
    die(f"block_id spans multiple source_keys: {chk[chk>1].to_dict()}")
# no donor may span two histologies
chk = obs.groupby("donor_id", observed=True)["histology"].nunique()
if (chk > 1).any():
    die(f"donor_id spans multiple histologies: {chk[chk>1].to_dict()}")

for c in ["source_key", "study", "donor_id", "library_id", "histology", "histology_source",
          "histology_confidence", "tissue_class", "sampling_site", "platform", "modality",
          "enrichment", "treatment", "block_id", "level1_compartment", "level2_celltype",
          "celltype_source", "patient_id", "block_id_native", "level2_celltype_native",
          "technical_replicate_group", "tnm", "sex", "smoking"]:
    obs[c] = obs[c].astype(str).astype("category")
for c in ["level3_celltype_native", "donor_id_native", "stage"]:
    obs[c] = obs[c].astype(str)
for c in ["composition_safe", "counts_are_estimates"]:
    obs[c] = obs[c].astype(bool)

# ---------------------------------------------------------------- 5. matrix
log("vstacking")
X = sp.vstack([X_l, X_a, X_b], format="csr")
del X_l, X_a, X_b
log(f"merged X {X.shape}  nnz {X.nnz:,}  dtype {X.dtype}")
obs["n_counts_intersect"] = np.asarray(X.sum(axis=1)).ravel()
obs["n_genes_intersect"] = np.diff(X.indptr).astype(np.int32)

# ---------------------------------------------------------------- 6. assertions
log("ASSERTIONS")
fail = []

# (2) cell_id uniqueness
ids = obs.index
dup = ids.duplicated(keep=False)
if dup.any():
    d = pd.Series(ids[dup]).value_counts()
    die(f"cell_id collisions: {int(dup.sum())} rows, {len(d)} ids, e.g. {d.head(10).to_dict()}")
log(f"  cell_id unique: {len(ids):,} / {len(ids):,} PASS")
form_ok = pd.Series(ids).str.count(r"\|").eq(2)
log(f"  cell_id has exactly two pipes on {int(form_ok.sum()):,} of {len(ids):,} rows "
    f"({'PASS' if form_ok.all() else 'FAIL'})")
pref = pd.Series(ids).str.split("|").str[0].value_counts().to_dict()
log(f"  cell_id source prefixes: {pref}")

# (3) required columns present and non-empty
viol = []
for c in REQUIRED:
    if c not in obs.columns:
        viol.append((c, "ABSENT", len(obs))); continue
    s = obs[c]
    n_na = int(s.isna().sum())
    n_empty = int((s.astype(str).str.strip() == "").sum()) if s.dtype.name in ("object", "category") else 0
    if n_na or n_empty:
        viol.append((c, "empty/NA", n_na + n_empty))
if viol:
    for c, k, n in viol:
        by = obs.loc[obs[c].astype(str).str.strip() == "", "source_key"].value_counts().to_dict() if c in obs else {}
        log(f"  VIOLATION required column {c}: {n:,} cells ({k}); by source_key {by}")
else:
    log("  all required columns present and non-empty on every cell PASS")

# spec invariants
cs = obs["composition_safe"].values
wt = (obs["enrichment"].astype(str) == "whole_tissue").values
log(f"  composition_safe == (enrichment == whole_tissue): {'PASS' if (cs == wt).all() else 'FAIL'} "
    f"({int(cs.sum()):,} safe cells)")
log(f"  donor_id never empty: {'PASS' if (obs.donor_id.astype(str).str.strip() != '').all() else 'FAIL'}")
log(f"  platform/modality never empty: "
    f"{'PASS' if (obs.platform.astype(str).str.strip() != '').all() and (obs.modality.astype(str).str.strip() != '').all() else 'FAIL'}")
# substring trap re-check on the merged object
mal = obs["level3_celltype_native"].isin(["Cancer cells", "tumor"])
log(f"  malignant native labels all in epithelial_malignant: "
    f"{'PASS' if (obs.loc[mal, 'level1_compartment'] == 'epithelial_malignant').all() else 'FAIL'} ({int(mal.sum()):,} cells)")
mast = obs["level2_celltype"].astype(str) == "mast"
log(f"  mast cells all in other: "
    f"{'PASS' if (obs.loc[mast, 'level1_compartment'] == 'other').all() else 'FAIL'} ({int(mast.sum()):,} cells)")
log(f"  n_counts > 0 on every cell: {'PASS' if (obs.n_counts.values > 0).all() else 'FAIL'}")
log(f"  cells with zero counts in the 1,906-gene intersection: "
    f"{int((obs.n_counts_intersect.values == 0).sum()):,}")

# ---------------------------------------------------------------- 7. write
var = panel.set_index("gene_symbol").loc[inter].copy()
var.index.name = "gene_symbol"
var["ensembl_id_luca"] = var_l.set_index(pd.Index(var_l["panel_symbol"].astype(str))).loc[inter, "gene_id"].values
var["ensembl_id_emtab6149"] = var_a.loc[inter, "ensembl_id_in_source"].astype(str).values
var["gene_id_resolution_emtab6149"] = var_a.loc[inter, "gene_id_resolution"].astype(str).values
var["ensembl_id_emtab6653"] = var_b.loc[inter, "ensembl_id"].astype(str).values
var["ensembl_id_agrees_all_arms"] = (
    (var["ensembl_id_luca"] == var["ensembl_id_emtab6149"]) &
    (var["ensembl_id_luca"] == var["ensembl_id_emtab6653"])
)
log(f"  Ensembl id agrees across all three arms on {int(var.ensembl_id_agrees_all_arms.sum())} of {len(var)} genes")

A = ad.AnnData(X=X, obs=obs, var=var)
A.uns["stage"] = "1a"
A.uns["built_utc"] = pd.Timestamp.utcnow().isoformat()
A.uns["script"] = "atlas/steps/stage1a_merge_panel_object.py"
A.uns["arms"] = ["arm_luca.h5ad", "arm_emtab6149.h5ad", "arm_emtab6653.h5ad"]
A.uns["X_is"] = ("raw counts for LuCA_census and EMTAB6149; Alevin EM count estimates for "
                 "EMTAB6653 and for the two LuCA Smart-seq2 studies. float32. "
                 "Read obs.counts_are_estimates before any negative-binomial model.")
A.uns["gene_join"] = "gene_symbol (Ensembl ids drift between arms; see var.ensembl_id_* and gene_id_resolution_emtab6149)"
A.uns["panel_genes_lost"] = [str(g) for g in lost["gene_symbol"]]
A.uns["panel_genes_lost_curated"] = [str(g) for g in lost.loc[lost.curated_not_hvg_only, "gene_symbol"]]
A.uns["block_id_rule"] = ("<study>_<platform>, with _<source_key> appended for every block of a study "
                          "that spans more than one accession. Arm-native strings kept in obs.block_id_native.")
A.uns["donor_id_note"] = ("EMTAB6653 donors were relabelled from Lambrechts_Thienpont_2018_6/7/8 to "
                          "_P6/_P7/_P8 to match the EMTAB6149 convention. Rename only, partition unchanged. "
                          "obs.donor_id_native keeps the arm value.")
A.uns["unlabelled_cells"] = ("EMTAB6653 (32,341 cells) carries no cell type labels: the deposit has none. "
                             "level1_compartment, level2_celltype and level3_celltype_native are empty there. "
                             "Every compartment query must filter level1_compartment != ''.")
A.uns["level2_map"] = json.dumps(LEVEL2_MAP)

log("writing " + OUT_H5)
A.write_h5ad(OUT_H5)
log("written")

# ---------------------------------------------------------------- 8. summary tables
rows = []
o = obs
lab = o[o["level1_compartment"].astype(str) != ""]

rows.append(dict(table="overall", n_cells=len(o), n_donors=o.donor_id.nunique(),
                 n_patients=o.patient_id.nunique(), n_libraries=o.library_id.nunique(),
                 n_studies=o.study.nunique(), n_blocks=o.block_id.nunique(),
                 n_cells_composition_safe=int(o.composition_safe.sum()),
                 n_cells_labelled=len(lab), n_genes=len(inter)))

for h, g in o.groupby("histology", observed=True):
    gl = g[g["level1_compartment"].astype(str) != ""]
    rows.append(dict(table="histology", histology=h, n_cells=len(g), n_donors=g.donor_id.nunique(),
                     n_patients=g.patient_id.nunique(), n_libraries=g.library_id.nunique(),
                     n_blocks=g.block_id.nunique(),
                     n_cells_composition_safe=int(g.composition_safe.sum()),
                     n_cells_labelled=len(gl)))

for (h, c), g in lab.groupby(["histology", "level1_compartment"], observed=True):
    rows.append(dict(table="histology_x_compartment", histology=h, level1_compartment=c,
                     n_cells=len(g), n_donors=g.donor_id.nunique(),
                     n_cells_composition_safe=int(g.composition_safe.sum()),
                     n_donors_composition_safe=g.loc[g.composition_safe, "donor_id"].nunique()))

blk = o.groupby("block_id", observed=True)
both = set()
for b, g in blk:
    hs = set(g.histology.astype(str))
    if {"LUAD", "LUSC"} <= hs:
        both.add(b)
for b, g in blk:
    rows.append(dict(table="block", block_id=b, study=g.study.iloc[0], platform=g.platform.iloc[0],
                     source_key=g.source_key.iloc[0], n_cells=len(g), n_donors=g.donor_id.nunique(),
                     n_donors_LUAD=g.loc[g.histology == "LUAD", "donor_id"].nunique(),
                     n_donors_LUSC=g.loc[g.histology == "LUSC", "donor_id"].nunique(),
                     n_cells_LUAD=int((g.histology == "LUAD").sum()),
                     n_cells_LUSC=int((g.histology == "LUSC").sum()),
                     n_cells_composition_safe=int(g.composition_safe.sum()),
                     n_cells_labelled=int((g.level1_compartment.astype(str) != "").sum()),
                     holds_both_histologies=b in both))

wp = lab[lab.block_id.astype(str).isin(both) & lab.histology.isin(["LUAD", "LUSC"])]
for (b, c, h), g in wp.groupby(["block_id", "level1_compartment", "histology"], observed=True):
    rows.append(dict(table="within_protocol_block_x_compartment", block_id=b, level1_compartment=c,
                     histology=h, n_cells=len(g), n_donors=g.donor_id.nunique(),
                     n_cells_composition_safe=int(g.composition_safe.sum()),
                     n_donors_composition_safe=g.loc[g.composition_safe, "donor_id"].nunique(),
                     holds_both_histologies=True))

cols = ["table", "block_id", "study", "platform", "source_key", "histology", "level1_compartment",
        "n_donors", "n_donors_LUAD", "n_donors_LUSC", "n_patients", "n_libraries", "n_studies",
        "n_blocks", "n_genes", "n_cells", "n_cells_LUAD", "n_cells_LUSC",
        "n_cells_composition_safe", "n_donors_composition_safe", "n_cells_labelled",
        "holds_both_histologies"]
S = pd.DataFrame(rows)
for c in cols:
    if c not in S.columns:
        S[c] = np.nan
S = S[cols]
S.to_csv(OUT_SUM, index=False)
log(f"summary written {OUT_SUM} ({len(S)} rows)")

# ---------------------------------------------------------------- 9. printed tables
pd.set_option("display.width", 250)
print("\n=== CELLS BY HISTOLOGY x COMPARTMENT (labelled cells only) ===")
print(pd.crosstab(lab.level1_compartment, lab.histology, margins=True).to_string())
print("\n=== DONORS BY HISTOLOGY x COMPARTMENT ===")
print(lab.groupby(["level1_compartment", "histology"], observed=True).donor_id.nunique()
      .unstack(fill_value=0).to_string())
print("\n=== COMPOSITION-SAFE CELLS BY HISTOLOGY x COMPARTMENT ===")
sfe = lab[lab.composition_safe]
print(pd.crosstab(sfe.level1_compartment, sfe.histology, margins=True).to_string())
print("\n=== BLOCK TABLE ===")
bt = S[S.table == "block"][["block_id", "source_key", "platform", "n_cells", "n_donors",
                            "n_donors_LUAD", "n_donors_LUSC", "n_cells_LUAD", "n_cells_LUSC",
                            "n_cells_composition_safe", "n_cells_labelled",
                            "holds_both_histologies"]].sort_values("n_cells", ascending=False)
print(bt.to_string(index=False))
print("\n=== WITHIN-PROTOCOL BLOCKS (both LUAD and LUSC): donors per compartment ===")
w = S[S.table == "within_protocol_block_x_compartment"]
print(w.pivot_table(index=["block_id", "level1_compartment"], columns="histology",
                    values="n_donors", fill_value=0).to_string())
print("\n=== WITHIN-PROTOCOL BLOCKS: composition-safe donors per compartment ===")
print(w.pivot_table(index=["block_id", "level1_compartment"], columns="histology",
                    values="n_donors_composition_safe", fill_value=0).to_string())

log("done")
