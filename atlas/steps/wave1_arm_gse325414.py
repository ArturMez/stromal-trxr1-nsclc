#!/usr/bin/env python3
"""
Atlas stage 03, step arm_gse325414: GSE325414, BD Rhapsody with CITE, 156,467 cells, 25 donors.

Paper items: part of the atlas behind Figure 2, Figure 3, Supplementary Figures 2 to 4 and Supplementary Data 1. No random numbers are drawn.

WHY IT IS WORTH BUILDING.
  1. It is merge-safe. Of the object's 1,906 panel genes it is missing only three (EEF1AKNMT, FLACC1,
     NICN1), together 0.005 per cent of fibroblast panel counts and none of them in a curated set, so
     unlike GSE246215 it could be merged into the object without moving any number.
  2. It carries author annotation on every compartment, 43 leaf labels, so it adds myeloid and T cells
     as well as fibroblasts.
  3. It carries resection tumor and resection uninvolved lung from the same donors, so tumor and
     adjacent normal tissue can be compared within donor.

HISTOLOGY. The deposited cell metadata carries no LUAD/LUSC label; it comes from the GEO series
matrix, per donor, through wave1/gse325414_donor_meta.csv. The donor table and the two panel inputs
below (the row map and the panel triplets in scratch/recon) are written by wave1_gse325414_prepare.py.

CELL TYPE MAP is an explicit dict over all 43 author leaf labels with a completeness assertion. No
substring matching: the label 'Tcell_CD4-_CD8-' contains 'Tcell' but is not CD4 or CD8,
and 'MastCells' must land in compartment other, not myeloid.

TISSUE CLASS is inferred from the author sample-type code and marked as inferred: RESRT resection
tumor, RESRU resection uninvolved lung, INXBX index biopsy. The native code is kept verbatim.

RUN: python atlas/steps/wave1_arm_gse325414.py   (normally through atlas/03_arms_wave1.py)
OUT -> output/RNA/meta_dataset/wave1/arm_gse325414.h5ad
"""
import os, time, json, gzip
import numpy as np, pandas as pd, h5py
from scipy.sparse import coo_matrix
from anndata.io import write_elem

T0 = time.time()
import os as _os  # noqa: E402
import sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # atlas/
from common import ROOT, SCRNA, SCRATCH, DATA, RSCRIPT  # noqa: E402,F401  paths from config.json
SRC = f"{SCRNA}/GSE325414"
TRI = f"{SCRATCH}/recon/g325414_panel.tsv.gz"
ROWMAP = f"{SCRATCH}/recon/g325414_rowmap.json"
OBJ = f"{ROOT}/output/RNA/meta_dataset/stage1a/panel_object.h5ad"
OUT = f"{ROOT}/output/RNA/meta_dataset/wave1"
os.makedirs(OUT, exist_ok=True)


def log(*a):
    """Print a message with the elapsed time."""
    print(f"[{time.time()-T0:7.1f}s]", *a, flush=True)


CT = {
 "EpithelialCells_TumorCells": ("epithelial_malignant", "epithelial_other"),
 "Fibroblasts": ("fibroblast", "fibroblast"),
 "SmoothMuscleCells": ("mural", "SMC"),
 "EndothelialCells": ("endothelial", "endothelial"),
 "Neutrophils": ("myeloid", "neutrophil"),
 "Macrophage_monocyte": ("myeloid", "macrophage"),
 "DendriticCell_cDC": ("myeloid", "DC"),
 "DendriticCell_pDC": ("myeloid", "pDC"),
 "MastCells": ("other", "mast"),
 "NKcells_CD56dim": ("T_NK", "NK"), "NKcells_CD56bright": ("T_NK", "NK"),
 "ILCs": ("T_NK", "T_other"),
 "Tcell_CD4_Treg": ("T_NK", "Treg"),
 "Tcell_CD4+_CD8+": ("T_NK", "T_other"), "Tcell_CD4-_CD8-": ("T_NK", "T_other"),
 "Unidentified_NonHematopoietic": ("other", ""),
 "Bcells_Naive": ("B_plasma", "B_cell"),
}
for s in ["CM", "Naive", "PM", "TE", "TFH", "TM"]:
    CT[f"Tcell_CD4_{s}"] = ("T_NK", "CD4_T")
for s in ["CM", "EE", "EM", "Naive", "TE", "TM"]:
    CT[f"Tcell_CD8_{s}"] = ("T_NK", "CD8_T")
for s in ["DN", "DN_IgA", "DN_IgG", "DN_IgM", "IgA", "IgD", "IgE", "IgG", "IgMD"]:
    CT[f"Bcell_MBC_{s}"] = ("B_plasma", "B_cell")
for s in ["IgA", "IgD", "IgE", "IgG", "IgM"]:
    CT[f"Bcell_PC_{s}"] = ("B_plasma", "plasma_cell")
TCLASS = {"RESRT": "primary_tumour", "RESRU": "adjacent_normal", "INXBX": "primary_tumour"}

with h5py.File(OBJ, "r") as h:
    PANEL = [x.decode() if isinstance(x, bytes) else str(x) for x in h["var/gene_symbol"][:]]
rowmap = {int(k): v for k, v in json.load(open(ROWMAP)).items()}
present = [g for g in PANEL if g in set(rowmap.values())]
gpos = {g: i for i, g in enumerate(present)}
row2col = {r: gpos[s] for r, s in rowmap.items()}
log(f"panel genes in this source: {len(present)} / {len(PANEL)}")

bcs = [l.strip() for l in gzip.open(f"{SRC}/GSE325414_barcodes.tsv.gz", "rt")]
meta = pd.read_csv(f"{SRC}/GSE325414_metadata_individual_cells.txt.gz", sep="\t", low_memory=False)
meta = meta.set_index("cell")
log(f"barcodes {len(bcs):,}   metadata rows {len(meta):,}   overlap {len(set(bcs)&set(meta.index)):,}")
seen = set(meta["sub.pop.level4"].astype(str).unique())
missing = seen - set(CT)
assert not missing, f"author labels with no explicit map entry: {sorted(missing)}"
log(f"cell-type map covers all {len(seen)} author leaf labels")

log("reading triplets ...")
T = pd.read_csv(TRI, sep="\t", header=None, names=["r", "c", "v"],
                dtype={"r": np.int32, "c": np.int32, "v": np.float32})
log(f"  {len(T):,} entries")
gcol = np.array([row2col[r] for r in T.r.values], dtype=np.int32)
X = coo_matrix((T.v.values, (T.c.values - 1, gcol)),
               shape=(len(bcs), len(present))).tocsr()
del T, gcol
log(f"X built {X.shape}, nnz {X.nnz:,}")

keep = np.array([b in meta.index for b in bcs])
X = X[keep]
bck = [b for b in bcs if b in meta.index]
m = meta.loc[bck]
log(f"cells with metadata: {X.shape[0]:,}")

lin = m["sub.pop.level4"].astype(str).map(lambda x: CT[x][0])
l2 = m["sub.pop.level4"].astype(str).map(lambda x: CT[x][1])
obs = pd.DataFrame(index=[f"GSE325414|{s}|{b}" for s, b in zip(m["sample"].astype(str), bck)])
obs["source_key"] = "GSE325414"
obs["study"] = "Jimenez_Bruno_2026"   # Jimenez et al., Clin Cancer Res 2026, PMID 42307634
obs["donor_id"] = ["Jimenez_Bruno_2026_" + str(d) for d in m.donor]
obs["donor_id_native"] = m.donor.astype(str).values
obs["patient_id"] = obs["donor_id"].values
obs["library_id"] = m["sample"].astype(str).values
DMETA = pd.read_csv(f"{ROOT}/output/RNA/meta_dataset/wave1/gse325414_donor_meta.csv").set_index("donor")
HMAP = {"Adenocarcinoma": "LUAD", "Mucinous adenocarcinoma": "LUAD",
        "Squamous cell carcinoma": "LUSC"}
assert set(DMETA["histologic subtype"]) <= set(HMAP), "unmapped histologic subtype"
hh = DMETA["histologic subtype"].map(HMAP)
obs["histology"] = hh.reindex(m.donor.astype(str).values).values
obs["histology_source"] = "GEO"
obs["histology_confidence"] = ["inferred" if x == "Mucinous adenocarcinoma" else "high"
                              for x in DMETA["histologic subtype"].reindex(m.donor.astype(str).values)]
obs["tissue_class"] = m["sample.type.3"].astype(str).map(TCLASS).values
obs["sampling_site"] = "NA"
obs["platform"] = "BD Rhapsody"
obs["modality"] = "single_cell"
obs["enrichment"] = "whole_tissue"
obs["composition_safe"] = True
obs["treatment"] = "treated"
obs["block_id"] = "Jimenez_Bruno_2026_BD Rhapsody"
obs["level1_compartment"] = lin.values
obs["level2_celltype"] = l2.values
obs["level2_celltype_map_native"] = l2.values
obs["level3_celltype_native"] = m["sub.pop.level4"].astype(str).values
obs["celltype_source"] = "author"
obs["celltype_confidence"] = np.nan
obs["is_malignant"] = pd.array([pd.NA] * len(obs), dtype="boolean")
obs["malignancy_score"] = np.nan
obs["n_counts"] = np.asarray(X.sum(axis=1)).ravel().astype("int64")
obs["n_genes"] = np.diff(X.indptr).astype("int32")
obs["counts_are_estimates"] = False
obs["doublet_score"] = np.nan
obs["sex"] = DMETA["Sex"].reindex(m.donor.astype(str).values).str.lower().values
obs["age"] = DMETA["age"].reindex(m.donor.astype(str).values).values; obs["stage"] = DMETA["nsclc clinical_stage"].reindex(m.donor.astype(str).values).astype(str).values; obs["smoking"] = "unknown"
obs["tnm"] = ""; obs["technical_replicate_group"] = ""
obs["macrophage_strict"] = pd.array(
    [x == "Macrophage_monocyte" for x in m["sub.pop.level4"].astype(str)], dtype="boolean")
obs["macrophage_strict"] = pd.array([False] * len(obs), dtype="boolean")   # author merged mono+macro

drop = obs.n_counts.values == 0
if drop.any():
    log(f"dropping {int(drop.sum()):,} cells with zero panel counts")
    X = X[~drop]; obs = obs[~drop]
assert obs.index.is_unique and (obs.index.str.count(r"\|") == 2).all()
assert (obs.block_id == obs.study + "_" + obs.platform).all()
log(f"obs {obs.shape}; median panel counts {int(np.median(obs.n_counts)):,}, "
    f"median panel genes {int(np.median(obs.n_genes)):,}")
log("compartments:\n" + obs.level1_compartment.value_counts().to_string())
log("fibroblasts by tissue_class:\n" +
    obs[obs.level1_compartment == "fibroblast"].tissue_class.value_counts().to_string())

var = pd.DataFrame(index=pd.Index(present, name="gene_symbol"))
var["present_in_source"] = True
uns = {"accession": "GSE325414", "study_full": "BD Rhapsody CITE NSCLC, 25 donors",
       "X_is": "integer counts, panel genes only",
       "merge_safe": "Missing only EEF1AKNMT, FLACC1, NICN1 of the 1,906 panel genes.",
       "histology_note": ("From the GEO series matrix, downloaded 31 Aug 2026: 15 adenocarcinoma, 2 mucinous "
                          "adenocarcinoma mapped to LUAD and flagged inferred, 8 squamous. 25 donors."),
       "treatment_note": ("Study group PEF Treatment, treatment 'Ablation after diagnostic biopsy'. The "
                          "RESRT and RESRU samples are POST-ablation resections; INXBX is the pre-treatment "
                          "index biopsy. Any analysis sensitive to treatment should restrict to INXBX."),
       "tissue_class_note": "Inferred from sample.type.3: RESRT primary_tumour, RESRU adjacent_normal, "
                            "INXBX primary_tumour. Native code kept in library_id and level3.",
       "macrophage_note": "Authors merged macrophages and monocytes into one label, so level2 is "
                          "macrophage with macrophage_strict False on every cell.",
       "celltype_map": json.dumps({k: list(v) for k, v in CT.items()})}

path = f"{OUT}/arm_gse325414.h5ad"
with h5py.File(path, "w") as fo:
    fo.attrs["encoding-type"] = "anndata"; fo.attrs["encoding-version"] = "0.1.0"
    write_elem(fo, "obs", obs); write_elem(fo, "var", var); write_elem(fo, "uns", uns)
    for k in ["layers", "obsm", "obsp", "varm", "varp"]:
        write_elem(fo, k, {})
    write_elem(fo, "X", X)
log(f"WROTE {path} ({os.path.getsize(path):,} bytes)")
pos = np.asarray(X[:, gpos["TXNRD1"]].todense()).ravel() > 0
for L in ["fibroblast", "myeloid", "T_NK", "epithelial_malignant"]:
    s = obs.level1_compartment.values == L
    log(f"  TXNRD1+ in {L}: {100*pos[s].mean():.1f}% of {int(s.sum()):,}")
