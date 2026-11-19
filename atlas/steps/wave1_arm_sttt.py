#!/usr/bin/env python3
"""
Atlas stage 03, step arm_sttt: STTT / Zhang 2022, a large squamous single-cell source.

Paper items: part of the atlas behind Figure 2, Figure 3, Supplementary Figures 2 to 4 and Supplementary Data 1. No random numbers are drawn.

Zhang L et al., "Integrated single-cell RNA sequencing analysis reveals distinct cellular and
transcriptional modules associated with survival in lung cancer", Signal Transduction and Targeted
Therapy 2022, PMID 35027529.

The paper carries no data availability section and no accession. The data is served by the authors'
portal at lungcancer.chenlulab.com behind `api/download?name=` (see atlas/sources.tsv). A range-capable
GET works and resumes (00_download.sh uses curl -C -).

WHAT IT ADDS. 2,233 squamous TUMOR fibroblasts across 14 LUSC patients, 8 of them above 20 cells.
It is a single both-histology block, 14 LUSC against 22 LUAD on one platform, and it
carries matched normal tissue from the same donors.

INPUT is a dense CSV, 17,277 genes as rows by 293,432 cells as columns, streamed a line at a time so
the 10 GB uncompressed matrix never lands in memory. Only panel-gene rows are parsed.

CELL TYPE MAP is explicit over all 17 author labels with a completeness assertion. Two are ambiguous
from the label alone and are resolved by marker expression at the end of this script rather than by
assumption: `NE` and `Gran`.

RUN: python atlas/steps/wave1_arm_sttt.py   (normally through atlas/03_arms_wave1.py)
OUT -> output/RNA/meta_dataset/wave1/arm_sttt.h5ad
"""
import os, gzip, time, json
import numpy as np, pandas as pd, h5py
from scipy.sparse import coo_matrix
from anndata.io import write_elem

T0 = time.time()
import os as _os  # noqa: E402
import sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # atlas/
from common import ROOT, SCRNA, SCRATCH, DATA, RSCRIPT  # noqa: E402,F401  paths from config.json
SRC = f"{SCRNA}/STTT_Zhang2022"
OBJ = f"{ROOT}/output/RNA/meta_dataset/stage1a/panel_object.h5ad"
OUT = f"{ROOT}/output/RNA/meta_dataset/wave1"
os.makedirs(OUT, exist_ok=True)

CT = {
 "Fib":   ("fibroblast", "fibroblast"),
 "EC":    ("endothelial", "endothelial"),
 "Mo":    ("myeloid", "macrophage"),
 "DC":    ("myeloid", "DC"),
 "Gran":  ("myeloid", "neutrophil"),
 "Mast":  ("other", "mast"),
 "CD4":   ("T_NK", "CD4_T"),
 "CD8":   ("T_NK", "CD8_T"),
 "Tregs": ("T_NK", "Treg"),
 "NK":    ("T_NK", "NK"),
 "B":     ("B_plasma", "B_cell"),
 "ATII":  ("epithelial_malignant", "AT2"),
 "Basal": ("epithelial_malignant", "basal"),
 "Club":  ("epithelial_malignant", "club"),
 "Cilia": ("epithelial_malignant", "ciliated"),
 "Epi":   ("epithelial_malignant", "epithelial_other"),
 "NE":    ("epithelial_malignant", "epithelial_other"),   # neuroendocrine; verified below
}
HMAP = {"LUAD": "LUAD", "LUSC": "LUSC", "LUAD_Normal": "LUAD", "LUSC_Normal": "LUSC",
        "LULC": "large_cell", "LULC_Normal": "large_cell",
        "Pleiomorphic": "pleomorphic", "Pleiomorphic Normal": "pleomorphic"}
TCLASS = {"Tumor": "primary_tumour", "Normal": "adjacent_normal"}
MARKERS = {"NE": ["CHGA", "CALCA", "SYP", "ASCL1"], "Gran": ["FCGR3B", "CSF3R", "S100A8", "CXCR2"]}


def log(*a):
    """Print a message with the elapsed time."""
    print(f"[{time.time()-T0:7.1f}s]", *a, flush=True)


with h5py.File(OBJ, "r") as h:
    PANEL = [x.decode() if isinstance(x, bytes) else str(x) for x in h["var/gene_symbol"][:]]
PSET = set(PANEL)
meta = pd.read_csv(f"{SRC}/STTT_meta.csv", low_memory=False).set_index("Cells")
assert set(meta.CellName.unique()) <= set(CT), f"unmapped: {set(meta.CellName.unique())-set(CT)}"
log(f"metadata {len(meta):,} cells; cell-type map covers all {meta.CellName.nunique()} labels")

need = PSET | {g for gl in MARKERS.values() for g in gl}
rows, cols, vals, kept = [], [], [], []
with gzip.open(f"{SRC}/raw_data.csv.gz", "rt") as f:
    header = next(f).rstrip("\n").split(",")
    cells = [c.strip('"') for c in header[1:]]
    log(f"matrix: {len(cells):,} cells")
    for ln, line in enumerate(f):
        i = line.index(",")
        g = line[:i].strip('"')
        if g not in need:
            continue
        v = np.fromstring(line[i + 1:], sep=",", dtype=np.float32)
        nz = np.nonzero(v)[0]
        if len(nz):
            r = len(kept)
            rows.append(np.full(len(nz), r, np.int32)); cols.append(nz.astype(np.int32))
            vals.append(v[nz])
        kept.append(g)
        if len(kept) % 250 == 0:
            log(f"  kept {len(kept)} genes, scanned {ln+1} rows")
log(f"parsed {len(kept)} needed genes of {len(need)} sought")

gidx = {g: i for i, g in enumerate(kept)}
X = coo_matrix((np.concatenate(vals), (np.concatenate(cols), np.concatenate(rows))),
               shape=(len(cells), len(kept)), dtype=np.float32).tocsr()
del rows, cols, vals
log(f"X {X.shape} nnz {X.nnz:,}")

# resolve the two ambiguous labels from markers before committing the map
lab = meta.reindex(cells).CellName.astype(str).values
for k, ms in MARKERS.items():
    ms = [m for m in ms if m in gidx]
    if not ms: continue
    sub = X[:, [gidx[m] for m in ms]]
    tot = np.asarray(X.sum(axis=1)).ravel(); tot[tot == 0] = 1
    frac = np.asarray(sub.sum(axis=1)).ravel() / tot
    log(f"  {k} marker check ({ms}): in {k} cells {1e4*frac[lab==k].mean():.1f} vs "
        f"elsewhere {1e4*frac[lab!=k].mean():.1f} (x1e4 of counts)")

pcols = [gidx[g] for g in PANEL if g in gidx]
present = [g for g in PANEL if g in gidx]
log(f"panel genes present: {len(present)} / {len(PANEL)}; absent {sorted(set(PANEL)-set(present))[:12]}")
X = X[:, pcols]

m = meta.reindex(cells)
obs = pd.DataFrame(index=[f"STTT|{s}|{c.split('_',1)[1]}" for s, c in zip(m.SampleID.astype(str), cells)])
obs["source_key"] = "STTT"
obs["study"] = "Zhang_Chen_2022"
obs["donor_id"] = ["Zhang_Chen_2022_" + str(p) for p in m.PatientID]
obs["donor_id_native"] = m.PatientID.astype(str).values
obs["patient_id"] = obs["donor_id"].values
obs["library_id"] = m.SampleID.astype(str).values
obs["histology"] = m.Disease.astype(str).map(HMAP).values
obs["histology_source"] = "paper_supplement"
obs["histology_confidence"] = "high"
obs["tissue_class"] = m.Tissue.astype(str).map(TCLASS).values
obs["sampling_site"] = "NA"
obs["platform"] = "10x 3' v2"
obs["modality"] = "single_cell"
obs["enrichment"] = "whole_tissue"
obs["composition_safe"] = True
obs["treatment"] = "unknown"
obs["block_id"] = "Zhang_Chen_2022_10x 3' v2"
obs["level1_compartment"] = [CT[x][0] for x in m.CellName.astype(str)]
obs["level2_celltype"] = [CT[x][1] for x in m.CellName.astype(str)]
obs["level2_celltype_map_native"] = obs["level2_celltype"].values
obs["level3_celltype_native"] = m.CellName.astype(str).values
obs["celltype_source"] = "author"
obs["celltype_confidence"] = np.nan
obs["is_malignant"] = pd.array([pd.NA] * len(obs), dtype="boolean")
obs["malignancy_score"] = np.nan
obs["n_counts"] = np.asarray(X.sum(axis=1)).ravel().astype("int64")
obs["n_genes"] = np.diff(X.indptr).astype("int32")
obs["counts_are_estimates"] = False
obs["doublet_score"] = np.nan
obs["sex"] = m.Gender.astype(str).str.lower().values
obs["age"] = m.Age.values
obs["stage"] = m.Stage.astype(str).values
obs["smoking"] = "unknown"
obs["tnm"] = (m["T"].astype(str) + m["N"].astype(str) + m["M"].astype(str)).values
obs["technical_replicate_group"] = "batch" + m.Batch.astype(str).values
obs["macrophage_strict"] = pd.array([False] * len(obs), dtype="boolean")

drop = obs.n_counts.values == 0
if drop.any():
    log(f"dropping {int(drop.sum()):,} cells with zero panel counts")
    X = X[~drop]; obs = obs[~drop]
assert obs.index.is_unique and (obs.index.str.count(r"\|") == 2).all()
assert (obs.block_id == obs.study + "_" + obs.platform).all()
log(f"obs {obs.shape}; median panel counts {int(np.median(obs.n_counts)):,}, "
    f"median panel genes {int(np.median(obs.n_genes)):,}")
log("compartments:\n" + obs.level1_compartment.value_counts().to_string())
log("fibroblasts by histology x tissue:\n" + pd.crosstab(
    obs[obs.level1_compartment=="fibroblast"].histology,
    obs[obs.level1_compartment=="fibroblast"].tissue_class).to_string())

var = pd.DataFrame(index=pd.Index(present, name="gene_symbol")); var["present_in_source"] = True
uns = {"accession": "portal lungcancer.chenlulab.com (no repository accession published)",
       "study_full": "Zhang L et al., Signal Transduct Target Ther 2022, PMID 35027529",
       "X_is": "integer counts, panel genes only",
       "provenance_note": ("The paper has no data availability section. Files come from the authors' "
                           "portal endpoint api/download?name=, 484,663,258 bytes, gzip verified."),
       "source_gene_space": "17,277 gene symbols in the deposited matrix",
       "celltype_map": json.dumps({k: list(v) for k, v in CT.items()})}
path = f"{OUT}/arm_sttt.h5ad"
with h5py.File(path, "w") as fo:
    fo.attrs["encoding-type"] = "anndata"; fo.attrs["encoding-version"] = "0.1.0"
    write_elem(fo, "obs", obs); write_elem(fo, "var", var); write_elem(fo, "uns", uns)
    for k in ["layers", "obsm", "obsp", "varm", "varp"]: write_elem(fo, k, {})
    write_elem(fo, "X", X)
log(f"WROTE {path} ({os.path.getsize(path):,} bytes)")
