#!/usr/bin/env python3
"""
Atlas stage 03, step arm_gse246215: GSE246215, the Gao 2024 pan-cancer fibroblast atlas, NSCLC arm only.

Paper items: part of the atlas behind Figure 2, Figure 3, Supplementary Figures 2 to 4 and Supplementary Data 1. No random numbers are drawn.

WHAT IT IS. 13 donors, every one with at least 20 fibroblasts, all deposited as fibroblasts so no cell
typing is needed, all primary tumor so tissue_class is known; its `Dataset` column reads "this study"
for all 7,244 NSCLC cells, so none of it duplicates a study already in the panel object.

IT IS BUILT AS A SIDECAR, NOT MERGED INTO THE PANEL OBJECT. The matrix lacks 18 of the object's 1,906
panel genes, genuinely rather than by symbol renaming: it carries neither CCN1 nor its old symbol CYR61,
neither CCN2 nor CTGF, and so on. Merging would shrink the object's gene axis, and because every
expression proportion divides by the panel total, that would change every number. The
18 genes carry 1.31 per cent of all fibroblast panel counts, CCN1 and CCN2 alone carrying most of it.
Since every analysis is computed WITHIN a block and combined by block votes, a new block does
not need to share a gene axis with the others. It needs to carry the gene set being scored (the
full-gene build, stage 04, carries a presence mask for this).

Sixteen of the 18 are recovered anyway: the alias map below resolves the ones that are only renames.

RUN: python atlas/steps/wave1_arm_gse246215.py   (normally through atlas/03_arms_wave1.py)
FROM: GSE246215 (Fibroblast counts csv, metadata) and the gene axis of stage1a/panel_object.h5ad.
OUT -> output/RNA/meta_dataset/wave1/arm_gse246215.h5ad
"""
import os, time, json
import numpy as np
import pandas as pd
import h5py
from scipy.sparse import csr_matrix
from anndata.io import write_elem

T0 = time.time()
import os as _os  # noqa: E402
import sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # atlas/
from common import ROOT, SCRNA, SCRATCH, DATA, RSCRIPT  # noqa: E402,F401  paths from config.json
SRC = f"{SCRNA}/GSE246215"
OBJ = f"{ROOT}/output/RNA/meta_dataset/stage1a/panel_object.h5ad"
OUT = f"{ROOT}/output/RNA/meta_dataset/wave1"
os.makedirs(OUT, exist_ok=True)

# previous HGNC symbols, applied only where the alias is in this matrix AND is not itself a panel gene
ALIAS = {"CEP43": "FGFR1OP", "CIBAR2": "FAM92B", "CILK1": "ICK", "HOATZ": "C11orf88",
         "PABIR3": "FAM122C", "PEDS1": "TMEM189", "TAMALIN": "GRASP", "TLCD3A": "FAM57A",
         "SLC66A1": "PQLC2", "SLC49A4": "DIRC2", "POGLUT3": "KDELC2", "CRPPA": "ISPD",
         "ECRG4": "C2orf40", "BNIP5": "C6orf222", "FLACC1": "ALS2CR11", "H1-1": "HIST1H1A",
         "H1-2": "HIST1H1C", "H2AC20": "HIST2H2AA3", "H4-16": "HIST4H4", "H4C3": "HIST1H4C",
         "KRT86": "KRTHB6", "VARS1": "VARS", "CCN1": "CYR61", "CCN2": "CTGF",
         "PALM2AKAP2": "PALM2-AKAP2"}


def log(*a):
    """Print a message with the elapsed time."""
    print(f"[{time.time()-T0:7.1f}s]", *a, flush=True)


with h5py.File(OBJ, "r") as h:
    PANEL = [x.decode() if isinstance(x, bytes) else str(x) for x in h["var/gene_symbol"][:]]
PSET = set(PANEL)
log(f"object panel: {len(PANEL)} genes")

meta = pd.read_csv(f"{SRC}/GSE246215_Fibroblast_metadata.csv.gz")
meta = meta[meta.CancerType_short == "NSCLC"].copy()
log(f"GSE246215 NSCLC fibroblasts: {len(meta):,} cells, {meta.PatientID.nunique()} donors, "
    f"provenance {sorted(meta.Dataset.unique())}")
assert set(meta.Dataset.unique()) == {"this study"}, "unexpected provenance, would duplicate a study"
assert set(meta.Phenotype.unique()) == {"Tumor"}
keep_cells = list(meta.CellName)

# ---- read the counts matrix: genes as rows, cells as columns, keeping only what we need
header = pd.read_csv(f"{SRC}/GSE246215_Fibroblast_counts.csv.gz", nrows=0)
allcols = list(header.columns)
log(f"matrix columns: {len(allcols):,} (first is {allcols[0]})")
want = [allcols[0]] + [c for c in keep_cells if c in set(allcols)]
log(f"NSCLC columns found in matrix: {len(want)-1:,} of {len(keep_cells):,}")

blocks, kept_genes = [], []
n_read = 0
for chunk in pd.read_csv(f"{SRC}/GSE246215_Fibroblast_counts.csv.gz", usecols=want,
                         index_col=0, chunksize=4000, dtype={allcols[0]: str}):
    n_read += len(chunk)
    sym = pd.Index(chunk.index.astype(str))
    hit = sym.isin(PSET) | sym.isin(set(ALIAS.values()))
    if hit.any():
        sub = chunk.loc[hit]
        blocks.append(sub.astype(np.float32))
        kept_genes += list(sub.index.astype(str))
    log(f"  rows read {n_read:,}   kept so far {len(kept_genes)}")
M = pd.concat(blocks, axis=0)
del blocks
log(f"matrix subset: {M.shape[0]} genes x {M.shape[1]:,} cells")

# rename aliases back to the panel symbol, then order as the panel
rev = {v: k for k, v in ALIAS.items()}
M.index = [rev.get(g, g) for g in M.index.astype(str)]
M = M[~M.index.duplicated(keep="first")]
present = [g for g in PANEL if g in set(M.index)]
lost = [g for g in PANEL if g not in set(M.index)]
log(f"panel genes present after alias resolution: {len(present)} / {len(PANEL)}; LOST {len(lost)}: {lost}")
M = M.loc[present]

X = csr_matrix(M.T.values.astype(np.float32))          # cells x genes
cells = list(M.columns)
del M
log(f"X built {X.shape}, nnz {X.nnz:,}")

m = meta.set_index("CellName").loc[cells]
barcode = [c.split("_", 1)[1] if "_" in c else c for c in cells]
lib = [c.split("_", 1)[0] if "_" in c else c for c in cells]
obs = pd.DataFrame(index=[f"GSE246215|{l}|{b}" for l, b in zip(lib, barcode)])
obs["source_key"] = "GSE246215"
obs["study"] = "Gao_Zhang_2024"
obs["donor_id"] = ["Gao_Zhang_2024_" + str(p) for p in m.PatientID]
obs["donor_id_native"] = m.PatientID.astype(str).values
obs["patient_id"] = obs["donor_id"].values
obs["library_id"] = m.SampleID.astype(str).values
obs["histology"] = "NSCLC_unspecified"
obs["histology_source"] = "not_determinable"
obs["histology_confidence"] = "unknown"
obs["tissue_class"] = "primary_tumour"
obs["sampling_site"] = "NA"
obs["platform"] = "10x 5' v1"
obs["modality"] = "single_cell"
obs["enrichment"] = "stromal_sort"
obs["composition_safe"] = False
obs["treatment"] = "unknown"
obs["block_id"] = "Gao_Zhang_2024_10x 5' v1"
obs["level1_compartment"] = "fibroblast"
obs["level2_celltype"] = "fibroblast"
obs["level2_celltype_map_native"] = "fibroblast"
obs["level3_celltype_native"] = ["Gao_" + str(c) for c in m.Cluster]
obs["celltype_source"] = "author"
obs["celltype_confidence"] = np.nan
obs["is_malignant"] = pd.array([False] * len(obs), dtype="boolean")
obs["malignancy_score"] = np.nan
obs["n_counts"] = np.asarray(X.sum(axis=1)).ravel().astype("int64")
obs["n_genes"] = np.diff(X.indptr).astype("int32")
obs["counts_are_estimates"] = False
obs["doublet_score"] = np.nan
obs["sex"] = "unknown"; obs["age"] = np.nan; obs["stage"] = ""; obs["smoking"] = "unknown"
obs["tnm"] = ""; obs["technical_replicate_group"] = ""
obs["macrophage_strict"] = pd.array([False] * len(obs), dtype="boolean")

assert obs.index.is_unique, "cell_id collision"
assert (obs.index.str.count(r"\|") == 2).all(), "cell_id must have exactly two pipes"
assert (obs.block_id == obs.study + "_" + obs.platform).all()
assert (obs.n_counts.values > 0).all(), "cells with zero panel counts"
log(f"obs built {obs.shape}; median panel counts {int(np.median(obs.n_counts)):,}, "
    f"median panel genes {int(np.median(obs.n_genes)):,}")

var = pd.DataFrame(index=pd.Index(present, name="gene_symbol"))
var["present_in_source"] = True
var["resolved_via_alias"] = [g in ALIAS and ALIAS[g] not in PSET for g in present]
uns = {"accession": "GSE246215", "study_full": "Gao et al. 2024 pan-cancer fibroblast atlas, NSCLC arm",
       "X_is": "integer UMI counts stored as float32, panel genes only",
       "sidecar_note": (f"NOT merged into panel_object.h5ad. The matrix genuinely lacks {len(lost)} "
                        "panel genes, so merging would shrink the object gene axis and change every "
                        "expression proportion in the project. Used as an additional BLOCK."),
       "panel_genes_absent": json.dumps(lost),
       "enrichment_note": ("The deposited object contains fibroblasts only, computationally selected "
                           "by the authors. Labelled stromal_sort so composition_safe is False; no "
                           "composition query may use it."),
       "histology_note": "Gao's tables carry no LUAD/LUSC split, so histology is NSCLC_unspecified."}

path = f"{OUT}/arm_gse246215.h5ad"
with h5py.File(path, "w") as fo:
    fo.attrs["encoding-type"] = "anndata"; fo.attrs["encoding-version"] = "0.1.0"
    write_elem(fo, "obs", obs); write_elem(fo, "var", var); write_elem(fo, "uns", uns)
    for k in ["layers", "obsm", "obsp", "varm", "varp"]:
        write_elem(fo, k, {})
    write_elem(fo, "X", X)
log(f"WROTE {path} ({os.path.getsize(path):,} bytes)")
