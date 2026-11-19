#!/usr/bin/env python
"""Atlas stage 01, step verify_panel_object: independent verification of panel_object.h5ad.

Paper items: none directly; it checks the three-arm panel object (1,146,594 cells x 1,906 genes) of the
atlas behind Figure 2, Figure 3, Supplementary Figures 2 to 4 and Supplementary Data 1. Every number
is re-derived from the written file. Prints per-arm cell counts, donors and blocks by histology, and
compartment tables; no random numbers are drawn.
Run: python atlas/steps/stage1a_verify_panel_object.py   (normally through atlas/01_arms_luca_lambrechts.py)
"""
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad, h5py, time
from anndata.io import read_elem

T0 = time.time()
import os as _os  # noqa: E402
import sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # atlas/
from common import MD  # noqa: E402  paths from config.json, see atlas/common.py
S1A = MD + "/stage1a/"
def log(*a):
    """Print a message with the elapsed time."""
    print(f"[{time.time()-T0:6.1f}s]", *a, flush=True)

A = ad.read_h5ad(S1A + "panel_object.h5ad")
o = A.obs
log(f"shape {A.shape}  nnz {A.X.nnz:,}  dtype {A.X.dtype}  min {A.X.data.min()}  max {A.X.data.max()}")
log(f"cell_id unique {o.index.is_unique}  n {len(o):,}")
log(f"var unique {A.var_names.is_unique}  n {A.n_vars}")

print("\n--- 1. per-arm cell counts (must equal the arm reports) ---")
print(o.source_key.value_counts().to_string())

print("\n--- 2. donors / patients / libraries / blocks by histology ---")
g = o.groupby("histology", observed=True)
tab = pd.DataFrame({
    "cells": g.size(),
    "donors": g.donor_id.nunique(),
    "patients": g.patient_id.nunique(),
    "libraries": g.library_id.nunique(),
    "blocks": g.block_id.nunique(),
    "studies": g.study.nunique(),
    "cells_composition_safe": g.composition_safe.sum(),
    "cells_labelled": g.apply(lambda d: (d.level1_compartment.astype(str) != "").sum(), include_groups=False),
}).sort_values("cells", ascending=False)
print(tab.to_string())
print(f"TOTAL cells {len(o):,}  donors {o.donor_id.nunique()}  patients {o.patient_id.nunique()} "
      f" libraries {o.library_id.nunique()}  studies {o.study.nunique()}  blocks {o.block_id.nunique()}")

print("\n--- 3. LUAD vs LUSC donor counts, and the composition-safe subset ---")
for h in ["LUAD", "LUSC"]:
    m = o.histology == h
    s = m & o.composition_safe
    print(f"  {h}: {o.loc[m,'donor_id'].nunique()} donors / {int(m.sum()):,} cells; "
          f"composition_safe {o.loc[s,'donor_id'].nunique()} donors / {int(s.sum()):,} cells")

print("\n--- 4. GENE COLUMN ALIGNMENT, checked biologically per arm ---")
def cp10k(mask, genes):
    """For each gene: detection percent and mean log1p of counts per 10,000 in the masked cells."""
    X = A[mask, :].X
    tot = np.asarray(X.sum(axis=1)).ravel(); tot[tot == 0] = 1
    out = {}
    for gname in genes:
        j = A.var_names.get_loc(gname)
        col = np.asarray(X[:, j].todense()).ravel()
        out[gname] = (100 * (col > 0).mean(), float(np.log1p(col / tot * 1e4).mean()))
    return out

markers = ["PTPRC", "EPCAM", "COL1A1", "RGS5", "CLDN5", "TPSAB1", "CD3E", "LYZ", "MS4A1"]
lab = (o.level1_compartment.astype(str) != "").values
for comp in ["T_NK", "myeloid", "epithelial_malignant", "B_plasma", "fibroblast", "mural", "endothelial", "other"]:
    m = (o.level1_compartment.astype(str) == comp).values
    r = cp10k(m, markers)
    print(f"  {comp:22s} n={int(m.sum()):7,d} " + "  ".join(f"{k} {v[1]:.2f}" for k, v in r.items()))

print("\n--- 5. the same markers, per arm, to prove the symbol join aligned columns ---")
for sk in ["LuCA_census", "EMTAB6149", "EMTAB6653"]:
    for comp in ["epithelial_malignant", "T_NK", "fibroblast"]:
        m = ((o.source_key.astype(str) == sk) & (o.level1_compartment.astype(str) == comp)).values
        if m.sum() == 0: continue
        r = cp10k(m, ["PTPRC", "EPCAM", "COL1A1", "TXNIP", "MRC1"])
        print(f"  {sk:12s} {comp:20s} n={int(m.sum()):7,d} " + "  ".join(f"{k} {v[0]:5.1f}pct" for k, v in r.items()))
# TXNIP and MRC1 are the two curated genes whose Ensembl id drifts in EMTAB6149.
print("  (TXNIP and MRC1 are two of the 26 genes with a drifted Ensembl id in EMTAB6149;")
print("   a non-zero detection rate there is the proof the symbol join recovered them.)")

print("\n--- 6. squamous / adeno identity by histology, epithelial_malignant only ---")
for h in ["LUAD", "LUSC"]:
    m = ((o.histology.astype(str) == h) & (o.level1_compartment.astype(str) == "epithelial_malignant")).values
    r = cp10k(m, ["KRT5", "KRT6A", "TP63", "NKX2-1", "NAPSA", "SFTPC"])
    print(f"  {h} n={int(m.sum()):7,d} " + "  ".join(f"{k} {v[0]:5.1f}pct" for k, v in r.items()))

print("\n--- 7. redox core genes: detection by arm ---")
for sk in ["LuCA_census", "EMTAB6149", "EMTAB6653"]:
    m = (o.source_key.astype(str) == sk).values
    r = cp10k(m, ["TXNRD1", "NFE2L2", "KEAP1", "TXN", "PRDX1", "NQO1", "TXNIP"])
    print(f"  {sk:12s} n={int(m.sum()):9,d} " + "  ".join(f"{k} {v[0]:5.1f}pct" for k, v in r.items()))

print("\n--- 8. spot check against the arm files: column sums on LuCA cells ---")
with h5py.File(S1A + "arm_luca.h5ad", "r") as h:
    lv = read_elem(h["var"])
sym = pd.Index(lv["panel_symbol"].astype(str))
Al = ad.read_h5ad(S1A + "arm_luca.h5ad")
Al.var_names = sym
probe = ["TXNRD1", "EPCAM", "COL1A1", "RGS5", "KRT5", "NKX2-1", "PTPRC", "NQO1"]
m = (o.source_key.astype(str) == "LuCA_census").values
Xm = A[m, :]
ok = True
for gname in probe:
    a = float(Al[:, gname].X.sum()); b = float(Xm[:, gname].X.sum())
    same = abs(a - b) < 1e-3 * max(1.0, abs(a))
    ok &= same
    print(f"  {gname:8s} arm {a:15,.0f}   merged {b:15,.0f}   {'PASS' if same else 'FAIL'}")
print(f"  column sums identical on all probes: {'PASS' if ok else 'FAIL'}")
del Al

print("\n--- 9. schema invariants ---")
print(f"  composition_safe only where whole_tissue: "
      f"{bool((o.composition_safe.values == (o.enrichment.astype(str)=='whole_tissue').values).all())}")
print(f"  counts_are_estimates True on {int(o.counts_are_estimates.sum()):,} cells, "
      f"studies {sorted(o.loc[o.counts_are_estimates,'study'].astype(str).unique())}")
print(f"  platform levels: {sorted(o.platform.astype(str).unique())}")
print(f"  modality levels: {sorted(o.modality.astype(str).unique())}")
print(f"  enrichment levels: {dict(o.enrichment.value_counts())}")
print(f"  treatment levels: {dict(o.treatment.value_counts())}")
print(f"  tissue_class levels: {dict(o.tissue_class.value_counts())}")
print(f"  histology_source: {dict(o.histology_source.value_counts())}")
print(f"  histology_confidence: {dict(o.histology_confidence.value_counts())}")
print(f"  celltype_source: {dict(o.celltype_source.value_counts())}")
print(f"  level2_celltype: {sorted(o.level2_celltype.astype(str).unique())}")
print(f"  sampling_site: {dict(o.sampling_site.value_counts())}")
print(f"  blocks: {o.block_id.nunique()}")
print(f"  donor_id relabelled for EMTAB6653: "
      f"{sorted(o.loc[o.source_key.astype(str)=='EMTAB6653','donor_id'].astype(str).unique())} "
      f"was {sorted(o.loc[o.source_key.astype(str)=='EMTAB6653','donor_id_native'].astype(str).unique())}")
print(f"  Lambrechts donors overall: {sorted(o.loc[o.study.astype(str)=='Lambrechts_Thienpont_2018','donor_id'].astype(str).unique())}")

print("\n--- 10. depth by platform (n_counts over each arm's own full gene space) ---")
print(o.groupby("platform", observed=True)[["n_counts", "n_genes", "n_counts_intersect", "n_genes_intersect"]]
      .median().astype(int).to_string())
log("verify done")
