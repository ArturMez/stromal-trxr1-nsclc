#!/usr/bin/env python3
"""
Atlas stage 04, step extract_emtab6653: extract the E-MTAB-6653 cells at all measured genes.

Paper items: part of the full-gene atlas (fullgene_object.h5ad) behind Figure 2, Figure 3,
Supplementary Figures 2 to 4 and Supplementary Data 1.

23,163 of the target fibroblast, myeloid and T/NK cells, from one aggregated MatrixMarket of
23,979 Ensembl ids by 32,341 cells. No random numbers are drawn.

THE FLAG THAT TRAVELS WITH THIS ARM. These values are **Alevin EM estimates and are not integer UMI
counts**, which `stage1a_arm_emtab6653.py` establishes and which this script re-checks rather than
trusts. Every cell here carries `counts_are_estimates = True`.

IDENTIFIERS. Rows are unversioned Ensembl ids, and every gene set in the analysis is written in
symbols, so the ids are mapped through the E-MTAB-13526 CellRanger reference, which carries both for
33,538 genes, with LuCA's `var` as a second pass. Ensembl ids that resolve to no symbol are KEPT under
their Ensembl id rather than dropped, so no counts are silently lost and the merge step can still
align whatever does resolve. Columns read `<run>-<barcode>` and panel ids read
`EMTAB6653|<run>|<barcode>`, so the join key is the pair.

RUN: python3 atlas/steps/fullbuild_01c_extract_emtab6653.py   (normally through atlas/04_build_objects.py)
FROM: E-MTAB-6653/raw_counts/matrix_market, the E-MTAB-13526 CellRanger features, LuCA var, and
      stage1a/panel_object.h5ad (obs).
OUT -> output/RNA/meta_dataset/fullbuild/arm_emtab6653_full.h5ad
"""
import os
import numpy as np
import pandas as pd
import h5py

import os as _os  # noqa: E402
import sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # atlas/
from common import ROOT, SCRNA, SCRATCH, DATA, RSCRIPT  # noqa: E402,F401  paths from config.json
PANEL = f"{ROOT}/output/RNA/meta_dataset/stage1a/panel_object.h5ad"
MM = f"{SCRNA}/E-MTAB-6653/raw_counts/matrix_market"
REF = f"{SCRNA}/E-MTAB-13526/raw_counts/D1_1-features.tsv.gz"
LUCA = f"{SCRNA}/LuCA/luca_extended.h5ad"
OUT = f"{ROOT}/output/RNA/meta_dataset/fullbuild"
DST = f"{OUT}/arm_emtab6653_full.h5ad"
TARGET = ["fibroblast", "myeloid", "T_NK"]
CARRY = ["level1_compartment", "level2_celltype", "level2_celltype_native", "level3_celltype_native",
         "donor_id", "patient_id", "study", "block_id", "block_qc_pass", "histology", "tissue_class",
         "platform", "sampling_site", "n_counts", "n_genes", "source_key", "is_malignant_basis",
         "sex", "age", "smoking", "stage", "counts_are_estimates", "enrichment", "doublet_score",
         "celltype_confidence", "celltype_source", "composition_safe", "library_id"]


def dec(a):
    """Decode an array of bytes or strings to an array of str."""
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in a])


def col(h, path):
    """Read one dataset or categorical obs column of an h5ad file as an array of values."""
    g = h[path]
    if isinstance(g, h5py.Group) and "categories" in g:
        return dec(g["categories"][:])[g["codes"][:]]
    v = g[:]
    return dec(v) if v.dtype.kind in "SO" else v


def set_colorder(grp, cols):
    """Write the anndata `column-order` attribute (an empty list needs an explicit empty array)."""
    if cols:
        grp.attrs["column-order"] = np.array(list(cols), dtype=h5py.string_dtype())
    else:
        grp.attrs.create("column-order", np.array([], dtype="S1"))


print("reading panel_object obs ...", flush=True)
with h5py.File(PANEL, "r") as h:
    p_idx = dec(h["obs/_index"][:])
    P = pd.DataFrame({k: col(h, f"obs/{k}") for k in CARRY if k in h["obs"]})
P["cell_id"] = p_idx
parts = pd.Series(p_idx).str.split("|")
P["key"] = parts.str[1] + "-" + parts.str[2]
sel = P[(P.level1_compartment.isin(TARGET)) & (P.source_key == "EMTAB6653")].copy()
print(f"  target cells {len(sel):,}", flush=True)

rows = pd.read_csv(f"{MM}/E-MTAB-6653.aggregated_filtered_counts.mtx_rows", sep="\t", header=None)
cols_ = pd.read_csv(f"{MM}/E-MTAB-6653.aggregated_filtered_counts.mtx_cols", sep="\t", header=None)
ens = rows[0].astype(str).values
cellkey = cols_[0].astype(str).values
print(f"  matrix ids: {len(ens):,} genes x {len(cellkey):,} cells", flush=True)

ref = pd.read_csv(REF, sep="\t", header=None, compression="gzip")
m1 = dict(zip(ref[0].astype(str), ref[1].astype(str)))
with h5py.File(LUCA, "r") as h:
    m2 = dict(zip(col(h, "var/gene_id"), col(h, "var/feature_name")))
sym = np.array([m1.get(e, m2.get(e, e)) for e in ens])
n_unres = int(sum(1 for e, s in zip(ens, sym) if s == e))
print(f"  Ensembl ids resolved to a symbol: {len(ens)-n_unres:,}; kept under their Ensembl id: {n_unres:,}",
      flush=True)
uniq, inv = np.unique(sym, return_inverse=True)
print(f"  genes {len(ens):,} -> {len(uniq):,} unique symbols "
      f"({len(ens)-len(uniq)} duplicate columns summed)", flush=True)

pos = pd.Series(np.arange(len(cellkey)), index=cellkey)
pos = pos[~pos.index.duplicated()]
sel["cix"] = sel.key.map(pos)
if sel.cix.isna().any():
    raise SystemExit(f"{int(sel.cix.isna().sum())} target cells not found in mtx_cols")
sel = sel.sort_values("cix").reset_index(drop=True)
cix = sel.cix.astype(int).values
newcol = np.full(len(cellkey), -1, np.int64)
newcol[cix] = np.arange(len(cix))

print("reading the matrix in chunks ...", flush=True)
R, C, V = [], [], []
kept = 0
with open(f"{MM}/E-MTAB-6653.aggregated_filtered_counts.mtx") as f:
    line = f.readline()
    while True:
        line = f.readline()
        if not line.startswith("%"):
            break
    for chunk in pd.read_csv(f, sep=r"\s+", header=None, chunksize=5_000_000, engine="c"):
        t = chunk.values
        r = t[:, 0].astype(np.int64) - 1
        c = t[:, 1].astype(np.int64) - 1
        v = t[:, 2].astype("float32")
        keep = newcol[c] >= 0
        R.append(inv[r[keep]])
        C.append(newcol[c[keep]])
        V.append(v[keep])
        kept += int(keep.sum())
        print(f"    {kept:,} nonzeros kept", flush=True)
r = np.concatenate(R)
c = np.concatenate(C)
v = np.concatenate(V)
print(f"  values integral: {np.allclose(v, np.round(v))}  (expected False, Alevin EM estimates)",
      flush=True)

srt = np.lexsort((r, c))
c, r, v = c[srt], r[srt], v[srt]
counts = np.bincount(c, minlength=len(cix))
starts = np.concatenate([[0], np.cumsum(counts)])
cd, ci, ptr = [], [], [0]
dup = len(ens) - len(uniq)
for j in range(len(cix)):
    a, b = int(starts[j]), int(starts[j + 1])
    gg, vv = r[a:b], v[a:b]
    if dup and len(gg):
        uu, st = np.unique(gg, return_index=True)
        vv = np.add.reduceat(vv, st)
        gg = uu
    cd.append(vv)
    ci.append(gg)
    ptr.append(ptr[-1] + len(gg))
cd = np.concatenate(cd).astype("float32")
ci = np.concatenate(ci).astype("int32")

with h5py.File(DST, "w") as o:
    g = o.create_group("X")
    g.create_dataset("data", data=cd, chunks=(min(len(cd), 1 << 20),))
    g.create_dataset("indices", data=ci, chunks=(min(len(ci), 1 << 20),))
    g.create_dataset("indptr", data=np.asarray(ptr, "int64"))
    g.attrs["encoding-type"] = "csr_matrix"
    g.attrs["encoding-version"] = "0.1.0"
    g.attrs["shape"] = np.array([len(sel), len(uniq)], "int64")
    v_ = o.create_group("var")
    v_.attrs["_index"] = "gene_symbol"
    v_.attrs["encoding-type"] = "dataframe"
    v_.attrs["encoding-version"] = "0.2.0"
    set_colorder(v_, [])
    v_.create_dataset("gene_symbol", data=np.array(uniq, dtype=object), dtype=h5py.string_dtype())
    ob = o.create_group("obs")
    ob.attrs["_index"] = "_index"
    ob.attrs["encoding-type"] = "dataframe"
    ob.attrs["encoding-version"] = "0.2.0"
    cc = [c_ for c_ in CARRY if c_ in sel.columns]
    set_colorder(ob, cc)
    ob.create_dataset("_index", data=np.array(sel.cell_id.values, dtype=object),
                      dtype=h5py.string_dtype())
    for c_ in cc:
        x = sel[c_].values
        if x.dtype == object or str(x.dtype).startswith("str"):
            ob.create_dataset(c_, data=np.array(x.astype(str), dtype=object),
                              dtype=h5py.string_dtype())
        else:
            ob.create_dataset(c_, data=x)
    o.attrs["encoding-type"] = "anndata"
    o.attrs["encoding-version"] = "0.1.0"

print(f"\nwritten -> {DST}")
print(f"  cells {len(sel):,}   genes {len(uniq):,}   nonzeros {len(cd):,}   "
      f"{os.path.getsize(DST)/1e9:.2f} GB on disk")
print(f"  compartments: {sel.level1_compartment.value_counts().to_dict()}")
print(f"  counts_are_estimates True: {int(pd.Series(sel.counts_are_estimates).astype(str).isin(['True','1']).sum()):,}")
