#!/usr/bin/env python3
"""
Atlas stage 04, step extract_emtab6149: assemble the E-MTAB-6149 arm from the recovered counts.

Paper items: part of the full-gene atlas (fullgene_object.h5ad) behind Figure 2, Figure 3,
Supplementary Figures 2 to 4 and Supplementary Data 1. No random numbers are drawn.

WHY THIS ARM ENTERS A COUNTS-BASED OBJECT. The deposited values look like normalized expression, but
they are exactly invertible to integer counts, so the arm carries no estimate flag:

  1  The deposited values are `ln(1 + c * 1e4 / T)`, which is exactly invertible: within one cell the
     smallest positive value is c = 1, so the scale factor and then every integer count follow.
     `stage1a_export6149.R` does this and reports the deviation rather than asserting it. Run across
     ALL 22,180 genes rather than the 1,906-gene panel:
     **integers recovered for 52,698 of 52,698 cells, worst deviation 3.274e-11, minimum 1,
     maximum 28,114.** These are counts.
  2  The accession that genuinely carries non-integer values is E-MTAB-6653, whose cells are Alevin EM
     estimates and are flagged.

  The five patients `Lambrechts_Thienpont_2018_P1` to `P5` contribute 35,962 target cells (21,145 LUAD,
  11,732 large cell and 3,085 LUSC), and none of them reaches the object through any other source.

THE JOIN. The R export names cells `<barcode>_<library_index>`, and panel ids read
`EMTAB6149|<library_id>|<barcode>` where `library_id` ends in `_L<index>`. So the key is the barcode
plus the library index, and the script REPORTS the join rate and stops if it is not total, because a
partial join here silently drops patients.

RUN: python3 atlas/steps/fullbuild_01d_extract_emtab6149.py   (normally through atlas/04_build_objects.py)
     (after the R export has written panel_counts.mtx into fullbuild/e6149/)
OUT -> output/RNA/meta_dataset/fullbuild/arm_emtab6149_full.h5ad
"""
import os
import re
import numpy as np
import pandas as pd
import h5py

import os as _os  # noqa: E402
import sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # atlas/
from common import ROOT, SCRNA, SCRATCH, DATA, RSCRIPT  # noqa: E402,F401  paths from config.json
PANEL = f"{ROOT}/output/RNA/meta_dataset/stage1a/panel_object.h5ad"
E = f"{ROOT}/output/RNA/meta_dataset/fullbuild/e6149"
DST = f"{ROOT}/output/RNA/meta_dataset/fullbuild/arm_emtab6149_full.h5ad"
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


for f in ("panel_counts.mtx", "panel_genes_present.txt", "panel_cells.txt", "matrix_features.tsv"):
    if not os.path.exists(f"{E}/{f}"):
        raise SystemExit(f"the R export has not written {f} yet")

print("reading panel_object obs ...", flush=True)
with h5py.File(PANEL, "r") as h:
    p_idx = dec(h["obs/_index"][:])
    P = pd.DataFrame({k: col(h, f"obs/{k}") for k in CARRY if k in h["obs"]})
P["cell_id"] = p_idx
parts = pd.Series(p_idx).str.split("|")
P["bc"] = parts.str[2]
P["lidx"] = P.library_id.astype(str).str.extract(r"_L0*(\d+)$")[0]
P["key"] = P.bc + "_" + P.lidx
sel = P[(P.level1_compartment.isin(TARGET)) & (P.source_key == "EMTAB6149")].copy()
print(f"  target cells {len(sel):,} across {sel.patient_id.nunique()} patients", flush=True)

cells = pd.read_csv(f"{E}/panel_cells.txt", header=None)[0].astype(str).values
ens = pd.read_csv(f"{E}/panel_genes_present.txt", header=None)[0].astype(str).values
feat = pd.read_csv(f"{E}/matrix_features.tsv", sep="\t")
symcol = [c for c in feat.columns if "gene" in c.lower() or "symbol" in c.lower() or "name" in c.lower()]
idcol = [c for c in feat.columns if "ensembl" in c.lower() or c.lower() in ("id", "gene_id")]
m = {}
if symcol and idcol:
    m = dict(zip(feat[idcol[0]].astype(str), feat[symcol[0]].astype(str)))
sym = np.array([m.get(e, e) if isinstance(m.get(e, e), str) and m.get(e, e) not in ("nan", "") else e
                for e in ens])
n_unres = int((sym == ens).sum())
print(f"  Ensembl ids resolved to a symbol: {len(ens)-n_unres:,}; kept as Ensembl: {n_unres:,}",
      flush=True)
uniq, inv = np.unique(sym, return_inverse=True)
dup = len(ens) - len(uniq)
print(f"  genes {len(ens):,} -> {len(uniq):,} unique symbols ({dup} duplicate columns summed)",
      flush=True)

pos = pd.Series(np.arange(len(cells)), index=cells)
pos = pos[~pos.index.duplicated()]
sel["cix"] = sel.key.map(pos)
rate = 100 * float(sel.cix.notna().mean())
print(f"  JOIN RATE panel cells found in the export: {rate:.2f}%", flush=True)
if rate < 100:
    bad = sel[sel.cix.isna()]
    raise SystemExit(f"join incomplete, {len(bad):,} cells unmatched across "
                     f"{bad.patient_id.nunique()} patients; refusing to drop patients silently")
sel = sel.sort_values("cix").reset_index(drop=True)
cix = sel.cix.astype(int).values
newcol = np.full(len(cells), -1, np.int64)
newcol[cix] = np.arange(len(cix))

print("reading the recovered count matrix ...", flush=True)
R, C, V = [], [], []
with open(f"{E}/panel_counts.mtx") as f:
    f.readline()
    while True:
        line = f.readline()
        if not line.startswith("%"):
            break
    for chunk in pd.read_csv(f, sep=r"\s+", header=None, chunksize=10_000_000, engine="c"):
        t = chunk.values
        r = t[:, 0].astype(np.int64) - 1
        c = t[:, 1].astype(np.int64) - 1
        v = t[:, 2].astype("float32")
        keep = newcol[c] >= 0
        R.append(inv[r[keep]])
        C.append(newcol[c[keep]])
        V.append(v[keep])
        print(f"    kept {sum(len(x) for x in V):,} nonzeros", flush=True)
r = np.concatenate(R)
c = np.concatenate(C)
v = np.concatenate(V)
print(f"  values integral: {np.allclose(v, np.round(v))}  (expected True, counts were recovered)",
      flush=True)

srt = np.lexsort((r, c))
c, r, v = c[srt], r[srt], v[srt]
counts = np.bincount(c, minlength=len(cix))
starts = np.concatenate([[0], np.cumsum(counts)])
cd, ci, ptr = [], [], [0]
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
print(f"  patients {sel.patient_id.nunique()}: {sorted(sel.patient_id.unique())}")
print(f"  histology: {sel.histology.value_counts().to_dict()}")
