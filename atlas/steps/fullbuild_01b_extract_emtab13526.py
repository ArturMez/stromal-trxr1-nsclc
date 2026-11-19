#!/usr/bin/env python3
"""
Atlas stage 04, step extract_emtab13526: extract the De Zuani cells at all 33,514 genes.

Paper items: part of the full-gene atlas (fullgene_object.h5ad) behind Figure 2, Figure 3,
Supplementary Figures 2 to 4 and Supplementary Data 1.

E-MTAB-13526 is the second backbone: 756,848 of the 1,587,067 fibroblast, myeloid and T/NK cells the
full-gene object targets, in 81 CellRanger MTX libraries; all of them join by library and barcode. Panel ids
read `EMTAB13526|<library_id>|<barcode>`, and the 81 library names on disk match the 81 in the object
exactly, so the library is the join key and the barcode is unique inside it.

WHY THIS SOURCE IS WORTH THE TROUBLE. It carries 33,514 measured genes against LuCA's 17,811, so it is
what actually widens the gene axis. The genes LuCA never measured are structural zeros for LuCA cells
and real measurements here, which is why the merge step writes a per-source presence mask
rather than pretending the axis is shared.

FEATURES ARE ASSERTED IDENTICAL ACROSS LIBRARIES, not assumed. If any library deposits a different
feature list the script stops, because silently taking the intersection or the first list would
misalign every gene index in the arm.

MEMORY. One library at a time: its triplets are read, filtered to the wanted barcodes, converted to
CSR rows and appended to resizable datasets. Peak memory is one library, not one arm.

DUPLICATE SYMBOLS are summed rather than dropped, and the number affected is reported. CellRanger
references carry a few hundred symbols on more than one Ensembl id and dropping them loses counts.

RUN: python3 atlas/steps/fullbuild_01b_extract_emtab13526.py   (normally through atlas/04_build_objects.py)
FROM: E-MTAB-13526/raw_counts (81 libraries) and stage1a/panel_object.h5ad (obs). No random numbers are drawn.
OUT -> output/RNA/meta_dataset/fullbuild/arm_emtab13526_full.h5ad
"""
import os
import gzip
import numpy as np
import pandas as pd
import h5py

import os as _os  # noqa: E402
import sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # atlas/
from common import ROOT, SCRNA, SCRATCH, DATA, RSCRIPT  # noqa: E402,F401  paths from config.json
PANEL = f"{ROOT}/output/RNA/meta_dataset/stage1a/panel_object.h5ad"
RAW = f"{SCRNA}/E-MTAB-13526/raw_counts"
OUT = f"{ROOT}/output/RNA/meta_dataset/fullbuild"
DST = f"{OUT}/arm_emtab13526_full.h5ad"
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
    """Write the anndata  attribute (an empty list needs an explicit empty array)."""
    if cols:
        grp.attrs["column-order"] = np.array(list(cols), dtype=h5py.string_dtype())
    else:
        grp.attrs.create("column-order", np.array([], dtype="S1"))


def read_mtx(path):
    """Read a CellRanger MatrixMarket file as (rows, cols, vals, shape). Rows are features."""
    with gzip.open(path, "rt") as f:
        line = f.readline()
        while True:
            line = f.readline()
            if not line.startswith("%"):
                break
        nr, nc, nnz = (int(x) for x in line.split())
        t = pd.read_csv(f, sep=r"\s+", header=None, dtype=np.int64, engine="c").values
    return t[:, 0] - 1, t[:, 1] - 1, t[:, 2].astype("float32"), (nr, nc)


print("reading panel_object obs ...", flush=True)
with h5py.File(PANEL, "r") as h:
    p_idx = dec(h["obs/_index"][:])
    P = pd.DataFrame({k: col(h, f"obs/{k}") for k in CARRY if k in h["obs"]})
P["cell_id"] = p_idx
P["barcode"] = pd.Series(p_idx).str.split("|").str[-1]
sel = P[(P.level1_compartment.isin(TARGET)) & (P.source_key == "EMTAB13526")].copy()
libs = sorted(sel.library_id.unique())
print(f"  target cells {len(sel):,} across {len(libs)} libraries", flush=True)

feat0 = None
for L in libs:
    fp = f"{RAW}/{L}-features.tsv.gz"
    if not os.path.exists(fp):
        raise SystemExit(f"missing features for library {L}")
    f = pd.read_csv(fp, sep="\t", header=None, compression="gzip")
    sy = f[1].astype(str).values if f.shape[1] > 1 else f[0].astype(str).values
    if feat0 is None:
        feat0 = sy
    elif not np.array_equal(feat0, sy):
        raise SystemExit(f"library {L} deposits a DIFFERENT feature list; stopping rather than guessing")
print(f"  features identical across all {len(libs)} libraries: {len(feat0):,}", flush=True)

uniq, inv = np.unique(feat0, return_inverse=True)
n_dup = len(feat0) - len(uniq)
print(f"  genes {len(feat0):,} -> {len(uniq):,} unique symbols ({n_dup} duplicate columns summed)",
      flush=True)

order, ptr, nnz = [], [0], 0
with h5py.File(DST, "w") as o:
    g = o.create_group("X")
    d_out = g.create_dataset("data", (0,), maxshape=(None,), dtype="float32", chunks=(1 << 20,))
    i_out = g.create_dataset("indices", (0,), maxshape=(None,), dtype="int32", chunks=(1 << 20,))
    for k, L in enumerate(libs):
        want = sel[sel.library_id == L]
        bc = pd.read_csv(f"{RAW}/{L}-barcodes.tsv.gz", sep="\t", header=None,
                         compression="gzip")[0].astype(str).values
        pos = pd.Series(np.arange(len(bc)), index=bc)
        pos = pos[~pos.index.duplicated()]
        cix = want.barcode.map(pos)
        if cix.isna().any():
            raise SystemExit(f"library {L}: {int(cix.isna().sum())} barcodes not found")
        cix = cix.astype(int).values
        r, c, v, shape = read_mtx(f"{RAW}/{L}-matrix.mtx.gz")
        if shape[0] != len(feat0):
            raise SystemExit(f"library {L}: matrix has {shape[0]} rows, features has {len(feat0)}")
        newcol = np.full(shape[1], -1, np.int64)
        newcol[cix] = np.arange(len(cix))
        keep = newcol[c] >= 0
        r, v, cc = r[keep], v[keep], newcol[c[keep]]
        gidx = inv[r]
        srt = np.lexsort((gidx, cc))
        cc, gidx, v = cc[srt], gidx[srt], v[srt]
        counts = np.bincount(cc, minlength=len(cix))
        starts = np.concatenate([[0], np.cumsum(counts)])
        cd, ci, cp = [], [], []
        for j in range(len(cix)):
            a, b = int(starts[j]), int(starts[j + 1])
            gg, vv = gidx[a:b], v[a:b]
            if n_dup and len(gg):
                uu, st = np.unique(gg, return_index=True)
                vv = np.add.reduceat(vv, st)
                gg = uu
            cd.append(vv)
            ci.append(gg)
            cp.append(len(gg))
        cd = np.concatenate(cd).astype("float32") if cd else np.array([], "float32")
        ci = np.concatenate(ci).astype("int32") if ci else np.array([], "int32")
        d_out.resize((nnz + len(cd),))
        i_out.resize((nnz + len(ci),))
        d_out[nnz:] = cd
        i_out[nnz:] = ci
        nnz += len(cd)
        for n_ in cp:
            ptr.append(ptr[-1] + n_)
        order.append(want)
        print(f"    [{k+1}/{len(libs)}] {L}: {len(cix):,} cells, running nonzeros {nnz:,}", flush=True)
    sel2 = pd.concat(order).reset_index(drop=True)
    g.create_dataset("indptr", data=np.asarray(ptr, "int64"))
    g.attrs["encoding-type"] = "csr_matrix"
    g.attrs["encoding-version"] = "0.1.0"
    g.attrs["shape"] = np.array([len(sel2), len(uniq)], "int64")

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
    cols = [c for c in CARRY if c in sel2.columns]
    set_colorder(ob, cols)
    ob.create_dataset("_index", data=np.array(sel2.cell_id.values, dtype=object),
                      dtype=h5py.string_dtype())
    for c in cols:
        x = sel2[c].values
        if x.dtype == object or str(x.dtype).startswith("str"):
            ob.create_dataset(c, data=np.array(x.astype(str), dtype=object), dtype=h5py.string_dtype())
        else:
            ob.create_dataset(c, data=x)
    o.attrs["encoding-type"] = "anndata"
    o.attrs["encoding-version"] = "0.1.0"

print(f"\nwritten -> {DST}")
print(f"  cells {len(sel2):,}   genes {len(uniq):,}   nonzeros {nnz:,}   "
      f"{os.path.getsize(DST)/1e9:.1f} GB on disk")
print(f"  compartments: {sel2.level1_compartment.value_counts().to_dict()}")
print(f"  donors {sel2.donor_id.nunique()}   libraries {sel2.library_id.nunique()}")
