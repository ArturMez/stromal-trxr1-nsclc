#!/usr/bin/env python3
"""
Atlas stage 04, step extract_luca: extract the LuCA cells at all 17,811 genes.

Paper items: part of the full-gene atlas (fullgene_object.h5ad) behind Figure 2, Figure 3,
Supplementary Figures 2 to 4 and Supplementary Data 1.

LuCA is the backbone: 13 of the 14 studies in `panel_object.h5ad`, and 771,094 of the 1,587,067
fibroblast, myeloid and T/NK cells the full-gene object targets. Every one of those cells joins to
LuCA by barcode, so nothing is lost here. No random numbers are drawn.

WHAT THIS TAKES AND WHAT IT LEAVES. It takes the raw counts from `layers['count']`, because `X` in that
file is normalized and an object that mixes the two is silently wrong. It takes all 17,811 genes rather
than the 1,906-gene panel, which is the entire point of the rebuild. It leaves behind every compartment
other than the three targeted, and it leaves the LuCA annotations behind as well: the harmonized
`level1_compartment`, `level2_celltype`, `block_qc_pass`, `histology`, `tissue_class`, `donor_id` and
`block_id` come from `panel_object.h5ad`, so the rebuilt object stays answerable to the same QC and the
same block structure as the panel object.

WHY IT STREAMS. The output is roughly 900 million nonzeros, about 7 GB. Holding that in memory while
also holding the source is not practical, so rows are read in contiguous chunks, filtered
inside the chunk, and appended to resizable datasets in the output file. Memory stays flat at one chunk.

GENE IDENTITY. LuCA's `var` is indexed by Ensembl id with symbols in `var/feature_name`. Symbols are
what every gene set in the analysis is written in, so symbols are the axis. Where a symbol repeats,
which happens for a handful of Ensembl ids, the duplicate columns are summed rather than dropped, and
the count of affected symbols is reported.

RUN: python3 atlas/steps/fullbuild_01_extract_luca.py   (normally through atlas/04_build_objects.py)
FROM: LuCA/luca_extended.h5ad (raw counts in layers['count']) and stage1a/panel_object.h5ad (obs).
OUT -> output/RNA/meta_dataset/fullbuild/arm_luca_full.h5ad
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
SRC = f"{SCRNA}/LuCA/luca_extended.h5ad"
OUT = f"{ROOT}/output/RNA/meta_dataset/fullbuild"
DST = f"{OUT}/arm_luca_full.h5ad"
TARGET = ["fibroblast", "myeloid", "T_NK"]
CHUNK = 20000
# counts_are_estimates is NOT optional: the 17,891 Smart-seq2 cells in this arm carry non-integer
# values in layers['count'], because Smart-seq2 counts are length-scaled estimates, and they agree with
# this flag on 100 per cent of a 400-cell sample. An object that mixes them with true counts without a
# flag is silently wrong (the same reason E-MTAB-6149 is flagged).
CARRY = ["level1_compartment", "level2_celltype", "level2_celltype_native", "level3_celltype_native",
         "donor_id", "patient_id", "study", "block_id", "block_qc_pass", "histology", "tissue_class",
         "platform", "sampling_site", "n_counts", "n_genes", "source_key", "is_malignant_basis",
         "sex", "age", "smoking", "stage", "counts_are_estimates", "enrichment", "doublet_score",
         "celltype_confidence", "celltype_source", "composition_safe"]
os.makedirs(OUT, exist_ok=True)


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
    """Write the anndata `column-order` attribute. An EMPTY numpy object array has no HDF5
    equivalent, so the empty case is handled explicitly."""
    if cols:
        grp.attrs["column-order"] = np.array(list(cols), dtype=h5py.string_dtype())
    else:
        grp.attrs.create("column-order", np.array([], dtype="S1"))


def x_is_complete(path, n_rows):
    """True if a previous run already streamed the whole matrix, so it can be reused."""
    if not os.path.exists(path):
        return False
    try:
        with h5py.File(path, "r") as h:
            ip = h["X/indptr"]
            return len(ip) == n_rows + 1 and int(ip[-1]) == h["X/data"].shape[0]
    except Exception:
        return False


print("reading panel_object obs ...", flush=True)
with h5py.File(PANEL, "r") as h:
    p_idx = dec(h["obs/_index"][:])
    P = {k: col(h, f"obs/{k}") for k in CARRY if k in h["obs"]}
P = pd.DataFrame(P)
P["cell_id"] = p_idx
P["barcode"] = pd.Series(p_idx).str.split("|").str[-1]
sel = P[(P.level1_compartment.isin(TARGET)) & (P.source_key == "LuCA_census")].copy()
print(f"  target cells from LuCA: {len(sel):,}", flush=True)

print("reading LuCA var and index ...", flush=True)
with h5py.File(SRC, "r") as h:
    sym = col(h, "var/feature_name")
    src_bc = col(h, "obs/_index")
    n_src, n_gene = (int(x) for x in h["X"].attrs["shape"])
pos = pd.Series(np.arange(n_src), index=src_bc)
pos = pos[~pos.index.duplicated()]
sel["row"] = sel.barcode.map(pos)
miss = int(sel.row.isna().sum())
if miss:
    print(f"  WARNING {miss:,} cells did not map and are dropped", flush=True)
    sel = sel[sel.row.notna()]
sel["row"] = sel.row.astype(int)
sel = sel.sort_values("row").reset_index(drop=True)
want = np.zeros(n_src, bool)
want[sel.row.values] = True

uniq, inv = np.unique(sym, return_inverse=True)
n_dup = len(sym) - len(uniq)
print(f"  genes {n_gene:,} -> {len(uniq):,} unique symbols ({n_dup} duplicate columns summed)",
      flush=True)
collapse = n_dup > 0


def stream_counts():
    """Read the raw counts for the selected rows in contiguous chunks and append them to the
    output file. Memory stays flat at one chunk regardless of how large the result is."""
    with h5py.File(SRC, "r") as h, h5py.File(DST, "w") as o:
        indptr = h["layers/count/indptr"][:]
        D, I = h["layers/count/data"], h["layers/count/indices"]
        g = o.create_group("X")
        d_out = g.create_dataset("data", (0,), maxshape=(None,), dtype="float32", chunks=(1 << 20,))
        i_out = g.create_dataset("indices", (0,), maxshape=(None,), dtype="int32", chunks=(1 << 20,))
        ptr, nnz, done = [0], 0, 0
        for a in range(0, n_src, CHUNK):
            b = min(a + CHUNK, n_src)
            if not want[a:b].any():
                continue
            lo, hi = int(indptr[a]), int(indptr[b])
            dat, ind = D[lo:hi], I[lo:hi]
            loc = indptr[a:b + 1] - lo
            cd, ci = [], []
            for r in np.where(want[a:b])[0]:
                s0, s1 = int(loc[r]), int(loc[r + 1])
                gi, vv = ind[s0:s1], dat[s0:s1]
                if collapse:
                    gi2 = inv[gi]
                    order = np.argsort(gi2, kind="stable")
                    gi2, vv = gi2[order], vv[order]
                    uu, start = np.unique(gi2, return_index=True)
                    vv = np.add.reduceat(vv, start) if len(uu) else vv
                    gi = uu
                else:
                    gi = inv[gi]
                cd.append(vv)
                ci.append(gi)
                ptr.append(ptr[-1] + len(gi))
                done += 1
            cd = np.concatenate(cd).astype("float32")
            ci = np.concatenate(ci).astype("int32")
            d_out.resize((nnz + len(cd),))
            i_out.resize((nnz + len(ci),))
            d_out[nnz:] = cd
            i_out[nnz:] = ci
            nnz += len(cd)
            if done % 100000 < CHUNK:
                print(f"    {done:,} / {len(sel):,} cells, {nnz:,} nonzeros", flush=True)
        g.create_dataset("indptr", data=np.asarray(ptr, "int64"))
        g.attrs["encoding-type"] = "csr_matrix"
        g.attrs["encoding-version"] = "0.1.0"
        g.attrs["shape"] = np.array([len(sel), len(uniq)], "int64")
    return nnz


def write_metadata():
    """var and obs. Separate from the stream so a metadata bug never costs the matrix again."""
    with h5py.File(DST, "r+") as o:
        for gname in ("var", "obs"):
            if gname in o:
                del o[gname]
        v = o.create_group("var")
        v.attrs["_index"] = "gene_symbol"
        v.attrs["encoding-type"] = "dataframe"
        v.attrs["encoding-version"] = "0.2.0"
        set_colorder(v, [])
        v.create_dataset("gene_symbol", data=np.array(uniq, dtype=object),
                         dtype=h5py.string_dtype())
        ob = o.create_group("obs")
        ob.attrs["_index"] = "_index"
        ob.attrs["encoding-type"] = "dataframe"
        ob.attrs["encoding-version"] = "0.2.0"
        cols = [c for c in CARRY if c in sel.columns]
        set_colorder(ob, cols)
        ob.create_dataset("_index", data=np.array(sel.cell_id.values, dtype=object),
                          dtype=h5py.string_dtype())
        for c in cols:
            x = sel[c].values
            if x.dtype == object or str(x.dtype).startswith("str"):
                ob.create_dataset(c, data=np.array(x.astype(str), dtype=object),
                                  dtype=h5py.string_dtype())
            else:
                ob.create_dataset(c, data=x)
        o.attrs["encoding-type"] = "anndata"
        o.attrs["encoding-version"] = "0.1.0"


if x_is_complete(DST, len(sel)):
    with h5py.File(DST, "r") as o:
        nnz = int(o["X/data"].shape[0])
    print(f"  X already complete, reusing it and writing metadata only ({nnz:,} nonzeros)", flush=True)
else:
    print("streaming counts ...", flush=True)
    nnz = stream_counts()
write_metadata()

print(f"\nwritten -> {DST}")
print(f"  cells {len(sel):,}   genes {len(uniq):,}   nonzeros {nnz:,}   "
      f"{os.path.getsize(DST)/1e9:.1f} GB on disk")
print(f"  compartments: {sel.level1_compartment.value_counts().to_dict()}")
print(f"  donors {sel.donor_id.nunique()}   studies {sel.study.nunique()}")
