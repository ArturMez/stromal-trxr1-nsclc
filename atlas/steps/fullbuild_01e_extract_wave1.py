#!/usr/bin/env python3
"""
Atlas stage 04, step extract_wave1: the three wave-1 arms at all measured genes.

Paper items: part of the full-gene atlas (fullgene_object.h5ad) behind Figure 2, Figure 3,
Supplementary Figures 2 to 4 and Supplementary Data 1. No random numbers are drawn.

STTT, GSE325414 and GSE246215, re-extracted from their raw deposits at full transcriptome. Their
harmonized annotations already exist in `output/RNA/meta_dataset/wave1/arm_*.h5ad`, built at the
1,906-gene panel (stage 03), so this script takes the OBS from there and only the COUNTS from the
deposits, joining on the cell id the stage 03 scripts constructed. GSE136246 is not part of the atlas: it deposits
`human_normalized_gene_counts.txt.gz` rather than raw counts, it is shallow inDrop, and seven of its
twelve donors duplicate `Zilionis_Klein_2019`.

CELL IDS, inverted exactly as `wave1_arm_*.py` built them (a wrong id would silently drop cells):
  STTT        `STTT|{SampleID}|{column.split('_',1)[1]}`, SampleID from STTT_meta.csv indexed by column
  GSE325414   `GSE325414|{sample}|{barcode}`, both from GSE325414_metadata_individual_cells.txt.gz
  GSE246215   `GSE246215|{column.split('_',1)[0]}|{column.split('_',1)[1]}`
Each source reports its join rate and stops below 100 per cent.

WHY STTT NEEDS A DIFFERENT ALGORITHM. Its deposit is a DENSE csv of 17,277 genes by 293,432 cells, so
it is read gene-major while the output must be cell-major, and the whole thing does not fit in memory.
One pass reads row chunks, writes each chunk's triplets to a scratch file and accumulates the nonzero
count per cell; a second stage allocates the output once, then scatters the scratch files into it with
a per-cell write cursor. Peak memory is the output arrays plus one chunk, and the csv is parsed once.

A WARNING THAT TRAVELS WITH THESE ARMS. `n_genes` in the wave1 arms is the PANEL-restricted count,
median 116, while in `panel_object.h5ad` it is the full-transcriptome count, median 1,170. This script
therefore recomputes `n_genes` and `n_counts` from the full matrix and keeps the panel-restricted
values as `n_genes_panel` and `n_counts_panel`, so depth matching across arms compares like with like.

RUN: python3 atlas/steps/fullbuild_01e_extract_wave1.py [sttt|gse325414|gse246215|all]
     (normally through atlas/04_build_objects.py, which passes "all")
OUT -> output/RNA/meta_dataset/fullbuild/arm_{sttt,gse325414,gse246215}_full.h5ad
"""
import os
import sys
import gzip
import glob
import numpy as np
import pandas as pd
import h5py

import os as _os  # noqa: E402
import sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # atlas/
from common import ROOT, SCRNA, SCRATCH, DATA, RSCRIPT  # noqa: E402,F401  paths from config.json
W = f"{ROOT}/output/RNA/meta_dataset/wave1"
RAW = f"{SCRNA}"
OUT = f"{ROOT}/output/RNA/meta_dataset/fullbuild"
TMP = f"{OUT}/_scratch"
TARGET = ["fibroblast", "myeloid", "T_NK"]
os.makedirs(TMP, exist_ok=True)


def dec(a):
    """Decode an array of bytes or strings to an array of str."""
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in a])


def col(h, path):
    """Read an obs column stored plain, as a categorical, or as an anndata nullable
    {mask, values} group. The wave1 arms use the third form for is_malignant and
    macrophage_strict, and a reader that only knows the first two raises on them."""
    g = h[path]
    if isinstance(g, h5py.Group):
        if "categories" in g:
            return dec(g["categories"][:])[g["codes"][:]]
        if "values" in g:
            v = g["values"][:]
            if "mask" in g:
                m = g["mask"][:].astype(bool)
                out = np.array(v, dtype=object)
                out[m] = None
                return out
            return dec(v) if v.dtype.kind in "SO" else v
        raise TypeError(f"unrecognised obs group at {path}: {list(g.keys())}")
    v = g[:]
    return dec(v) if v.dtype.kind in "SO" else v


def set_colorder(grp, cols):
    """Write the anndata `column-order` attribute (an empty list needs an explicit empty array)."""
    if cols:
        grp.attrs["column-order"] = np.array(list(cols), dtype=h5py.string_dtype())
    else:
        grp.attrs.create("column-order", np.array([], dtype="S1"))


def load_obs(arm):
    """Read the obs table of a wave-1 arm, keep the target compartments, and rename the panel-restricted n_genes and n_counts to *_panel."""
    with h5py.File(f"{W}/arm_{arm}.h5ad", "r") as h:
        ids = dec(h["obs/_index"][:])
        d = {k: col(h, f"obs/{k}") for k in h["obs"].keys() if k != "_index"}
    D = pd.DataFrame(d)
    D["cell_id"] = ids
    D = D.rename(columns={"n_genes": "n_genes_panel", "n_counts": "n_counts_panel"})
    return D[D.level1_compartment.isin(TARGET)].copy()


def write_arm(dst, sel, genes, data, indices, indptr):
    """Write one full-gene arm (CSR matrix, var, obs) as an h5ad file; n_genes and n_counts are recomputed from the full matrix."""
    n_counts = np.add.reduceat(data, indptr[:-1]) if len(data) else np.zeros(len(sel))
    n_counts = np.where(np.diff(indptr) > 0, n_counts, 0)
    sel = sel.copy()
    sel["n_genes"] = np.diff(indptr).astype("int32")
    sel["n_counts"] = n_counts.astype("int64")
    with h5py.File(dst, "w") as o:
        g = o.create_group("X")
        g.create_dataset("data", data=data.astype("float32"),
                         chunks=(min(max(len(data), 1), 1 << 20),))
        g.create_dataset("indices", data=indices.astype("int32"),
                         chunks=(min(max(len(indices), 1), 1 << 20),))
        g.create_dataset("indptr", data=indptr.astype("int64"))
        g.attrs["encoding-type"] = "csr_matrix"
        g.attrs["encoding-version"] = "0.1.0"
        g.attrs["shape"] = np.array([len(sel), len(genes)], "int64")
        v = o.create_group("var")
        v.attrs["_index"] = "gene_symbol"
        v.attrs["encoding-type"] = "dataframe"
        v.attrs["encoding-version"] = "0.2.0"
        set_colorder(v, [])
        v.create_dataset("gene_symbol", data=np.array(genes, dtype=object),
                         dtype=h5py.string_dtype())
        ob = o.create_group("obs")
        ob.attrs["_index"] = "_index"
        ob.attrs["encoding-type"] = "dataframe"
        ob.attrs["encoding-version"] = "0.2.0"
        cols = [c for c in sel.columns if c != "cell_id"]
        set_colorder(ob, cols)
        ob.create_dataset("_index", data=np.array(sel.cell_id.values, dtype=object),
                          dtype=h5py.string_dtype())
        for c in cols:
            x = sel[c].values
            if x.dtype == object or x.dtype.kind in "OU":
                ob.create_dataset(c, data=np.array(x.astype(str), dtype=object),
                                  dtype=h5py.string_dtype())
            else:
                ob.create_dataset(c, data=x)
        o.attrs["encoding-type"] = "anndata"
        o.attrs["encoding-version"] = "0.1.0"
    print(f"\nwritten -> {dst}")
    print(f"  cells {len(sel):,}   genes {len(genes):,}   nonzeros {len(data):,}   "
          f"{os.path.getsize(dst)/1e9:.2f} GB")
    print(f"  compartments: {sel.level1_compartment.value_counts().to_dict()}")
    print(f"  median full-transcriptome n_genes {int(np.median(sel.n_genes)):,} "
          f"against panel median {int(np.median(pd.to_numeric(sel.n_genes_panel))):,}")


# --------------------------------------------------------------------------- STTT
def do_sttt():
    """STTT: read the dense gene-by-cell csv in row chunks via scratch files and scatter them into a cell-major CSR matrix."""
    src = f"{RAW}/STTT_Zhang2022"
    sel = load_obs("sttt")
    meta = pd.read_csv(f"{src}/STTT_meta.csv", low_memory=False).set_index("Cells")
    print(f"STTT target cells {len(sel):,}", flush=True)
    hdr = pd.read_csv(f"{src}/raw_data.csv.gz", nrows=0)
    cols_all = [c for c in hdr.columns if c != hdr.columns[0]]
    smap = meta["SampleID"].astype(str)
    key = pd.Series([f"STTT|{smap.get(c, '?')}|{c.split('_', 1)[1]}" if "_" in c else "?"
                     for c in cols_all], index=cols_all)
    inv = pd.Series(np.arange(len(sel)), index=sel.cell_id.values)
    want = key[key.isin(inv.index)]
    rate = 100 * len(want) / len(sel)
    print(f"  columns matched to target cells: {len(want):,} of {len(sel):,} ({rate:.2f}%)", flush=True)
    if len(want) < len(sel):
        raise SystemExit("join incomplete for STTT; refusing to drop cells silently")
    outrow = inv.reindex(want.values).values.astype(np.int64)
    usecols = [hdr.columns[0]] + list(want.index)
    ncell = len(sel)
    counts = np.zeros(ncell, np.int64)
    genes, ch = [], 0
    files = []
    for chunk in pd.read_csv(f"{src}/raw_data.csv.gz", usecols=usecols, chunksize=1500,
                             index_col=0, low_memory=False):
        genes.extend(chunk.index.astype(str).tolist())
        M = chunk[list(want.index)].to_numpy(dtype=np.float32, copy=False)
        gi, ci = np.nonzero(M)
        v = M[gi, ci]
        gi = gi + (len(genes) - M.shape[0])
        cc = outrow[ci]
        np.add.at(counts, cc, 1)
        f = f"{TMP}/sttt_{ch}.npz"
        np.savez(f, c=cc.astype(np.int32), g=gi.astype(np.int32), v=v)
        files.append(f)
        ch += 1
        print(f"    chunk {ch}: genes so far {len(genes):,}, nonzeros this chunk {len(v):,}", flush=True)
    genes = np.array(genes)
    uniq, ginv = np.unique(genes, return_inverse=True)
    indptr = np.zeros(ncell + 1, np.int64)
    np.cumsum(counts, out=indptr[1:])
    total = int(indptr[-1])
    print(f"  total nonzeros {total:,}; allocating output", flush=True)
    data = np.zeros(total, np.float32)
    indices = np.zeros(total, np.int32)
    cursor = indptr[:-1].copy()
    for f in files:
        z = np.load(f)
        c_, g_, v_ = z["c"], ginv[z["g"]].astype(np.int32), z["v"]
        # A cell appears many times inside one chunk, once per detected gene. Taking cursor[c_] for
        # all of them at once would give every entry of that cell the SAME destination, so they would
        # overwrite each other and the rest of the cell's slot would stay zero.
        # So: sort by cell, then add each entry's rank within its own run to the cursor.
        order = np.argsort(c_, kind="stable")
        cs, gs, vs = c_[order], g_[order], v_[order]
        change = np.r_[True, cs[1:] != cs[:-1]]
        idxs = np.flatnonzero(change)
        first = np.repeat(idxs, np.diff(np.r_[idxs, len(cs)]))
        pos_ = cursor[cs] + (np.arange(len(cs), dtype=np.int64) - first)
        data[pos_] = vs
        indices[pos_] = gs
        np.add.at(cursor, cs, 1)
        os.remove(f)
    if not np.array_equal(cursor, indptr[1:]):
        raise SystemExit("scatter did not fill every slot exactly once; refusing to write")
    for i in range(ncell):
        a, b = int(indptr[i]), int(indptr[i + 1])
        if b > a:
            o = np.argsort(indices[a:b], kind="stable")
            indices[a:b] = indices[a:b][o]
            data[a:b] = data[a:b][o]
    write_arm(f"{OUT}/arm_sttt_full.h5ad", sel, uniq, data, indices, indptr)


# --------------------------------------------------------------------------- GSE325414
def do_325414():
    """GSE325414: select the target cells from the Matrix Market deposit and write the arm."""
    src = f"{RAW}/GSE325414"
    sel = load_obs("gse325414")
    bcs = np.array([l.strip() for l in gzip.open(f"{src}/GSE325414_barcodes.tsv.gz", "rt")])
    feats = np.array([l.strip() for l in gzip.open(f"{src}/GSE325414_features.tsv.gz", "rt")])
    meta = pd.read_csv(f"{src}/GSE325414_metadata_individual_cells.txt.gz", sep="\t",
                       low_memory=False).set_index("cell")
    samp = meta["sample"].astype(str)
    key = np.array([f"GSE325414|{samp.get(b, '?')}|{b}" for b in bcs])
    inv = pd.Series(np.arange(len(sel)), index=sel.cell_id.values)
    outrow = inv.reindex(key).values
    got = int(np.isfinite(outrow.astype(float)).sum())
    print(f"GSE325414 target cells {len(sel):,}, matched barcodes {got:,}", flush=True)
    if got < len(sel):
        raise SystemExit("join incomplete for GSE325414")
    outrow = np.where(pd.isna(outrow), -1, outrow).astype(np.int64)
    uniq, ginv = np.unique(feats, return_inverse=True)
    R, C, V = [], [], []
    with gzip.open(f"{src}/GSE325414_matrix.mtx.gz", "rt") as f:
        f.readline()
        while True:
            line = f.readline()
            if not line.startswith("%"):
                break
        for chunk in pd.read_csv(f, sep=r"\s+", header=None, chunksize=10_000_000, engine="c"):
            t = chunk.values
            g_ = t[:, 0].astype(np.int64) - 1
            c_ = t[:, 1].astype(np.int64) - 1
            v_ = t[:, 2].astype("float32")
            keep = outrow[c_] >= 0
            R.append(ginv[g_[keep]].astype(np.int32))
            C.append(outrow[c_[keep]].astype(np.int32))
            V.append(v_[keep])
            print(f"    kept {sum(len(x) for x in V):,} nonzeros", flush=True)
    build_from_coo(f"{OUT}/arm_gse325414_full.h5ad", sel, uniq,
                   np.concatenate(C), np.concatenate(R), np.concatenate(V))


# --------------------------------------------------------------------------- GSE246215
def do_246215():
    """GSE246215: select the target cells from the dense fibroblast count csv and write the arm."""
    src = f"{RAW}/GSE246215"
    sel = load_obs("gse246215")
    hdr = pd.read_csv(f"{src}/GSE246215_Fibroblast_counts.csv.gz", nrows=0)
    cols_all = [c for c in hdr.columns if c != hdr.columns[0]]
    key = pd.Series([f"GSE246215|{c.split('_', 1)[0]}|{c.split('_', 1)[1]}" if "_" in c else "?"
                     for c in cols_all], index=cols_all)
    inv = pd.Series(np.arange(len(sel)), index=sel.cell_id.values)
    want = key[key.isin(inv.index)]
    print(f"GSE246215 target cells {len(sel):,}, matched columns {len(want):,}", flush=True)
    if len(want) < len(sel):
        raise SystemExit("join incomplete for GSE246215")
    outrow = inv.reindex(want.values).values.astype(np.int64)
    M = pd.read_csv(f"{src}/GSE246215_Fibroblast_counts.csv.gz",
                    usecols=[hdr.columns[0]] + list(want.index), index_col=0, low_memory=False)
    genes = M.index.astype(str).values
    uniq, ginv = np.unique(genes, return_inverse=True)
    A = M[list(want.index)].to_numpy(dtype=np.float32)
    gi, ci = np.nonzero(A)
    build_from_coo(f"{OUT}/arm_gse246215_full.h5ad", sel, uniq,
                   outrow[ci].astype(np.int32), ginv[gi].astype(np.int32), A[gi, ci])


def build_from_coo(dst, sel, genes, cells, gidx, vals):
    """Sort coordinate-format entries by cell and gene, build the CSR index pointer and write the arm."""
    n = len(sel)
    order = np.lexsort((gidx, cells))
    cells, gidx, vals = cells[order], gidx[order], vals[order]
    counts = np.bincount(cells, minlength=n)
    indptr = np.zeros(n + 1, np.int64)
    np.cumsum(counts, out=indptr[1:])
    write_arm(dst, sel, genes, vals, gidx, indptr)


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("all", "gse246215"):
        do_246215()
    if which in ("all", "gse325414"):
        do_325414()
    if which in ("all", "sttt"):
        do_sttt()
