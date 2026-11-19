#!/usr/bin/env python3
"""
Atlas stage 04, step merge_panel: merge the three wave-1 arms into the panel object, so block QC can
run on them.

Paper items: writes panel_object_v2.h5ad (2,500,510 cells x 1,721 genes, 24 blocks), on which the block
QC of the atlas behind Figure 2, Figure 3, Supplementary Figures 2 to 4 and Supplementary Data 1 is
computed. No random numbers are drawn.

The three wave-1 arms carry `block_qc_pass = False` until they have been assessed, which would exclude
all 301,843 of their cells from every analysis that applies the block QC. This merge puts them on the
panel gene axis so the block QC can assess them.

WHY THE QC CANNOT RUN ON THE FULL-TRANSCRIPTOME OBJECT INSTEAD. `blockqc_02` asks, in effect, which
genes are detected everywhere else and zero here. On the 50,552-gene union axis LuCA has 32,741
structural zeros and STTT 33,275, purely from measuring fewer genes, so every block would fail for the
wrong reason. The QC has to run where it was calibrated: the panel axis, where every arm has near
complete coverage.

THE GENE AXIS SHRINKS. Intersecting the 1,906 panel genes with the three arms gives **1,721**, losing
185. The redox, CAF, complement and lineage marker genes scored in the analyses all survive (TXNRD1,
the NFE2L2 target genes, AKR1C1/AKR1C2/ALDH1A1, the collagens, C3/CFD/C1S/C7, CD3E, NKG7, NCAM1). The
panel-axis object is used for the block QC; the analyses read the full-gene object.

GSE136246 IS NOT INCLUDED: normalized counts rather than raw, shallow inDrop, and
seven of its twelve donors duplicate `Zilionis_Klein_2019`. The merged object therefore carries the
same seven arms as `fullgene_object.h5ad`, so the two stay in step.

TWO THINGS FIXED IN PASSING
  `block_id` and `level1_compartment` are written as CATEGORICALS, because `blockqc_01` and
  `blockqc_02` read `/categories` and `/codes` directly and would fail on plain strings.
  `n_genes` IS MADE COMPARABLE. In `panel_object` it is the full-transcriptome count, median 1,170;
  in the wave1 arms it is the panel-restricted count, median 116. Depth matching uses it with a gap of
  150 and would pair almost anything across that boundary. Here `n_genes` is the full-transcriptome
  count for every cell, taken from `panel_object` for the four original arms and from the matching
  `fullbuild/arm_*_full.h5ad` for the three wave-1 ones, and `n_genes_panel` is the count on this
  object's own 1,721-gene axis.

RUN: python3 atlas/steps/fullbuild_05_merge_panel.py   (normally through atlas/04_build_objects.py)
FROM: stage1a/panel_object.h5ad, wave1/arm_*.h5ad and fullbuild/arm_*_full.h5ad (wave-1 arms).
OUT -> output/RNA/meta_dataset/fullbuild/panel_object_v2.h5ad
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
W = f"{ROOT}/output/RNA/meta_dataset/wave1"
FB = f"{ROOT}/output/RNA/meta_dataset/fullbuild"
DST = f"{FB}/panel_object_v2.h5ad"
NEW = [("STTT", f"{W}/arm_sttt.h5ad", f"{FB}/arm_sttt_full.h5ad"),
       ("GSE325414", f"{W}/arm_gse325414.h5ad", f"{FB}/arm_gse325414_full.h5ad"),
       ("GSE246215", f"{W}/arm_gse246215.h5ad", f"{FB}/arm_gse246215_full.h5ad")]
CHUNK = 50000
CAT_AS_CATEGORICAL = ["block_id", "level1_compartment", "level2_celltype", "histology",
                      "tissue_class", "source_key", "study", "platform", "donor_id", "patient_id"]


def dec(a):
    """Decode an array of bytes or strings to an array of str."""
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in a])


def col(h, path):
    """Read one obs column of an h5ad file (categorical, nullable or plain) as an array."""
    g = h[path]
    if isinstance(g, h5py.Group):
        if "categories" in g:
            return dec(g["categories"][:])[g["codes"][:]]
        if "values" in g:
            v = g["values"][:]
            if "mask" in g:
                out = np.array(v, dtype=object)
                out[g["mask"][:].astype(bool)] = None
                return out
            return dec(v) if v.dtype.kind in "SO" else v
        raise TypeError(path)
    v = g[:]
    return dec(v) if v.dtype.kind in "SO" else v


def gene_axis(f):
    """Return the gene symbols of an h5ad file (its var index)."""
    with h5py.File(f, "r") as h:
        gk = h["var"].attrs.get("_index", "gene_symbol")
        gk = gk.decode() if isinstance(gk, bytes) else str(gk)
        return dec(h[f"var/{gk}"][:])


print("computing the intersection gene axis ...", flush=True)
pg = gene_axis(PANEL)
inter = set(pg)
for _, p, _ in NEW:
    inter &= set(gene_axis(p))
axis = np.array([g for g in pg if g in inter])          # keep panel order
print(f"  panel {len(pg)} -> merged {len(axis)} (loses {len(pg)-len(axis)})", flush=True)

sources = [("panel", PANEL, None)] + NEW
obscols, have, kind = None, {}, {}
for _, p, _ in sources:
    with h5py.File(p, "r") as h:
        c = set(h["obs"].keys()) - {"_index"}
        for k in c:
            have.setdefault(k, set()).add(p)
            if k not in kind:
                g = h[f"obs/{k}"]
                if isinstance(g, h5py.Group):
                    kind[k] = "str"
                else:
                    kind[k] = "num" if g.dtype.kind in "fiub" else "str"
    obscols = c if obscols is None else (obscols | c)
obscols = sorted(obscols)
# A numeric column missing from some source must NOT be filled with a string: that would make the
# concatenation object dtype and write numeric columns as categoricals. Numeric columns with a gap
# are promoted to float64 and filled with NaN; only genuine string columns are filled with "".
partial = {k for k in obscols if len(have[k]) < len(sources)}
print(f"  obs columns, union across sources: {len(obscols)}", flush=True)
print(f"  present in every source: {len(obscols)-len(partial)}; "
      f"partial: {len(partial)} ({sum(1 for k in partial if kind[k]=='num')} numeric, "
      f"filled with NaN as float64)", flush=True)

ids_acc, obs_acc, ptr, nnz = [], {c: [] for c in obscols}, [0], 0
ngenes_full_acc = []
with h5py.File(DST, "w") as o:
    g = o.create_group("X")
    d_out = g.create_dataset("data", (0,), maxshape=(None,), dtype="float32", chunks=(1 << 20,))
    i_out = g.create_dataset("indices", (0,), maxshape=(None,), dtype="int32", chunks=(1 << 20,))
    for name, p, fullp in sources:
        src_axis = gene_axis(p)
        pos = pd.Series(np.arange(len(axis)), index=axis)
        remap = pos.reindex(src_axis).values                    # NaN where the gene is dropped
        keepg = np.isfinite(remap.astype(float))
        remap_i = np.where(keepg, np.nan_to_num(remap), -1).astype(np.int64)
        with h5py.File(p, "r") as h:
            ip = h["X/indptr"][:]
            nc = len(ip) - 1
            cid = dec(h["obs/_index"][:])
            ids_acc.append(cid)
            for c in obscols:
                if c in h["obs"]:
                    v = col(h, f"obs/{c}")
                    obs_acc[c].append(np.asarray(v, float) if (kind[c] == "num" and c in partial)
                                      else v)
                else:
                    obs_acc[c].append(np.full(nc, np.nan) if kind[c] == "num"
                                      else np.array([""] * nc, dtype=object))
            # full-transcriptome n_genes: from panel_object directly, else from the full arm
            if fullp is None:
                ngenes_full_acc.append(np.asarray(col(h, "obs/n_genes"), float))
            else:
                with h5py.File(fullp, "r") as hf:
                    fid = dec(hf["obs/_index"][:])
                    fng = np.asarray(hf["obs/n_genes"][:], float)
                s = pd.Series(fng, index=fid)
                ngenes_full_acc.append(s.reindex(cid).values)
            for a in range(0, nc, CHUNK):
                b = min(a + CHUNK, nc)
                lo, hi = int(ip[a]), int(ip[b])
                dat = h["X/data"][lo:hi]
                ind = remap_i[h["X/indices"][lo:hi]]
                loc = ip[a:b + 1] - lo
                cd, ci = [], []
                for r in range(b - a):
                    s0, s1 = int(loc[r]), int(loc[r + 1])
                    gg, vv = ind[s0:s1], dat[s0:s1]
                    m = gg >= 0
                    gg, vv = gg[m], vv[m]
                    srt = np.argsort(gg, kind="stable")
                    cd.append(vv[srt])
                    ci.append(gg[srt])
                    ptr.append(ptr[-1] + len(gg))
                cd = np.concatenate(cd) if cd else np.array([], "float32")
                ci = np.concatenate(ci) if ci else np.array([], "int64")
                d_out.resize((nnz + len(cd),))
                i_out.resize((nnz + len(ci),))
                d_out[nnz:] = cd.astype("float32")
                i_out[nnz:] = ci.astype("int32")
                nnz += len(cd)
            print(f"  merged {name:11s} {nc:>9,} cells, running nonzeros {nnz:,}", flush=True)
    total = len(ptr) - 1
    g.create_dataset("indptr", data=np.asarray(ptr, "int64"))
    g.attrs["encoding-type"] = "csr_matrix"
    g.attrs["encoding-version"] = "0.1.0"
    g.attrs["shape"] = np.array([total, len(axis)], "int64")

    v = o.create_group("var")
    v.attrs.update({"_index": "gene_symbol", "encoding-type": "dataframe",
                    "encoding-version": "0.2.0"})
    v.attrs.create("column-order", np.array([], dtype="S1"))
    v.create_dataset("gene_symbol", data=np.array(axis, dtype=object), dtype=h5py.string_dtype())

    ob = o.create_group("obs")
    ob.attrs.update({"_index": "_index", "encoding-type": "dataframe",
                     "encoding-version": "0.2.0"})
    ob.create_dataset("_index", data=np.array(np.concatenate(ids_acc), dtype=object),
                      dtype=h5py.string_dtype())
    written = []
    for c in obscols:
        x = np.concatenate([np.asarray(a, dtype=object) if np.asarray(a).dtype.kind in "OU"
                            else np.asarray(a) for a in obs_acc[c]])
        if c == "n_genes":
            continue
        if kind[c] == "num":
            ob.create_dataset(c, data=np.asarray(x, float) if c in partial else x)
            written.append(c)
            continue
        if c in CAT_AS_CATEGORICAL or x.dtype.kind in "OU":
            xs = pd.Series(x).astype(str).fillna("")
            cats, codes = np.unique(xs.values), None
            codes = pd.Categorical(xs, categories=cats).codes.astype("int32")
            grp = ob.create_group(c)
            grp.attrs.update({"encoding-type": "categorical", "encoding-version": "0.2.0",
                              "ordered": False})
            grp.create_dataset("categories", data=np.array(cats, dtype=object),
                               dtype=h5py.string_dtype())
            grp.create_dataset("codes", data=codes)
        else:
            ob.create_dataset(c, data=x)
        written.append(c)
    ng_full = np.concatenate([np.asarray(a, float) for a in ngenes_full_acc])
    ob.create_dataset("n_genes", data=np.nan_to_num(ng_full).astype("int32"))
    ob.create_dataset("n_genes_panel", data=np.diff(np.asarray(ptr, "int64")).astype("int32"))
    written += ["n_genes", "n_genes_panel"]
    ob.attrs["column-order"] = np.array(sorted(written), dtype=h5py.string_dtype())
    # stage2_merge_02_blockqc.py writes its record into uns, so the group must exist
    u = o.create_group("uns")
    u.attrs.update({"encoding-type": "dict", "encoding-version": "0.1.0"})
    o.attrs.update({"encoding-type": "anndata", "encoding-version": "0.1.0"})

print(f"\nwritten -> {DST}")
print(f"  cells {total:,}   genes {len(axis):,}   nonzeros {nnz:,}   "
      f"{os.path.getsize(DST)/1e9:.2f} GB")
with h5py.File(DST, "r") as h:
    b = dec(h["obs/block_id/categories"][:])
    print(f"  blocks: {len(b)}")
    ngf = h["obs/n_genes"][:]
    ngp = h["obs/n_genes_panel"][:]
    print(f"  n_genes (full transcriptome) median {int(np.median(ngf)):,}, "
          f"n_genes_panel median {int(np.median(ngp)):,}")
    print(f"  cells with no full-transcriptome n_genes: {int((ngf == 0).sum()):,}")
