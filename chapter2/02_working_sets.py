#!/usr/bin/env python3
"""
Chapter 2, step 02: the fibroblast, macrophage, T/NK and cancer-cell working sets, extracted from the atlas object.

WHAT: for each cell set, the cells of

          tissue_class == "primary_tumour" AND block_qc_pass AND histology in {LUSC, LUAD}   ("primary_tumour" is the atlas value)
          AND level1_compartment in <set> [AND level2_celltype == <set>]
      AND not one of the redistributed copies of the Lambrechts 2018 patients (redox.atlas_fixes)

      over all genes of the object, with the metadata the tests need (batch, patient, arm, genes
      detected, level 1 and level 2 labels). The presets are:
          fib   level1_compartment fibroblast
          mac   level1_compartment myeloid, level2_celltype macrophage
          tnk   level1_compartment T_NK
          can   cancer cells (Supplementary Figure 2c): built only when named with --cells can, from the
                epithelial object written by atlas/06_cancer_cell_object.py (configuration key
                "epi_object"), keeping the cells the original authors labeled as cancer cells
      Every other script of chapter 2 reads these working sets (through ch2lib.load_set).
ITEMS: none directly. The working sets feed Figure 2 (a, c), Figure 3, Supplementary Figures 2 to 4 and
      Supplementary Data 1.
FROM: the atlas object named by the configuration key "atlas_object" (fullgene_object.h5ad, 1,888,910
      cells by 50,552 genes, counts in CSR form), read with h5py alone.
TO:   work_root/intermediate/chapter2/work_<cell set>/
          genes.npy        gene symbols, in the column order of the submatrices
          present.pkl      varm/gene_present, genes x 7 arms (is the gene measured in that arm)
          raw_LUSC.npz     sub (row in the object), data, indices, indptr, n, ngenes of the LUSC rows
          raw_LUAD.npz     the same for LUAD
          meta_LUSC.npz    block, patient, arm, ngenes_obs, l1, l2 of the LUSC rows
          meta_LUAD.npz    the same for LUAD
          extract_log.txt  the console log
          checksums.csv    one row per array: file, array, dtype, shape, SHA-256 of its bytes (+ .md)
WHY:  the later steps work on a small fixed subset of the 1.9 million cells. The selection leaves out the
      50,651 cells of the eight redistributed Lambrechts 2018 patients and writes the corrected study
      label (redox.atlas_fixes).

      The object is read without loading X: X/indptr is read whole, then X/data and X/indices are
      streamed once in windows of --window elements, copying out only the rows of the working set,
      so the peak memory is a few hundred MB to 2 GB rather than the 22 GB of the whole matrix.

RUN:  python3 chapter2/02_working_sets.py [--cells fib mac tnk can] [--window 16000000]
      Default: fib mac tnk. About 15 to 25 seconds per cell set.
"""
import argparse
import hashlib
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ch2lib as L  # noqa: E402
from redox import atlas_fixes, paths  # noqa: E402

ARM_MAP = {"LuCA_census": "LuCA", "EMTAB13526": "E-MTAB-13526", "EMTAB6653": "E-MTAB-6653",
           "EMTAB6149": "E-MTAB-6149", "STTT": "STTT", "GSE325414": "GSE325414",
           "GSE246215": "GSE246215"}


def atlas_object():
    """The atlas object from the configuration (key "atlas_object")."""
    named = paths.CFG.get("atlas_object")
    if not named:
        raise SystemExit(f'{paths.CFG["_config_path"]} has no "atlas_object" entry: add the path of '
                         "fullgene_object.h5ad (see README, Chapter 2)")
    p = Path(named).expanduser()
    if not p.exists():
        raise SystemExit(f"atlas object not found: {p}")
    return p


def strcol(ds):
    """an obs column of variable-length strings, as a numpy array of str"""
    v = ds[:]
    if v.dtype == object:
        return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in v])
    return v


def strcol_at(ds, idx):
    """the same at the given rows; the column is read whole and subset in memory, which is three
    orders of magnitude faster than h5py fancy indexing on a large index array"""
    v = ds[:][idx]
    if v.dtype == object:
        return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in v])
    return v


def object_for(cell):
    """The object a cell set is read from: the atlas object, or for a set of ch2lib.EXTRA_CELL_SETS the object its
    "object" key names in the configuration ("epi_object", the epithelial object built by
    atlas/06_cancer_cell_object.py)."""
    key = L.cell_spec(cell).get("object")
    if not key:
        return atlas_object()
    named = paths.CFG.get(key)
    if not named or not Path(named).expanduser().exists():
        raise SystemExit(f'{paths.CFG["_config_path"]} has no usable "{key}" entry: add the path of the epithelial '
                         "object written by atlas/06_cancer_cell_object.py")
    return Path(named).expanduser()


def extract(cell, h5ad, window):
    """Write the working-set files of one cell set (see the header) by streaming the object once; returns the folder."""
    spec = L.cell_spec(cell)
    W = L.work_dir(cell)
    log = L.Log()
    log(f"02_working_sets.py, cell set {cell} ({spec['label']})")
    log(f"  object    {h5ad.name}")
    log(f"  window    {window:,} elements ({window * 8 / 1e6:.0f} MB of buffers)")
    with h5py.File(h5ad, "r") as f:
        shape = tuple(int(x) for x in f["X"].attrs["shape"])
        nnz_total = int(f["X/data"].shape[0])
        log(f"  X         {shape[0]:,} cells x {shape[1]:,} genes, {nnz_total:,} stored nonzeros")
        gk = f["var"].attrs["_index"]
        gk = gk.decode() if isinstance(gk, bytes) else gk
        genes = strcol(f[f"var/{gk}"])
        gp = f["varm/gene_present"]
        arms = [a.decode() if isinstance(a, bytes) else a for a in gp.attrs["column-order"]]
        present = pd.DataFrame({a: gp[a][:] for a in arms}, index=genes)
        np.save(W / "genes.npy", genes)
        present.to_pickle(W / "present.pkl")
        if L.SPLITTER not in set(genes):
            raise SystemExit("TXNRD1 is not on the gene axis of this object")
        log(f"  genes     {len(genes):,} symbols, gene_present has {len(arms)} arms: {arms}")

        t0 = time.time()
        keep = (strcol(f["obs/tissue_class"]) == "primary_tumour")  # atlas value
        keep &= f["obs/block_qc_pass"][:].astype(bool)
        l1_all = strcol(f["obs/level1_compartment"])
        keep &= np.isin(l1_all, spec["l1"])
        if spec.get("l3_native"):    # the original authors' own labels (cancer cells)
            keep &= np.isin(strcol(f["obs/level3_celltype_native"]), spec["l3_native"])
        if spec["l2"]:
            keep &= (strcol(f["obs/level2_celltype"]) == spec["l2"])
        dup = atlas_fixes.duplicate(strcol(f["obs/patient_id"]))
        log(f"  left out  {int((keep & dup).sum()):,} cells of the redistributed Lambrechts 2018 patients "
            f"(redox.atlas_fixes)")
        keep &= ~dup
        hist_all = strcol(f["obs/histology"])
        sel = {h: np.where(keep & (hist_all == h))[0] for h in L.HISTS}
        del keep, hist_all
        for h in L.HISTS:
            log(f"  {h}      {len(sel[h]):,} working-set cells")
        log(f"  obs read and filtered in {time.time() - t0:.1f}s")

        l2_all = strcol(f["obs/level2_celltype"])
        for h in L.HISTS:
            idx = sel[h]
            src = strcol_at(f["obs/source_key"], idx)
            miss = sorted(set(src) - set(ARM_MAP))
            if miss:
                raise SystemExit(f"source_key values with no gene_present column: {miss}")
            np.savez(W / f"meta_{h}.npz",
                     block=atlas_fixes.relabel(strcol_at(f["obs/block_id"], idx)),
                     patient=atlas_fixes.relabel(strcol_at(f["obs/patient_id"], idx)),
                     arm=np.array([ARM_MAP[s] for s in src]),
                     ngenes_obs=np.asarray(f["obs/n_genes"][:][idx], float),
                     l1=l1_all[idx], l2=l2_all[idx])
            log(f"  wrote meta_{h}.npz")
        del l1_all, l2_all

        indptr = f["X/indptr"][:]
        res = {}
        for h in L.HISTS:
            idx = sel[h]
            lens = (indptr[idx + 1] - indptr[idx]).astype(np.int64)
            ptr = np.zeros(len(idx) + 1, dtype=np.int64)
            ptr[1:] = np.cumsum(lens)
            res[h] = dict(idx=idx, ptr=ptr, data=np.empty(int(ptr[-1]), np.float32),
                          indices=np.empty(int(ptr[-1]), np.int32))
            log(f"  {h} submatrix will hold {int(ptr[-1]):,} nonzeros")
        # one flat list of (source start, source end, histology, destination offset), sorted by source
        # start, so the windows are walked once
        jobs = []
        for h in L.HISTS:
            for k, r in enumerate(res[h]["idx"]):
                jobs.append((int(indptr[r]), int(indptr[r + 1]), h, int(res[h]["ptr"][k])))
        jobs.sort()
        Xd, Xi = f["X/data"], f["X/indices"]
        t0 = time.time()
        j0 = 0
        w0 = jobs[0][0] - (jobs[0][0] % window)
        nwin = 0
        while j0 < len(jobs) and w0 < nnz_total:
            w1 = min(w0 + window, nnz_total)
            j = j0
            touched = False
            while j < len(jobs) and jobs[j][0] < w1:
                s, e, h, off = jobs[j]
                lo, hi = max(s, w0), min(e, w1)
                if hi > lo:
                    if not touched:
                        bd = Xd[w0:w1]
                        bi = Xi[w0:w1]
                        touched = True
                        nwin += 1
                    res[h]["data"][off + lo - s: off + hi - s] = bd[lo - w0: hi - w0]
                    res[h]["indices"][off + lo - s: off + hi - s] = bi[lo - w0: hi - w0]
                j += 1
            if touched:
                del bd, bi
            while j0 < len(jobs) and jobs[j0][1] <= w1:
                j0 += 1
            w0 = w1
        log(f"  streamed {nwin:,} windows in {time.time() - t0:.1f}s (peak resident {L.rss_gb():.2f} GB)")

    for h in L.HISTS:
        o = res[h]
        np.savez(W / f"raw_{h}.npz", sub=o["idx"], data=o["data"], indices=o["indices"], indptr=o["ptr"],
                 n=np.array([len(o["idx"])]), ngenes=np.array([len(genes)]))
        log(f"  wrote raw_{h}.npz, {len(o['idx']):,} cells, {int(o['ptr'][-1]):,} nonzeros")
    log(f"DONE {cell} in {log.elapsed()}")
    log.write(W / "extract_log.txt")
    write_checksums(W)
    return W


def arrays_of(path):
    """(array name, numpy array) of one working-set file."""
    if path.suffix == ".npz":
        z = np.load(path, allow_pickle=True)
        return [(k, z[k]) for k in sorted(z.files)]
    if path.suffix == ".npy":
        return [("genes", np.load(path, allow_pickle=True))]
    df = pd.read_pickle(path)
    return [("index", np.asarray(df.index))] + [(c, df[c].to_numpy()) for c in df.columns]


def digest(a):
    """SHA-256 of an array's bytes (of its text for object arrays)."""
    a = np.ascontiguousarray(a)
    if a.dtype == object:
        return hashlib.sha256("\x1f".join(map(str, a.ravel())).encode()).hexdigest()
    return hashlib.sha256(a.tobytes()).hexdigest()


FILES = ("genes.npy", "present.pkl", "meta_LUSC.npz", "meta_LUAD.npz", "raw_LUSC.npz", "raw_LUAD.npz")


def write_checksums(W):
    """checksums.csv of a working set: one row per array of its six files."""
    rows = []
    for fn in FILES:
        for name, a in arrays_of(W / fn):
            rows.append(dict(file=fn, array=name, dtype=str(a.dtype), shape="x".join(map(str, a.shape)),
                             sha256=digest(a)))
    out = W / "checksums.csv"
    pd.DataFrame(rows).to_csv(out, index=False)
    L.sidecar(out, "SHA-256 of every array of the working set, for comparing two extractions.",
              [("file", "file of the working set"), ("array", "array name inside the file"),
               ("dtype", "numpy dtype"), ("shape", "array shape"), ("sha256", "SHA-256 of the array bytes")])


def main():
    """Extract the working sets of the cell sets named with --cells."""
    ap = argparse.ArgumentParser(description=__doc__.split("RUN:")[0].strip(),
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cells", nargs="*", default=list(L.CELL_SETS),
                    choices=list(L.CELL_SETS) + list(L.EXTRA_CELL_SETS))
    ap.add_argument("--window", type=int, default=16_000_000)
    a = ap.parse_args()
    for cell in a.cells:
        extract(cell, object_for(cell), a.window)


if __name__ == "__main__":
    main()
