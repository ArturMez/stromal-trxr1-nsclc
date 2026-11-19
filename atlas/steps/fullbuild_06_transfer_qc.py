#!/usr/bin/env python3
"""
Atlas stage 04, step transfer_qc: transfer block_qc_pass onto the full-transcriptome object, and
repair n_genes on panel_object_v2.

Paper items: sets which cells of fullgene_object.h5ad (the atlas behind Figure 2, Figure 3,
Supplementary Figures 2 to 4 and Supplementary Data 1) pass the block QC (23 of 24 blocks pass).
No random numbers are drawn.

The block QC has to run on the panel axis, where it was calibrated, so it ran on
`panel_object_v2.h5ad`. Its verdict is a property of the BLOCK, not of the gene axis, so it transfers
onto `fullgene_object.h5ad` by cell id without being recomputed.

RESULT BEING TRANSFERRED, from `stage2_merge_02_blockqc.py` on 24 blocks, same rule as the original 21
(a gene counts as a dropout when its count is 0 across a block with expected detection above 5 %; a
block fails at FAIL_N_DROPOUT = 10 dropout genes):
  Adams_Kaminski_2020 10x 3' v2   18 dropouts   FAILS
  every other block                0 or 1       passes
  the three wave-1 blocks          0, 0 and 1   PASS
Adams_Kaminski is the only failure at every threshold from 2 to 15, so the verdict does not depend on
where the line sits.

TWO REPAIRS MADE HERE
  `block_qc_assessed` becomes True for every cell, since all 24 blocks have now been assessed. Before
  this it was False for the three wave-1 arms, which meant "not assessed" and excluded them.
  `n_genes` in `panel_object_v2` is 0 for 155,300 cells. Those are cells of the new arms in
  compartments the full-gene arms do not carry, because the full-gene build covers only fibroblast,
  myeloid and T/NK while the panel merge carries every compartment for the QC's composition
  adjustment. A zero there would mislead depth matching, which would pair zeros with each other. They become -1 and a boolean `n_genes_is_full` records which cells have a real value.

RUN: python3 atlas/steps/fullbuild_06_transfer_qc.py   (normally through atlas/04_build_objects.py)
FROM: blockqc_v2 verdict stored in panel_object_v2.h5ad; modifies fullgene_object.h5ad (obs) and
      panel_object_v2.h5ad (n_genes, n_genes_is_full) in place.
"""
import numpy as np
import pandas as pd
import h5py

import os as _os  # noqa: E402
import sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # atlas/
from common import ROOT, SCRNA, SCRATCH, DATA, RSCRIPT  # noqa: E402,F401  paths from config.json
FB = f"{ROOT}/output/RNA/meta_dataset/fullbuild"
PANEL2 = f"{FB}/panel_object_v2.h5ad"
FULL = f"{FB}/fullgene_object.h5ad"


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
            return g["values"][:]
    v = g[:]
    return dec(v) if v.dtype.kind in "SO" else v


print("reading the QC verdict from panel_object_v2 ...", flush=True)
with h5py.File(PANEL2, "r") as h:
    pid = dec(h["obs/_index"][:])
    q = h["obs/block_qc_pass"][:].astype(bool)
    blk = col(h, "obs/block_id").astype(str)
verdict = pd.Series(q, index=pid)
byblock = pd.DataFrame({"block": blk, "pass": q}).groupby("block")["pass"].agg(["size", "mean"])
print(f"  {len(verdict):,} cells, {int(q.sum()):,} pass, {int((~q).sum()):,} fail", flush=True)
print(f"  failing blocks: {sorted(byblock[byblock['mean'] < 1].index.tolist())}", flush=True)

print("\ntransferring onto fullgene_object ...", flush=True)
with h5py.File(FULL, "r+") as h:
    fid = dec(h["obs/_index"][:])
    old = col(h, "obs/block_qc_pass")
    old_b = pd.Series(old).astype(str).isin(["True", "1", "true"]).values
    new = verdict.reindex(fid)
    miss = int(new.isna().sum())
    if miss:
        raise SystemExit(f"{miss:,} full-gene cells have no verdict; refusing to guess")
    new_b = new.values.astype(bool)
    for k in ("block_qc_pass", "block_qc_assessed"):
        if k in h["obs"]:
            del h["obs"][k]
    h["obs"].create_dataset("block_qc_pass", data=new_b, dtype=bool)
    h["obs"].create_dataset("block_qc_assessed", data=np.ones(len(fid), bool), dtype=bool)
    sk = col(h, "obs/source_key").astype(str)
    co = list(dec(h["obs"].attrs["column-order"])) if "column-order" in h["obs"].attrs else []
    for k in ("block_qc_pass", "block_qc_assessed"):
        if k not in co:
            co.append(k)
    h["obs"].attrs["column-order"] = np.array(co, dtype=h5py.string_dtype())
D = pd.DataFrame({"source": sk, "was": old_b, "now": new_b})
print(f"  cells changing from excluded to usable: {int((~D.was & D.now).sum()):,}", flush=True)
print(f"  cells changing from usable to excluded: {int((D.was & ~D.now).sum()):,}", flush=True)
print("\n" + "=" * 92)
print("block_qc_pass ON THE FULL-TRANSCRIPTOME OBJECT, AFTER THE TRANSFER")
print("VARIABLE: does the cell's BLOCK pass the composition- and depth-adjusted dropout test.")
print("UNIT: cell.  DENOMINATOR: all cells of that source in the object.")
print("=" * 92)
print(D.groupby("source").agg(cells=("now", "size"),
                              pct_pass_before=("was", lambda x: round(100 * x.mean(), 1)),
                              pct_pass_now=("now", lambda x: round(100 * x.mean(), 1))).to_string())
print(f"\nTOTAL usable: {int(D.now.sum()):,} of {len(D):,} cells "
      f"(was {int(D.was.sum()):,})")

print("\nrepairing n_genes in panel_object_v2 ...", flush=True)
with h5py.File(PANEL2, "r+") as h:
    ng = h["obs/n_genes"][:]
    isfull = ng > 0
    ng2 = np.where(isfull, ng, -1).astype("int32")
    del h["obs/n_genes"]
    h["obs"].create_dataset("n_genes", data=ng2)
    if "n_genes_is_full" in h["obs"]:
        del h["obs"]["n_genes_is_full"]
    h["obs"].create_dataset("n_genes_is_full", data=isfull, dtype=bool)
    co = list(dec(h["obs"].attrs["column-order"])) if "column-order" in h["obs"].attrs else []
    if "n_genes_is_full" not in co:
        co.append("n_genes_is_full")
    h["obs"].attrs["column-order"] = np.array(sorted(set(co)), dtype=h5py.string_dtype())
print(f"  {int((~isfull).sum()):,} cells set to -1 with n_genes_is_full False; "
      f"{int(isfull.sum()):,} carry a real full-transcriptome count")
