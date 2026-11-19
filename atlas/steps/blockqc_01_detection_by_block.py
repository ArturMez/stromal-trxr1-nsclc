"""
Block QC step 1: per-gene detection rate for every block, independent recomputation.

Reads output/RNA/meta_dataset/stage1a/panel_object.h5ad directly with h5py (CSR, chunked)
and writes:
  blockqc_detection_by_block.csv   genes x blocks detection rate
  blockqc_block_summary.csv        per-block dropout / inflation counts
  blockqc_cell_flags.npz           per-cell detection booleans for probe genes
No cells are modified here. Read only.
Block QC rule (stage 01): a gene counts as a dropout in a block when its count is 0 across the block
and the rest of the cohort detects it in more than 5 percent of cells; the rule itself is applied in
blockqc_02_rule_and_flag.py. Run by atlas/01_arms_luca_lambrechts.py (step blockqc_01), no options.
No random numbers are drawn.
"""
import json
import os
import numpy as np
import pandas as pd
import h5py

import os as _os  # noqa: E402
import sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # atlas/
from common import MD  # noqa: E402  paths from config.json, see atlas/common.py
H5 = MD + "/stage1a/panel_object.h5ad"
from common import MD  # noqa: E402  paths from config.json, see atlas/common.py
OUT = MD + "/stage1a/"
PROBES = ["TXNRD1", "PRDX2", "CD68", "IGKC", "NFE2L2", "TXN", "PRDX1", "ACTA2"]
CHUNK = 100_000

f = h5py.File(H5, "r")

genes = np.array([g.decode() if isinstance(g, bytes) else g for g in f["var/gene_symbol"][:]])
n_genes = genes.size
blocks = np.array([b.decode() for b in f["obs/block_id/categories"][:]])
bcodes = f["obs/block_id/codes"][:].astype(np.int64)
n_blocks = blocks.size
n_cells = bcodes.size
assert n_genes == 1906, n_genes
assert n_cells == 1_146_594, n_cells

indptr = f["X/indptr"][:]
X_data = f["X/data"]
X_ind = f["X/indices"]

det = np.zeros(n_blocks * n_genes, dtype=np.int64)     # detection counts, block-major
probe_idx = np.array([int(np.where(genes == p)[0][0]) for p in PROBES])
probe_hit = np.zeros((n_cells, len(PROBES)), dtype=bool)
neg_entries = 0
explicit_zero_entries = 0

for start in range(0, n_cells, CHUNK):
    stop = min(start + CHUNK, n_cells)
    lo, hi = int(indptr[start]), int(indptr[stop])
    d = X_data[lo:hi]
    idx = X_ind[lo:hi].astype(np.int64)
    nnz_per_row = np.diff(indptr[start:stop + 1]).astype(np.int64)
    rows = np.repeat(np.arange(start, stop, dtype=np.int64), nnz_per_row)

    neg_entries += int((d < 0).sum())
    explicit_zero_entries += int((d == 0).sum())
    keep = d > 0
    idx_k = idx[keep]
    rows_k = rows[keep]

    flat = bcodes[rows_k] * n_genes + idx_k
    det += np.bincount(flat, minlength=n_blocks * n_genes)

    for j, gi in enumerate(probe_idx):
        m = idx_k == gi
        probe_hit[rows_k[m], j] = True

det = det.reshape(n_blocks, n_genes)
block_n = np.bincount(bcodes, minlength=n_blocks).astype(np.int64)
rate = det / block_n[:, None]

# pooled "elsewhere" rate: leave-this-block-out, cell weighted
tot_det = det.sum(axis=0)
tot_n = block_n.sum()
elsewhere = (tot_det[None, :] - det) / (tot_n - block_n)[:, None]

# leave-one-out median across the other blocks (18 of the 19), per gene
loo_med = np.zeros_like(rate)
for b in range(n_blocks):
    other = np.delete(np.arange(n_blocks), b)
    loo_med[b] = np.median(rate[other], axis=0)

cohort_med = np.median(rate, axis=0)  # all-block median, for reference

np.savez_compressed(
    OUT + "blockqc_cell_flags.npz",
    probes=np.array(PROBES), probe_hit=probe_hit, bcodes=bcodes, blocks=blocks,
)

pd.DataFrame(rate.T, index=genes, columns=blocks).to_csv(OUT + "blockqc_detection_by_block.csv")

# ---- per block dropout / inflation counts -------------------------------------------------
ELSEWHERE_MIN = 0.05     # gene must be detected above 5 pct in the rest of the cohort
INFLATE_FOLD = 5.0       # fold over the leave-one-out cohort median
INFLATE_FLOOR = 0.05     # and reach 5 pct absolute in this block, so tiny rates cannot qualify

rows_out = []
detail = {}
for b in range(n_blocks):
    drop = (det[b] == 0) & (elsewhere[b] > ELSEWHERE_MIN)
    # binomial expectation for an exact zero under the elsewhere rate, given this block's n
    p = np.clip(elsewhere[b], 0, 1)
    exp_zero = float(np.sum((1.0 - p[elsewhere[b] > ELSEWHERE_MIN]) ** block_n[b]))
    infl = (rate[b] > INFLATE_FOLD * loo_med[b]) & (rate[b] >= INFLATE_FLOOR) & (loo_med[b] > 0)
    infl_any = (rate[b] >= INFLATE_FLOOR) & (
        (loo_med[b] == 0) | (rate[b] > INFLATE_FOLD * loo_med[b])
    )
    rows_out.append(dict(
        block=blocks[b], n_cells=int(block_n[b]),
        n_genes_detected=int((det[b] > 0).sum()),
        median_detection=float(np.median(rate[b])),
        n_dropout=int(drop.sum()),
        expected_dropout_binomial=round(exp_zero, 3),
        n_inflated_5x=int(infl.sum()),
        n_inflated_5x_incl_zero_median=int(infl_any.sum()),
        max_fold=float(np.nanmax(np.where(loo_med[b] > 0, rate[b] / np.maximum(loo_med[b], 1e-12), np.nan))),
    ))
    detail[blocks[b]] = dict(
        dropout_genes=sorted(genes[drop].tolist()),
        inflated_genes=sorted(genes[infl].tolist()),
    )

summary = pd.DataFrame(rows_out).sort_values("n_dropout", ascending=False)
summary.to_csv(OUT + "blockqc_block_summary.csv", index=False)

with open(OUT + "blockqc_failing_genes.json", "w") as fh:
    json.dump(detail, fh, indent=1)

pd.set_option("display.width", 250)
print("explicit zero entries in X:", explicit_zero_entries, " negative entries:", neg_entries)
print()
print(summary.to_string(index=False))
print()

# probe genes, block by block
print("probe gene detection rate by block (pct)")
pr = pd.DataFrame(
    {p: rate[:, probe_idx[j]] * 100 for j, p in enumerate(PROBES)}, index=blocks
).round(2)
pr["n_cells"] = block_n
print(pr.to_string())
print()
print("cohort median across blocks (pct):")
print({p: round(float(cohort_med[probe_idx[j]] * 100), 2) for j, p in enumerate(PROBES)})
print("pooled whole-object (pct):")
print({p: round(float(tot_det[probe_idx[j]] / tot_n * 100), 2) for j, p in enumerate(PROBES)})
f.close()
