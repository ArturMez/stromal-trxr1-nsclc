"""
Block QC step 2: composition- and depth-adjusted block fault test, for all 19 blocks.

Why adjusted. A raw "gene is exactly zero here but detected elsewhere" test confounds three
different things:
  (a) a genuine block-level feature-mapping fault,
  (b) a block that contains no cells of the compartment that expresses the gene
      (Leader_Merad 10x 3' v1 is immune only: no epithelium, no stroma),
  (c) a block whose library depth or chemistry is far from the cohort (the two Smart-seq2
      blocks sit at 640k-870k counts per cell against 1-4k for the droplet blocks).
The test removes (b) and (c) before judging, so what remains is (a).

Expectation for gene g in block b:
    expected[b,g] = sum_c  w[b,c] * min(1, s[b,c] * r[-b,c,g])
    w[b,c]    block b's fraction of cells in level1 compartment c (unlabelled cells are their
              own group)
    r[-b,c,g] detection of gene g in compartment c pooled over every OTHER block. If no other
              block has that compartment (true for the unlabelled group, which is EMTAB6653
              alone) the cohort-marginal elsewhere rate is used instead, since compartment
              conditioning is not available there.
    s[b,c]    block b's sensitivity in compartment c: the MEDIAN over informative genes of
              observed/unscaled-expected. Median, so the block's own zeroed and inflated genes
              cannot move it.

Rule, applied identically to all 19 blocks:
    dropout gene   observed count exactly 0 across the whole block AND expected > 5 pct
    inflated gene  observed rate >= 5 pct AND observed rate > 5x expected
    BLOCK FAILS    dropout genes >= 10   (10 genes = 0.52 pct of the 1,906 gene panel)

The inflation count is recorded as supporting evidence but is deliberately NOT an independent
trigger: platform chemistry alone (full length Smart-seq2, BD Rhapsody, GEXSCOPE) moves
relative detection by more than 5-fold in blocks that show no other sign of fault, so an
inflation-only trigger is not specific. This is stated in the uns record as a known limit.

Read only. Writes CSV/JSON. The flag is written by blockqc_03_write_flag.py.
Reads stage1a/panel_object.h5ad; writes blockqc_det_block_comp_gene.npz, blockqc_adjusted_summary.csv and
blockqc_adjusted_detail.json to stage1a/. Run by atlas/01_arms_luca_lambrechts.py (step blockqc_02), no
options. No random numbers are drawn. The same rule is applied again to the 21 and 24 block objects by
stage2_merge_02_blockqc.py.
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
CHUNK = 100_000
EXPECTED_MIN = 0.05
INFLATE_FOLD = 5.0
INFLATE_FLOOR = 0.05
FAIL_N_DROPOUT = 10
MIN_CELLS_FOR_S = 50
S_INFORMATIVE = 0.02

f = h5py.File(H5, "r")
genes = np.array([g.decode() for g in f["var/gene_symbol"][:]])
n_genes = genes.size
blocks = np.array([b.decode() for b in f["obs/block_id/categories"][:]])
comps = np.array([c.decode() for c in f["obs/level1_compartment/categories"][:]])
comps_disp = np.array([c if c else "unlabelled" for c in comps])
bcodes = f["obs/block_id/codes"][:].astype(np.int64)
ccodes = f["obs/level1_compartment/codes"][:].astype(np.int64)
nB, nC, nCells = blocks.size, comps.size, bcodes.size
indptr = f["X/indptr"][:]
X_data, X_ind = f["X/data"], f["X/indices"]

det = np.zeros(nB * nC * n_genes, dtype=np.int64)
key = bcodes * nC + ccodes
for start in range(0, nCells, CHUNK):
    stop = min(start + CHUNK, nCells)
    lo, hi = int(indptr[start]), int(indptr[stop])
    d, idx = X_data[lo:hi], X_ind[lo:hi].astype(np.int64)
    rows = np.repeat(np.arange(start, stop, dtype=np.int64),
                     np.diff(indptr[start:stop + 1]).astype(np.int64))
    keep = d > 0
    det += np.bincount(key[rows[keep]] * n_genes + idx[keep], minlength=nB * nC * n_genes)

det = det.reshape(nB, nC, n_genes)
ncell = np.bincount(key, minlength=nB * nC).reshape(nB, nC).astype(np.int64)
assert ncell.sum() == nCells, "cell count lost"
assert det.sum() == X_data.shape[0], "nnz lost"
np.savez_compressed(OUT + "blockqc_det_block_comp_gene.npz",
                    det=det, ncell=ncell, blocks=blocks, comps=comps, genes=genes)

det_bg = det.sum(axis=1)                       # block x gene
n_b_all = ncell.sum(axis=1)                    # block
rate_bg = det_bg / n_b_all[:, None]
tot_c, totn_c = det.sum(axis=0), ncell.sum(axis=0)
tot_g, totn = det_bg.sum(axis=0), ncell.sum()

rows_out, detail = [], {}
for b in range(nB):
    n_b = n_b_all[b]
    marginal_else = (tot_g - det_bg[b]) / (totn - n_b)
    denom = totn_c - ncell[b]
    r_else = np.zeros((nC, n_genes))
    for c in range(nC):
        r_else[c] = (tot_c[c] - det[b, c]) / denom[c] if denom[c] > 0 else marginal_else

    s = np.full(nC, np.nan)
    for c in range(nC):
        if ncell[b, c] < MIN_CELLS_FOR_S:
            continue
        m = r_else[c] > S_INFORMATIVE
        if m.sum() >= 20:
            s[c] = np.median((det[b, c] / ncell[b, c])[m] / r_else[c][m])
    big = ncell[b] >= MIN_CELLS_FOR_S
    ok = big & np.isfinite(s)
    s_block = float(np.average(s[ok], weights=ncell[b][ok])) if ok.any() else 1.0
    s = np.where(np.isfinite(s), s, s_block)

    w = ncell[b] / n_b
    expected = np.zeros(n_genes)
    for c in range(nC):
        if w[c] > 0:
            expected += w[c] * np.minimum(1.0, s[c] * r_else[c])

    obs = rate_bg[b]
    drop = (det_bg[b] == 0) & (expected > EXPECTED_MIN)
    infl = (obs >= INFLATE_FLOOR) & (obs > INFLATE_FOLD * np.maximum(expected, 1e-12))
    fold = np.divide(obs, expected, out=np.zeros(n_genes), where=expected > 1e-6)

    rows_out.append(dict(
        block=blocks[b], n_cells=int(n_b), sensitivity=round(s_block, 3),
        n_dropout_adj=int(drop.sum()), n_inflated_adj=int(infl.sum()),
        worst_dropout_expected_pct=round(float(expected[drop].max() * 100), 1) if drop.any() else 0.0,
        max_fold_adj=round(float(fold[infl].max()), 1) if infl.any() else 0.0,
        fails=bool(drop.sum() >= FAIL_N_DROPOUT),
    ))
    detail[blocks[b]] = dict(
        n_cells=int(n_b), sensitivity_block=round(s_block, 3),
        sensitivity_by_compartment={comps_disp[c]: round(float(s[c]), 3)
                                    for c in range(nC) if ncell[b, c] > 0},
        n_dropout_adj=int(drop.sum()), n_inflated_adj=int(infl.sum()),
        dropout_genes_expected_pct={g: round(float(e * 100), 1) for g, e in
                                    sorted(zip(genes[drop], expected[drop]), key=lambda t: -t[1])},
        inflated_genes_obs_exp_fold={g: [round(float(o * 100), 1), round(float(e * 100), 1),
                                         round(float(fd), 1)]
                                     for g, o, e, fd in sorted(
                                         zip(genes[infl], obs[infl], expected[infl], fold[infl]),
                                         key=lambda t: -t[3])[:120]},
    )

summary = pd.DataFrame(rows_out).sort_values(
    ["n_dropout_adj", "n_inflated_adj"], ascending=False).reset_index(drop=True)
summary.to_csv(OUT + "blockqc_adjusted_summary.csv", index=False)
with open(OUT + "blockqc_adjusted_detail.json", "w") as fh:
    json.dump(detail, fh, indent=1)

pd.set_option("display.width", 250)
print("COMPOSITION- AND DEPTH-ADJUSTED BLOCK TEST, all 19 blocks")
print(summary.to_string(index=False))
print()
print("threshold invariance of the dropout rule:")
for k in [2, 3, 5, 10, 15, 20, 25]:
    fl = summary[summary.n_dropout_adj >= k].block.tolist()
    print(f"  fail if dropouts >= {k:3d}:  {[str(x) for x in fl]}")
print()

# ---- model free confirmation: same compartment, this block vs the rest --------------------
print("MODEL FREE CHECK for the failing block: detection inside the SAME compartment (pct)")
b = int(np.where(blocks == "Adams_Kaminski_2020_10x 3' v2")[0][0])
probe = [("CD68", "myeloid"), ("MS4A7", "myeloid"), ("CYBB", "myeloid"), ("IFI30", "myeloid"),
         ("IGKC", "B_plasma"), ("LAT", "T_NK"), ("SH2D1A", "T_NK"), ("PRDX2", "myeloid"),
         ("PRDX2", "T_NK"), ("BLOC1S1", "myeloid"), ("CHMP1B", "myeloid"),
         ("TXNRD1", "myeloid"), ("TXNRD1", "fibroblast"), ("TXNRD1", "epithelial_malignant")]
recs = []
for g, cname in probe:
    gi = int(np.where(genes == g)[0][0]); ci = int(np.where(comps == cname)[0][0])
    n_here = ncell[b, ci]; d_here = det[b, ci, gi]
    n_else = totn_c[ci] - n_here; d_else = tot_c[ci, gi] - d_here
    recs.append(dict(gene=g, compartment=cname, cells_in_block=int(n_here),
                     pct_in_block=round(100 * d_here / max(n_here, 1), 2),
                     cells_elsewhere=int(n_else),
                     pct_elsewhere=round(100 * d_else / max(n_else, 1), 2)))
print(pd.DataFrame(recs).to_string(index=False))
f.close()
