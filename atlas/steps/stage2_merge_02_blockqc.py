#!/usr/bin/env python
"""
Atlas stage 02 (merge_blockqc) and stage 04 (blockqc_v2): rerun the block QC dropout test on a draft
object, all blocks, and write obs.block_qc_pass plus a refreshed uns['block_qc'].

Paper items: part of the atlas behind Figure 2, Figure 3, Supplementary Figures 2 to 4 and Supplementary Data 1.
Run: python atlas/steps/stage2_merge_02_blockqc.py   (normally through atlas/02_arm_emtab13526.py on the
21-block draft, and through atlas/04_build_objects.py on the 24-block panel_object_v2; the environment
variables BLOCKQC_DRAFT and BLOCKQC_OUT select the object and output folder). No random numbers are drawn.

This is the identical rule and the identical estimator as blockqc_02_rule_and_flag.py /
blockqc_03_write_flag.py, with three differences:
  - it runs on the merged object that contains EMTAB13526 (and, in stage 04, the wave-1 arms), so
    the new blocks are judged on the same rule rather than trusted;
  - every number in the uns record, including the ones inside the human-readable failing
    reason, is recomputed from the object being described;
  - the unadjusted dropout count is reported alongside, as supporting information only.

Rule, unchanged:
    dropout gene   zero counts across the whole block AND adjusted expected detection > 5 pct
    inflated gene  observed >= 5 pct AND observed > 5x adjusted expected  (recorded, not a trigger)
    BLOCK FAILS    dropout genes >= 10
"""
import datetime as dt
import json
import sys
import time

import os
import h5py
import numpy as np
import pandas as pd

T0 = time.time()
import os as _os  # noqa: E402
import sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # atlas/
from common import MD  # noqa: E402  paths from config.json, see atlas/common.py
BASE = MD + "/"
DRAFT = os.environ.get("BLOCKQC_DRAFT", BASE + "stage2/panel_object_draft.h5ad")
OUT = os.environ.get("BLOCKQC_OUT", BASE + "stage2/")
CHUNK = 100_000
EXPECTED_MIN = 0.05
INFLATE_FOLD = 5.0
INFLATE_FLOOR = 0.05
FAIL_N_DROPOUT = 10
MIN_CELLS_FOR_S = 50
S_INFORMATIVE = 0.02


def log(*a):
    """Print a message with the elapsed time."""
    print(f"[{time.time()-T0:7.1f}s]", *a, flush=True)


f = h5py.File(DRAFT, "r+")
genes = np.array([g.decode() for g in f["var/gene_symbol"][:]])
n_genes = genes.size
blocks = np.array([b.decode() for b in f["obs/block_id/categories"][:]])
comps = np.array([c.decode() for c in f["obs/level1_compartment/categories"][:]])
comps_disp = np.array([c if c else "unlabelled" for c in comps])
bcodes = f["obs/block_id/codes"][:].astype(np.int64)
ccodes = f["obs/level1_compartment/codes"][:].astype(np.int64)
nB, nC, nCells = blocks.size, comps.size, bcodes.size
log(f"draft {nCells:,} cells, {n_genes} genes, {nB} blocks, {nC} compartment groups")

indptr = f["X/indptr"][:].astype(np.int64)
X_data, X_ind = f["X/data"], f["X/indices"]

det = np.zeros(nB * nC * n_genes, dtype=np.int64)
key = bcodes * nC + ccodes
for start in range(0, nCells, CHUNK):
    stop = min(start + CHUNK, nCells)
    lo, hi = int(indptr[start]), int(indptr[stop])
    d, idx = X_data[lo:hi], X_ind[lo:hi].astype(np.int64)
    rows = np.repeat(np.arange(start, stop, dtype=np.int64),
                     np.diff(indptr[start:stop + 1]))
    keep = d > 0
    det += np.bincount(key[rows[keep]] * n_genes + idx[keep], minlength=nB * nC * n_genes)
    if (start // CHUNK) % 5 == 0:
        log(f"  scanned {stop:,} / {nCells:,}")

det = det.reshape(nB, nC, n_genes)
ncell = np.bincount(key, minlength=nB * nC).reshape(nB, nC).astype(np.int64)
assert ncell.sum() == nCells, "cell count lost"
assert det.sum() == X_data.shape[0], f"nnz lost: {det.sum()} vs {X_data.shape[0]}"
log("scan complete, det.sum() == nnz and ncell.sum() == n_cells PASS")
np.savez_compressed(OUT + "stage2_blockqc_det_block_comp_gene.npz",
                    det=det, ncell=ncell, blocks=blocks, comps=comps, genes=genes)

det_bg = det.sum(axis=1)
n_b_all = ncell.sum(axis=1)
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

    obs_r = rate_bg[b]
    drop = (det_bg[b] == 0) & (expected > EXPECTED_MIN)
    infl = (obs_r >= INFLATE_FLOOR) & (obs_r > INFLATE_FOLD * np.maximum(expected, 1e-12))
    fold = np.divide(obs_r, expected, out=np.zeros(n_genes), where=expected > 1e-6)
    # unadjusted variant, supporting information only
    drop_raw = (det_bg[b] == 0) & (marginal_else > EXPECTED_MIN)

    rows_out.append(dict(
        block=blocks[b], n_cells=int(n_b), sensitivity=round(s_block, 3),
        n_dropout_adj=int(drop.sum()), n_dropout_unadjusted=int(drop_raw.sum()),
        n_inflated_adj=int(infl.sum()),
        worst_dropout_expected_pct=round(float(expected[drop].max() * 100), 1) if drop.any() else 0.0,
        max_fold_adj=round(float(fold[infl].max()), 1) if infl.any() else 0.0,
        fails=bool(drop.sum() >= FAIL_N_DROPOUT),
    ))
    detail[blocks[b]] = dict(
        n_cells=int(n_b), sensitivity_block=round(s_block, 3),
        cells_by_compartment={comps_disp[c]: int(ncell[b, c]) for c in range(nC) if ncell[b, c] > 0},
        sensitivity_by_compartment={comps_disp[c]: round(float(s[c]), 3)
                                    for c in range(nC) if ncell[b, c] > 0},
        n_dropout_adj=int(drop.sum()), n_dropout_unadjusted=int(drop_raw.sum()),
        n_inflated_adj=int(infl.sum()),
        dropout_genes_expected_pct={g: round(float(e * 100), 1) for g, e in
                                    sorted(zip(genes[drop], expected[drop]), key=lambda t: -t[1])},
        inflated_genes_obs_exp_fold={g: [round(float(o * 100), 1), round(float(e * 100), 1),
                                         round(float(fd), 1)]
                                     for g, o, e, fd in sorted(
                                         zip(genes[infl], obs_r[infl], expected[infl], fold[infl]),
                                         key=lambda t: -t[3])[:120]},
    )

summary = pd.DataFrame(rows_out).sort_values(
    ["n_dropout_adj", "n_inflated_adj"], ascending=False).reset_index(drop=True)
summary.to_csv(OUT + "stage2_blockqc_adjusted_summary.csv", index=False)
with open(OUT + "stage2_blockqc_adjusted_detail.json", "w") as fh:
    json.dump(detail, fh, indent=1)

pd.set_option("display.width", 250)
print()
print(f"COMPOSITION- AND DEPTH-ADJUSTED BLOCK TEST, all {nB} blocks")
print(summary.to_string(index=False))
print()
print("threshold invariance of the dropout rule:")
for k in [2, 3, 5, 10, 15, 20, 25]:
    fl = summary[summary.n_dropout_adj >= k].block.tolist()
    print(f"  fail if dropouts >= {k:3d}:  {[str(x) for x in fl]}")
print()

failing = summary.loc[summary.n_dropout_adj >= FAIL_N_DROPOUT, "block"].tolist()
print("failing blocks:", failing)

# ---------------------------------------------------------------- model-free numbers
def pct(b, cname, gname):
    """Cells and detection percent of a gene in a compartment, in this block and in all other blocks (4-tuple)."""
    gi = int(np.where(genes == gname)[0][0])
    ci = int(np.where(comps == cname)[0][0])
    n_here = int(ncell[b, ci]); d_here = int(det[b, ci, gi])
    n_else = int(totn_c[ci] - n_here); d_else = int(tot_c[ci, gi] - d_here)
    return (n_here, 100 * d_here / max(n_here, 1), n_else, 100 * d_else / max(n_else, 1))


print()
print("MODEL-FREE CHECK, same compartment, this block against the rest of the cohort")
probe = [("CD68", "myeloid"), ("IGKC", "B_plasma"), ("LAT", "T_NK"), ("PRDX2", "myeloid"),
         ("TXNRD1", "myeloid"), ("TXNRD1", "fibroblast"), ("TXNRD1", "epithelial_malignant")]
recs = []
for b in range(nB):
    if blocks[b] not in failing:
        continue
    for gname, cname in probe:
        nh, ph, ne, pe = pct(b, cname, gname)
        recs.append(dict(block=blocks[b], gene=gname, compartment=cname, cells_in_block=nh,
                         pct_in_block=round(ph, 2), cells_elsewhere=ne,
                         pct_elsewhere=round(pe, 2)))
if recs:
    print(pd.DataFrame(recs).to_string(index=False))

# ---------------------------------------------------------------- uns text, recomputed
RULE = (
    "block_qc_pass is False for every cell of a block that fails the block-level feature "
    "quantification test. A block FAILS if it has 10 or more dropout genes, where a dropout "
    "gene is a panel gene with EXACTLY ZERO counts summed over every cell of the block while "
    f"its composition- and depth-adjusted expected detection in that block exceeds 5 percent. "
    f"10 genes is {100*10/n_genes:.2f} percent of the {n_genes} gene panel. Threshold applied "
    f"identically to all {nB} blocks."
)
METHOD = (
    "Expected detection for gene g in block b is sum over level1 compartments c of "
    "w[b,c] * min(1, s[b,c] * r[-b,c,g]). w[b,c] is block b's fraction of cells in compartment "
    "c, so a block that lacks a compartment is not penalised for lacking that compartment's "
    "genes. r[-b,c,g] is detection of g in compartment c pooled over all OTHER blocks; where no "
    "other block has that compartment the cohort marginal is substituted. s[b,c] is the block's "
    "sensitivity in compartment c, the median over genes with elsewhere detection above 2 "
    "percent of observed/unscaled-expected, so library depth and chemistry are divided out and "
    "the block's own zeroed and inflated genes cannot move it. Cells with no compartment label "
    "form their own group."
)
inv_lines = []
for k in [3, 5, 10, 15, 20]:
    fl = summary[summary.n_dropout_adj >= k].block.tolist()
    inv_lines.append(f"{k}: {sorted(str(x) for x in fl)}")
srt = summary.sort_values("n_dropout_adj", ascending=False)
INVARIANCE = (
    "Verdict against threshold, recomputed on this object. Dropout cut: failing blocks. "
    + "; ".join(inv_lines) + ". "
    + f"Highest dropout count {int(srt.n_dropout_adj.iloc[0])} ({srt.block.iloc[0]}), "
    + f"next highest {int(srt.n_dropout_adj.iloc[1])} ({srt.block.iloc[1]}), "
    + f"{int((summary.n_dropout_adj == 0).sum())} of {nB} blocks have 0."
)
LIMITS = (
    "Gene inflation is recorded per block but is deliberately NOT an independent trigger. "
    "Platform chemistry alone moves relative detection by more than 5-fold in blocks with no "
    "other sign of fault (full length Smart-seq2, BD Rhapsody, GEXSCOPE), so an inflation-only "
    "rule is not specific. A block whose fault inflated genes without zeroing any would "
    "therefore not be caught by this rule. The dropout arm is the discriminating test because "
    "an exact zero across a whole block is not reachable by sampling."
)

reasons = []
for blk in failing:
    b = int(np.where(blocks == blk)[0][0])
    d = detail[blk]
    dg = d["dropout_genes_expected_pct"]
    top = list(dg.items())[:8]
    least = min(dg.values()) if dg else 0.0
    mf = []
    for gname, cname in [("CD68", "myeloid"), ("IGKC", "B_plasma")]:
        nh, ph, ne, pe = pct(b, cname, gname)
        mf.append(f"{gname} is detected in {ph:.2f} percent of the block's {nh:,} {cname} cells "
                  f"against {pe:.2f} percent of the {ne:,} {cname} cells in the rest of the cohort")
    tx = []
    for cname in ["myeloid", "fibroblast"]:
        nh, ph, ne, pe = pct(b, cname, "TXNRD1")
        tx.append(f"{ph:.2f} percent in this block's {nh:,} {cname} cells against {pe:.2f} "
                  f"percent in the {ne:,} elsewhere")
    reasons.append(
        f"{blk}: {d['n_dropout_adj']} panel genes carry exactly zero counts across all "
        f"{d['n_cells']:,} cells of the block while their composition- and depth-adjusted "
        f"expected detection is above 5 percent. Worst offenders (gene, expected detection): "
        + ", ".join(f"{g} {v}pct" for g, v in top)
        + f". The least extreme dropout gene still has expected detection {least} percent over "
        f"{d['n_cells']:,} cells, which sampling cannot produce. Same-compartment, model-free "
        f"confirmation: " + "; ".join(mf) + ". "
        f"In the same block {d['n_inflated_adj']} genes are inflated more than 5-fold over "
        f"expectation (max {summary.loc[summary.block == blk, 'max_fold_adj'].iloc[0]}x). "
        f"TXNRD1 sits at " + " and ".join(tx) + ". "
        f"Simultaneous complete absence of broadly expressed genes and many-fold inflation of "
        f"others is a block-level feature quantification fault, not biology."
    )

# ---------------------------------------------------------------- write the flag
pass_by_block = np.array([b not in failing for b in blocks], dtype=bool)
flag = pass_by_block[bcodes]
assert flag.shape[0] == nCells

if "block_qc_pass" in f["obs"]:
    del f["obs"]["block_qc_pass"]
ds = f["obs"].create_dataset("block_qc_pass", data=flag, dtype=bool, compression="gzip")
ds.attrs["encoding-type"] = "array"
ds.attrs["encoding-version"] = "0.2.0"
co = [c.decode() if isinstance(c, bytes) else str(c) for c in f["obs"].attrs["column-order"]]
if "block_qc_pass" not in co:
    co.append("block_qc_pass")
f["obs"].attrs["column-order"] = np.array(co, dtype=object)

if "block_qc" in f["uns"]:
    del f["uns"]["block_qc"]
g = f["uns"].create_group("block_qc")
g.attrs["encoding-type"] = "dict"
g.attrs["encoding-version"] = "0.1.0"


def sstr(name, val):
    """Write a string value into uns['block_qc'] with the anndata string encoding."""
    d = g.create_dataset(name, data=str(val))
    d.attrs["encoding-type"] = "string"
    d.attrs["encoding-version"] = "0.2.0"


def sarr(name, vals):
    """Write a string array into uns['block_qc'] with the anndata string-array encoding."""
    d = g.create_dataset(name, data=np.array([str(v) for v in vals], dtype=object),
                         dtype=h5py.special_dtype(vlen=str))
    d.attrs["encoding-type"] = "string-array"
    d.attrs["encoding-version"] = "0.2.0"


def narr(name, vals, dtype):
    """Write a numeric array of the given dtype into uns['block_qc']."""
    d = g.create_dataset(name, data=np.array(vals, dtype=dtype))
    d.attrs["encoding-type"] = "array"
    d.attrs["encoding-version"] = "0.2.0"


S = summary.set_index("block").loc[blocks]
sstr("rule", RULE)
sstr("method", METHOD)
sstr("threshold_invariance", INVARIANCE)
sstr("known_limits", LIMITS)
narr("threshold_dropout_genes", FAIL_N_DROPOUT, "int64")
narr("threshold_expected_detection_pct", 5.0, "float64")
narr("inflation_fold_recorded", 5.0, "float64")
narr("inflation_floor_pct", 5.0, "float64")
sarr("blocks", blocks)
narr("n_cells_by_block", S["n_cells"].values, "int64")
narr("n_dropout_genes_by_block", S["n_dropout_adj"].values, "int64")
narr("n_dropout_genes_unadjusted_by_block", S["n_dropout_unadjusted"].values, "int64")
narr("n_inflated_genes_by_block", S["n_inflated_adj"].values, "int64")
narr("sensitivity_by_block", S["sensitivity"].values, "float64")
narr("block_qc_pass_by_block", pass_by_block, bool)
sarr("compartments", comps_disp)
narr("n_cells_by_block_by_compartment", ncell, "int64")
sarr("failing_blocks", failing)
sarr("failing_reason", reasons)
sarr("dropout_genes_failing_block",
     list(detail[failing[0]]["dropout_genes_expected_pct"].keys()) if failing else [])
sstr("applies_to", "flag only. No cell was removed and no count altered. Downstream "
                   "analyses must filter obs.block_qc_pass themselves.")
sstr("computed_on", f"the object that carries this record: {nCells} cells, {n_genes} genes, "
                    f"{nB} blocks. Every number here, including the numbers inside "
                    f"failing_reason, was recomputed from it.")
sstr("script", "atlas/steps/stage2_merge_02_blockqc.py")
sstr("written_utc", dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"))
f.close()

print()
print("cells flagged block_qc_pass False:", int((~flag).sum()), "of", flag.shape[0])
log("step 2 done")
