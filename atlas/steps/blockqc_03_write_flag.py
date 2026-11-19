"""
Block QC step 3: write obs.block_qc_pass and uns['block_qc'] into panel_object.h5ad.

This is a FLAG, not a filter. No cell is removed, no count is altered. Only two new objects are
added to the file: the obs column and the uns record. Downstream analyses opt in by filtering
on obs.block_qc_pass.

Rule written to the object: a gene counts as a dropout in a block when its count is 0 across the whole
block while its composition- and depth-adjusted expected detection there exceeds 5 percent; a block fails
(block_qc_pass False for all its cells) at 10 dropout genes. Reads blockqc_adjusted_summary.csv and
blockqc_adjusted_detail.json (step 2); modifies stage1a/panel_object.h5ad in place after saving a backup
(panel_object.h5ad.prebockqc.bak). Run by atlas/01_arms_luca_lambrechts.py (step blockqc_03), no options.
No random numbers are drawn.
"""
import datetime as dt
import json
import shutil
import os
import numpy as np
import pandas as pd
import h5py

import os as _os  # noqa: E402
import sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # atlas/
from common import MD  # noqa: E402  paths from config.json, see atlas/common.py
H5 = MD + "/stage1a/panel_object.h5ad"
BAK = H5 + ".prebockqc.bak"
from common import MD  # noqa: E402  paths from config.json, see atlas/common.py
OUT = MD + "/stage1a/"
FAIL_N_DROPOUT = 10

summary = pd.read_csv(OUT + "blockqc_adjusted_summary.csv")
detail = json.load(open(OUT + "blockqc_adjusted_detail.json"))

RULE = (
    "block_qc_pass is False for every cell of a block that fails the block-level feature "
    "quantification test. A block FAILS if it has 10 or more dropout genes, where a dropout "
    "gene is a panel gene with EXACTLY ZERO counts summed over every cell of the block while "
    "its composition- and depth-adjusted expected detection in that block exceeds 5 percent. "
    "10 genes is 0.52 percent of the 1,906 gene panel. Threshold applied identically to all 19 "
    "blocks."
)
METHOD = (
    "Expected detection for gene g in block b is sum over level1 compartments c of "
    "w[b,c] * min(1, s[b,c] * r[-b,c,g]). w[b,c] is block b's fraction of cells in compartment "
    "c, so a block that lacks a compartment is not penalised for lacking that compartment's "
    "genes. r[-b,c,g] is detection of g in compartment c pooled over all OTHER blocks; where no "
    "other block has that compartment (the unlabelled group, EMTAB6653 alone) the cohort "
    "marginal is substituted. s[b,c] is the block's sensitivity in compartment c, the median "
    "over genes with elsewhere detection above 2 percent of observed/unscaled-expected, so "
    "library depth and chemistry are divided out and the block's own zeroed and inflated genes "
    "cannot move it. Unlabelled cells form their own compartment group."
)
INVARIANCE = (
    "The verdict does not depend on the cut. Any dropout threshold from 3 to 20 genes selects "
    "the same single block. Adams_Kaminski_2020_10x 3' v2 has 20 dropout genes; the next "
    "highest block in the cohort has 2 (Mayr_Schiller_2020_Drop-seq) and 16 blocks have 0."
)
LIMITS = (
    "Gene inflation is recorded per block but is deliberately NOT an independent trigger. "
    "Platform chemistry alone moves relative detection by more than 5-fold in blocks with no "
    "other sign of fault (full length Smart-seq2, BD Rhapsody, GEXSCOPE), so an inflation-only "
    "rule is not specific. A block whose fault inflated genes without zeroing any would "
    "therefore not be caught by this rule. The dropout arm is the discriminating test because "
    "an exact zero across a whole block is not reachable by sampling: for the failing block the "
    "least extreme dropout gene has expected detection 5.1 percent over 131,092 cells."
)

fail_mask = summary.n_dropout_adj >= FAIL_N_DROPOUT
failing = summary.loc[fail_mask, "block"].tolist()
reasons = []
for blk in failing:
    d = detail[blk]
    dg = d["dropout_genes_expected_pct"]
    top = list(dg.items())[:8]
    reasons.append(
        f"{blk}: {d['n_dropout_adj']} panel genes carry exactly zero counts across all "
        f"{d['n_cells']:,} cells of the block while their composition- and depth-adjusted "
        f"expected detection is above 5 percent. Worst offenders (gene, expected detection): "
        + ", ".join(f"{g} {v}pct" for g, v in top)
        + f". Same-compartment, model-free confirmation: CD68 is detected in 0.00 percent of "
        f"the block's 91,418 myeloid cells against 72.61 percent of the 276,703 myeloid cells "
        f"in the rest of the cohort; IGKC 0.00 percent of 3,752 B_plasma cells against 88.41 "
        f"percent elsewhere. In the same block {d['n_inflated_adj']} genes are inflated more "
        f"than 5-fold over expectation (max {summary.loc[summary.block == blk, 'max_fold_adj'].iloc[0]}x), "
        f"and the inflated set is dominated by very long genes (AGBL4, DLG2, ASTN2, MAGI2, "
        f"ANKS1B, CACNB2, IMMP2L), a pattern consistent with intron-inclusive quantification "
        f"against a different annotation. TXNRD1 sits at 75.02 percent in this block's myeloid "
        f"cells against 21.18 percent in myeloid cells elsewhere and 91.58 percent in its "
        f"fibroblasts against 15.42 percent elsewhere. Simultaneous complete absence of "
        f"broadly expressed genes and many-fold inflation of others is a block-level feature "
        f"quantification fault, not biology."
    )

if not os.path.exists(BAK):
    print("backing up to", BAK)
    shutil.copy2(H5, BAK)
    print("backup done")

with h5py.File(H5, "r+") as f:
    blocks = np.array([b.decode() for b in f["obs/block_id/categories"][:]])
    bcodes = f["obs/block_id/codes"][:]
    order = list(summary.set_index("block").loc[blocks].index)
    assert order == list(blocks)
    pass_by_block = np.array([b not in failing for b in blocks], dtype=bool)
    flag = pass_by_block[bcodes]
    assert flag.shape[0] == bcodes.shape[0]

    if "block_qc_pass" in f["obs"]:
        del f["obs"]["block_qc_pass"]
    ds = f["obs"].create_dataset("block_qc_pass", data=flag, dtype=bool, compression="gzip")
    ds.attrs["encoding-type"] = "array"
    ds.attrs["encoding-version"] = "0.2.0"
    co = list(f["obs"].attrs["column-order"])
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

    sstr("rule", RULE)
    sstr("method", METHOD)
    sstr("threshold_invariance", INVARIANCE)
    sstr("known_limits", LIMITS)
    narr("threshold_dropout_genes", FAIL_N_DROPOUT, "int64")
    narr("threshold_expected_detection_pct", 5.0, "float64")
    narr("inflation_fold_recorded", 5.0, "float64")
    narr("inflation_floor_pct", 5.0, "float64")
    sarr("blocks", blocks)
    narr("n_cells_by_block", summary.set_index("block").loc[blocks, "n_cells"].values, "int64")
    narr("n_dropout_genes_by_block",
         summary.set_index("block").loc[blocks, "n_dropout_adj"].values, "int64")
    narr("n_inflated_genes_by_block",
         summary.set_index("block").loc[blocks, "n_inflated_adj"].values, "int64")
    narr("block_qc_pass_by_block", pass_by_block, bool)
    sarr("failing_blocks", failing)
    sarr("failing_reason", reasons)
    sarr("dropout_genes_failing_block",
         list(detail[failing[0]]["dropout_genes_expected_pct"].keys()) if failing else [])
    sstr("applies_to", "flag only. No cell was removed and no count altered. Downstream "
                       "analyses must filter obs.block_qc_pass themselves.")
    sstr("script", "atlas/steps/blockqc_01_detection_by_block.py, "
                   "blockqc_02_rule_and_flag.py, blockqc_03_write_flag.py")
    sstr("written_utc", dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"))

print("failing blocks:", failing)
print("cells flagged False:", int((~flag).sum()), "of", flag.shape[0])

# verify by reopening through anndata
import anndata as ad
a = ad.read_h5ad(H5, backed="r")
print()
print("reopened with anndata:", a.shape)
print("obs.block_qc_pass dtype:", a.obs.block_qc_pass.dtype)
print(a.obs.groupby("block_id", observed=True).block_qc_pass.agg(["size", "mean"]).to_string())
print()
print("uns['block_qc'] keys:", sorted(a.uns["block_qc"].keys()))
print("rule:", a.uns["block_qc"]["rule"][:200], "...")
