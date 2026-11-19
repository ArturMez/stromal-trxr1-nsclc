#!/usr/bin/env python
"""
Atlas stage 02, step merge_verify: verify the draft against the pre-merge backup and against the arm,
then move it into place as stage1a/panel_object.h5ad.

Paper items: part of the atlas behind Figure 2, Figure 3, Supplementary Figures 2 to 4 and Supplementary Data 1.
Run: python atlas/steps/stage2_merge_03_verify.py   (normally through atlas/02_arm_emtab13526.py).
No random numbers are drawn.

Nothing here trusts the merge_add step. Every check re-reads from disk.

  A. shape, gene axis, cell_id uniqueness and form
  B. the first 1,146,594 rows of obs are IDENTICAL to the backup, column by column, values
     and NA pattern, including block_qc_pass
  C. the first 164,040,429 nonzeros of X are IDENTICAL to the backup, value for value, and
     the first 1,146,595 indptr entries are identical
  D. the appended 896,773 rows of X equal arm_emtab13526.h5ad restricted to the merged gene
     axis, compared entry by entry over every cell, not a sample
  E. obs.n_counts_intersect and n_genes_intersect recomputed from the written X for every cell
  F. schema invariants on the merged object

Only if every check passes is the draft moved into stage1a/panel_object.h5ad.
"""
import os
import sys
import time

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem

T0 = time.time()
import os as _os  # noqa: E402
import sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # atlas/
from common import MD  # noqa: E402  paths from config.json, see atlas/common.py
BASE = MD + "/"
S1A = BASE + "stage1a/"
S2 = BASE + "stage2/"
BAK = S1A + "panel_object_prestage2.h5ad"
ARM = S2 + "arm_emtab13526.h5ad"
DRAFT = S2 + "panel_object_draft.h5ad"
TARGET = S1A + "panel_object.h5ad"
CELL_CHUNK = 100_000
NNZ_CHUNK = 20_000_000

FAILS = []


def log(*a):
    """Print a message with the elapsed time."""
    print(f"[{time.time()-T0:7.1f}s]", *a, flush=True)


def check(name, ok, extra=""):
    """Print PASS or FAIL for one check and record failures."""
    print(f"  {'PASS' if ok else 'FAIL'}  {name}{('  ' + extra) if extra else ''}", flush=True)
    if not ok:
        FAILS.append(name)


fb = h5py.File(BAK, "r")
fd = h5py.File(DRAFT, "r")
fa = h5py.File(ARM, "r")

n_bak = int(fb["X"].attrs["shape"][0])
n_arm = int(fa["X"].attrs["shape"][0])
n_new = int(fd["X"].attrs["shape"][0])
n_genes = int(fd["X"].attrs["shape"][1])

print("\n=== A. shape and axes ===")
check("row count is backup + arm", n_new == n_bak + n_arm, f"{n_bak:,} + {n_arm:,} = {n_new:,}")
gb = np.array([g.decode() for g in fb["var/gene_symbol"][:]])
gd = np.array([g.decode() for g in fd["var/gene_symbol"][:]])
check("gene axis identical to the backup, same order", list(gb) == list(gd), f"{len(gd)} genes")

idb = np.array([x.decode() for x in fb["obs/_index"][:]])
idd = np.array([x.decode() for x in fd["obs/_index"][:]])
ida = np.array([x.decode() for x in fa["obs/_index"][:]])
check("existing cell_ids unchanged and in the same order", np.array_equal(idb, idd[:n_bak]))
check("appended cell_ids are the arm's, in the arm's order", np.array_equal(ida, idd[n_bak:]))
u, c = np.unique(idd, return_counts=True)
check("cell_id globally unique", u.size == n_new, f"{u.size:,} distinct of {n_new:,}")
if u.size != n_new:
    print("    collisions:", u[c > 1][:10])
check("cell_id has exactly two pipes on every row",
      bool(np.all(np.char.count(idd, "|") == 2)))
pref = pd.Series(np.char.partition(idd, "|")[:, 0]).value_counts().to_dict()
print("    source prefixes:", pref)

print("\n=== B. pre-existing obs values ===")
ob = read_elem(fb["obs"])
od = read_elem(fd["obs"])
check("obs column set identical", list(ob.columns) == list(od.columns),
      f"{len(od.columns)} columns")
head = od.iloc[:n_bak]
n_changed_cols = 0
for c in ob.columns:
    a = ob[c]
    b = head[c]
    if str(a.dtype) == "category":
        a = a.astype(str)
        b = b.astype(str)
    va, vb = a.values, b.values
    if hasattr(va, "to_numpy"):
        va = va.to_numpy(dtype=object)
    if hasattr(vb, "to_numpy"):
        vb = vb.to_numpy(dtype=object)
    na_a = pd.isna(pd.Series(va))
    na_b = pd.isna(pd.Series(vb))
    same_na = bool((na_a == na_b).all())
    m = (~na_a).values
    if va.dtype.kind == "f" or vb.dtype.kind == "f":
        same_val = bool(np.array_equal(np.asarray(va, dtype="float64")[m],
                                       np.asarray(vb, dtype="float64")[m]))
    else:
        same_val = bool(np.array_equal(np.asarray(va, dtype=object)[m],
                                       np.asarray(vb, dtype=object)[m]))
    if not (same_na and same_val):
        n_changed_cols += 1
        diff = int((~((np.asarray(va, dtype=object) == np.asarray(vb, dtype=object)) |
                      (na_a.values & na_b.values))).sum())
        print(f"    CHANGED column {c}: {diff:,} of {n_bak:,} rows differ")
check("no pre-existing cell's obs value changed in any of the 40 columns", n_changed_cols == 0,
      f"{len(ob.columns)} columns x {n_bak:,} rows compared")
qb = ob["block_qc_pass"].values
qd = head["block_qc_pass"].values
check("block_qc_pass unchanged for every pre-existing cell after the recompute",
      bool(np.array_equal(qb, qd)),
      f"{int((~qd).sum()):,} pre-existing cells still flagged False")

print("\n=== C. pre-existing X ===")
nnz_b = fb["X/data"].shape[0]
ipb = fb["X/indptr"][:].astype(np.int64)
ipd = fd["X/indptr"][:].astype(np.int64)
check("first indptr block identical", bool(np.array_equal(ipb, ipd[:n_bak + 1])),
      f"{n_bak+1:,} entries, closes at {ipb[-1]:,}")
ok_d = ok_i = True
for s in range(0, nnz_b, NNZ_CHUNK):
    e = min(s + NNZ_CHUNK, nnz_b)
    ok_d &= bool(np.array_equal(fb["X/data"][s:e], fd["X/data"][s:e]))
    ok_i &= bool(np.array_equal(fb["X/indices"][s:e], fd["X/indices"][s:e]))
check("X data of the pre-existing cells identical", ok_d, f"{nnz_b:,} nonzeros compared")
check("X indices of the pre-existing cells identical", ok_i)

print("\n=== D. appended X against the arm, entry by entry, all cells ===")
ag = np.array([x.decode() for x in fa["var/_index"][:]])
pos = {g: i for i, g in enumerate(ag)}
lut = np.full(len(ag), -1, dtype=np.int64)
for j, g in enumerate(gd):
    lut[pos[g]] = j
ipa = fa["X/indptr"][:].astype(np.int64)
bad_rows = 0
n_cmp = 0
nc_re = np.zeros(n_arm)
ng_re = np.zeros(n_arm, dtype=np.int64)
for start in range(0, n_arm, CELL_CHUNK):
    stop = min(start + CELL_CHUNK, n_arm)
    la, ha = int(ipa[start]), int(ipa[stop])
    d_a = fa["X/data"][la:ha].astype(np.float32)
    i_a = fa["X/indices"][la:ha]
    rows_a = np.repeat(np.arange(start, stop, dtype=np.int64), np.diff(ipa[start:stop + 1]))
    newcol = lut[i_a]
    keep = newcol >= 0
    exp_d, exp_c, exp_r = d_a[keep], newcol[keep], rows_a[keep]

    ld, hd = int(ipd[n_bak + start]), int(ipd[n_bak + stop])
    got_d = fd["X/data"][ld:hd]
    got_c = fd["X/indices"][ld:hd].astype(np.int64)
    got_r = np.repeat(np.arange(start, stop, dtype=np.int64),
                      np.diff(ipd[n_bak + start:n_bak + stop + 1]))
    if not (exp_d.size == got_d.size and np.array_equal(exp_c, got_c)
            and np.array_equal(exp_d, got_d) and np.array_equal(exp_r, got_r)):
        bad = np.unique(np.concatenate([
            exp_r[:min(exp_r.size, got_r.size)][
                exp_c[:min(exp_c.size, got_c.size)] != got_c[:min(exp_c.size, got_c.size)]],
        ])) if exp_c.size and got_c.size else np.array([start])
        bad_rows += max(1, bad.size)
        print(f"    MISMATCH in cells {start}-{stop}: exp nnz {exp_d.size} got {got_d.size}")
    n_cmp += exp_d.size
    nc_re[start:stop] = np.bincount(exp_r - start, weights=exp_d.astype(np.float64),
                                    minlength=stop - start)
    ng_re[start:stop] = np.bincount(exp_r - start, minlength=stop - start)
check("every appended row equals the arm restricted to the merged gene axis", bad_rows == 0,
      f"{n_arm:,} cells, {n_cmp:,} nonzeros compared")
check("total nnz accounted", ipd[-1] == nnz_b + n_cmp,
      f"{int(ipd[-1]):,} = {nnz_b:,} + {n_cmp:,}")

print("\n=== E. depth columns recomputed from the written X ===")
nci = od["n_counts_intersect"].values
ngi = od["n_genes_intersect"].values
check("n_counts_intersect matches for the appended cells",
      bool(np.allclose(nci[n_bak:], nc_re.astype(np.float32))))
check("n_genes_intersect matches for the appended cells",
      bool(np.array_equal(ngi[n_bak:], ng_re.astype(np.int32))))
row_nnz = np.diff(ipd)
check("n_genes_intersect equals the row nnz of X on every cell",
      bool(np.array_equal(ngi.astype(np.int64), row_nnz)))
check("no appended cell has zero panel counts", bool((nci[n_bak:] > 0).all()),
      f"min {float(nci[n_bak:].min()):.0f}")

print("\n=== F. schema invariants on the merged object ===")
cs = od["composition_safe"].values.astype(bool)
wt = (od["enrichment"].astype(str) == "whole_tissue").values
check("composition_safe == (enrichment == whole_tissue)", bool((cs == wt).all()),
      f"{int(cs.sum()):,} composition-safe cells")
new = (od["source_key"].astype(str) == "EMTAB13526").values
check("composition_safe is False for every EMTAB13526 cell", bool(not cs[new].any()),
      f"{int(new.sum()):,} cells")
check("every EMTAB13526 cell is in a sorted fraction",
      bool(set(od.loc[new, "enrichment"].astype(str).unique())
           <= {"CD45pos", "CD235a_neg", "CD45pos_MDSC", "MDSC"}),
      str(od.loc[new, "enrichment"].astype(str).value_counts().to_dict()))
chk = od.groupby("donor_id", observed=True)["histology"].nunique()
check("no donor spans two histologies", bool((chk <= 1).all()))
chk = od.groupby("block_id", observed=True)["source_key"].nunique()
check("no block spans two source_keys", bool((chk <= 1).all()))
chk = od.groupby("library_id", observed=True)["donor_id"].nunique()
check("no library spans two donors", bool((chk <= 1).all()))
COMP = {"epithelial_malignant", "T_NK", "myeloid", "B_plasma", "fibroblast", "mural",
        "endothelial", "other", ""}
check("level1_compartment inside the allowed vocabulary",
      set(od["level1_compartment"].astype(str).unique()) <= COMP)
mast = od["level2_celltype"].astype(str) == "mast"
check("every mast cell is in compartment other",
      bool((od.loc[mast, "level1_compartment"].astype(str) == "other").all()),
      f"{int(mast.sum()):,} mast cells")
sub = od["level3_celltype_native"].astype(str).str.contains("t cell", case=False, regex=False)
notT = sub & (od["level1_compartment"].astype(str) != "T_NK")
check("the substring trap does not fire: labels containing 't cell' that are not T_NK are "
      "correctly filed", bool((od.loc[notT, "level1_compartment"].astype(str) == "other").all()),
      f"{int(sub.sum()):,} matching labels, {int(notT.sum()):,} of them not T_NK")
check("trailing-space author label survives the merge",
      int((od["level3_celltype_native"].astype(str)
           == "Cycling anti-inflammatory alveolar macrophages ").sum()) == 1929)
check("n_counts positive on every cell", bool((od["n_counts"].values > 0).all()))
u = read_elem(fd["uns"])
for k in ["fibroblast_limitation", "structural_compartment_limitation",
          "sorted_fraction_limitation", "platform_note", "immune_contrast"]:
    check(f"uns['emtab13526']['{k}'] recorded", k in u.get("emtab13526", {}))
bq = u["block_qc"]
check("uns['block_qc'] lists all 21 blocks", len(bq["blocks"]) == 21)
check("uns['block_qc'] cell counts match the object",
      bool(np.array_equal(np.asarray(bq["n_cells_by_block"]),
                          od.groupby("block_id", observed=True).size()
                          .reindex([str(b) for b in bq["blocks"]]).values)))
check("uns['block_qc'] per-compartment matrix sums to the object's cell count",
      int(np.asarray(bq["n_cells_by_block_by_compartment"]).sum()) == n_new)
flagged = od["block_qc_pass"].values.astype(bool)
fb_blocks = [str(x) for x in bq["failing_blocks"]]
check("block_qc_pass is False exactly on the failing blocks",
      bool(np.array_equal(~flagged, od["block_qc_pass"].index.isin(
          od.index[od["block_id"].astype(str).isin(fb_blocks)]))),
      f"failing {fb_blocks}, {int((~flagged).sum()):,} cells")

fb.close(); fd.close(); fa.close()

print()
if FAILS:
    print("*** STOP: " + str(len(FAILS)) + " checks FAILED: " + "; ".join(FAILS))
    print("*** the draft was NOT moved into place; stage1a/panel_object.h5ad is untouched")
    sys.exit(1)

print("all checks PASS")
log(f"moving {DRAFT} -> {TARGET}")
os.replace(DRAFT, TARGET)
log(f"in place: {TARGET} ({os.path.getsize(TARGET):,} bytes)")

import anndata as ad
A = ad.read_h5ad(TARGET, backed="r")
print("reopened with anndata:", A.shape)
print("obs columns:", len(A.obs.columns))
print("source_key:", A.obs.source_key.value_counts().to_dict())
log("step 3 done")
