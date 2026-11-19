#!/usr/bin/env python3
"""
Atlas stage 04, step state: recompute the state of the full-gene build from disk; exit non-zero on drift.

Paper items: none directly; it checks that fullgene_object.h5ad (the atlas behind Figure 2, Figure 3,
Supplementary Figures 2 to 4 and Supplementary Data 1) agrees with the arm registry
atlas/data/arm_registry.tsv. The state is recomputed from the files, not restated from text. It reads
nothing but the objects and the registry, and it fails loudly if they disagree.
Random numbers: check 6 samples 400 cells with numpy default_rng(1).

WHAT IT CHECKS
  1  Every arm the registry marks `included` exists on disk with the recorded cell count.
  2  Every arm the registry marks `EXCLUDED` is ABSENT from the merged object, so an exclusion cannot
     be quietly undone by someone re-running a build script.
  3  The merged object's cell total equals the sum of the included arms.
  4  The presence mask covers exactly the included arms.
  5  block_qc_pass and block_qc_assessed exist, and the failing blocks are named.
  6  The counts_are_estimates flag still matches integrality on a sample.

RUN: python3 atlas/steps/fullbuild_state.py   (normally through atlas/04_build_objects.py)
     exit 0 = state matches the registry; exit 1 = drift, and the drift is printed.
"""
import os
import sys
import numpy as np
import pandas as pd
import h5py

import os as _os  # noqa: E402
import sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # atlas/
from common import ROOT, SCRNA, SCRATCH, DATA, RSCRIPT  # noqa: E402,F401  paths from config.json
FB = f"{ROOT}/output/RNA/meta_dataset/fullbuild"
REG = f"{DATA}/arm_registry.tsv"  # table shipped in atlas/data/
FULL = f"{FB}/fullgene_object.h5ad"
ARMFILE = {"LuCA": "arm_luca_full.h5ad", "E-MTAB-13526": "arm_emtab13526_full.h5ad",
           "E-MTAB-6653": "arm_emtab6653_full.h5ad", "E-MTAB-6149": "arm_emtab6149_full.h5ad",
           "STTT": "arm_sttt_full.h5ad", "GSE325414": "arm_gse325414_full.h5ad",
           "GSE246215": "arm_gse246215_full.h5ad"}
drift = []


def dec(a):
    """Decode an array of bytes or strings to an array of str."""
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in a])


def col(h, p):
    """Read one obs column of an h5ad file (categorical, nullable or plain) as an array."""
    g = h[p]
    if isinstance(g, h5py.Group):
        if "categories" in g:
            return dec(g["categories"][:])[g["codes"][:]]
        if "values" in g:
            return g["values"][:]
    v = g[:]
    return dec(v) if v.dtype.kind in "SO" else v


def say(ok, msg):
    """Print one check line and record a drift item when the check fails."""
    print(f"  [{'ok  ' if ok else 'DRIFT'}] {msg}")
    if not ok:
        drift.append(msg)


R = pd.read_csv(REG, sep="\t")
inc = R[R.status == "included"]
exc = R[R.status == "EXCLUDED"]
print(f"registry: {len(inc)} included, {len(exc)} excluded\n")

print("1. included arms present with the recorded cell count")
total_arm = 0
for _, r in inc.iterrows():
    f = f"{FB}/{ARMFILE.get(r.arm_id, '')}"
    if not os.path.exists(f):
        say(False, f"{r.arm_id}: file missing")
        continue
    with h5py.File(f, "r") as h:
        n = int(h["X"].attrs["shape"][0])
    total_arm += n
    say(n == int(r.cells_full), f"{r.arm_id}: {n:,} cells, registry says {int(r.cells_full):,}")

print("\n2. excluded arms absent from the merged object")
if not os.path.exists(FULL):
    say(False, "fullgene_object.h5ad missing")
    sys.exit(1)
with h5py.File(FULL, "r") as h:
    sk = pd.Series(col(h, "obs/source_key")).astype(str)
    shape = tuple(int(x) for x in h["X"].attrs["shape"])
    arms_in_mask = sorted(k for k in h["varm/gene_present"].keys() if k != "gene_symbol")
    obs_keys = set(h["obs"].keys())
    q = col(h, "obs/block_qc_pass").astype(bool) if "block_qc_pass" in obs_keys else None
    asd = col(h, "obs/block_qc_assessed").astype(bool) if "block_qc_assessed" in obs_keys else None
    blk = pd.Series(col(h, "obs/block_id")).astype(str)
    est = pd.Series(col(h, "obs/counts_are_estimates")).astype(str).isin(["True", "1", "true"]).values
    ip = h["X/indptr"][:]
present_keys = set(sk.unique())
for _, r in exc.iterrows():
    say(r.arm_id not in present_keys, f"{r.arm_id} excluded: {r.reason[:70]}")

print("\n3. merged total equals the sum of included arms")
say(shape[0] == total_arm, f"merged {shape[0]:,} cells against arm sum {total_arm:,}")

print("\n4. presence mask covers exactly the included arms")
say(arms_in_mask == sorted(inc.arm_id.tolist()),
    f"mask arms {arms_in_mask}")

print("\n5. block QC")
say(q is not None and asd is not None, "block_qc_pass and block_qc_assessed present")
if q is not None:
    fail_blocks = sorted(pd.DataFrame({"b": blk, "q": q}).groupby("b")["q"].mean().pipe(
        lambda s: s[s < 1]).index.tolist())
    say(bool(asd.all()), f"all cells assessed: {bool(asd.all())}")
    print(f"         failing blocks: {fail_blocks}")
    print(f"         usable {int(q.sum()):,} of {len(q):,} cells")

print("\n6. estimates flag still matches integrality")
rng = np.random.default_rng(1)
pick = rng.choice(shape[0], 400, replace=False)
dis = 0
with h5py.File(FULL, "r") as h:
    for r in np.sort(pick):
        a, b = int(ip[r]), int(ip[r + 1])
        if b <= a:
            continue
        d = h["X/data"][a:b]
        if (not np.allclose(d, np.round(d))) != bool(est[r]):
            dis += 1
say(dis == 0, f"{dis} of 400 sampled cells disagree with counts_are_estimates")

# No check of malignancy attribution is made: the full-gene object holds no epithelial cell.

print("\n" + "=" * 88)
if drift:
    print(f"DRIFT: {len(drift)} item(s). The registry and the objects disagree.")
    for d in drift:
        print(f"  - {d}")
    sys.exit(1)
print("state matches the registry")
sys.exit(0)
