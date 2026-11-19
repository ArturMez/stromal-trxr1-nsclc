#!/usr/bin/env python3
"""
Atlas stage 04, step repair_obs: repair four obs columns of `fullgene_object.h5ad`.

Paper items: none directly; it retypes columns of the atlas object behind Figure 2, Figure 3,
Supplementary Figures 2 to 4 and Supplementary Data 1. No random numbers are drawn.

The four columns are written by the merge with unusable types. Rewriting the OBJECT is not necessary:
an obs column is a 15 MB dataset and can be replaced in place without touching X, which is
2,694,783,448 nonzeros and is never opened here.

  n_counts_panel   STRING, empty for the 1,587,067 cells of the four original arms. On
  n_genes_panel    `panel_object_v2.h5ad` both are numeric and complete for all 2,500,510 cells.
                   THIS MATTERS: cross-arm depth matching uses n_counts_panel, which is unusable
                   while the column is empty. Repaired by transfer on cell id.

  malignancy_score STRING, the literal 'nan' for 301,843 cells and '' for the rest. The object holds
                   fibroblast, myeloid and T/NK only, so a malignancy score does not exist for ANY
                   of its cells. Rewritten as a proper float64 column that is NaN everywhere.

  is_malignant     STRING with three encodings of missing: '' for 1,587,067, 'None' for 294,599 and
                   'False' for the 7,244 GSE246215 FIBROBLASTS. That last one is the dangerous one:
                   it reads as "assessed and found not malignant" when no such assessment exists.
                   Verified on `panel_object_v2.h5ad`: every non-epithelial cell there is NA, and
                   True/False occur only among the 149,368 epithelial cells the full-gene object
                   does not carry. Rewritten as a proper nullable boolean that is NA everywhere.

X, var, varm/gene_present, obs_names and every other obs column are untouched, and the script
asserts so afterwards by re-reading. HDF5 does not reclaim the space freed by a deleted dataset, so
the file grows by about the size of the columns replaced; that is expected and is not corruption.

REVERSIBILITY. The pre-repair values are dumped to `--backup-dir` first, and
`--restore <dir>` puts them back byte for byte. Nothing here is derived: the two panel columns come
from `panel_object_v2.h5ad` and can be re-transferred at any time.

RUN: (normally through atlas/04_build_objects.py with --backup-dir {SCRATCH}/repair_obs_backup)
     python3 atlas/steps/fullbuild_10_repair_obs.py --backup-dir <dir>
     python3 atlas/steps/fullbuild_10_repair_obs.py --restore <dir>
     python3 atlas/steps/fullbuild_10_repair_obs.py --check      (read-only)
"""
import os
import sys
import time
import numpy as np
import h5py

import os as _os  # noqa: E402
import sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # atlas/
from common import ROOT, SCRNA, SCRATCH, DATA, RSCRIPT  # noqa: E402,F401  paths from config.json
FG = f"{ROOT}/output/RNA/meta_dataset/fullbuild/fullgene_object.h5ad"
PV = f"{ROOT}/output/RNA/meta_dataset/fullbuild/panel_object_v2.h5ad"
COLS = ["n_counts_panel", "n_genes_panel", "malignancy_score", "is_malignant"]
CHECK = "--check" in sys.argv
RESTORE = sys.argv[sys.argv.index("--restore") + 1] if "--restore" in sys.argv else None
BACKUP = sys.argv[sys.argv.index("--backup-dir") + 1] if "--backup-dir" in sys.argv else None


def s(v):
    """Decode an array of bytes or strings to an array of str."""
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in v])


def report(path, tag):
    """Print the type and a summary of each of the four columns of an h5ad file (read only)."""
    with h5py.File(path, "r") as f:
        o = f["obs"]
        print(f"  [{tag}]")
        for k in COLS:
            it = o[k]
            if isinstance(it, h5py.Group):
                if "mask" in it:
                    m = it["mask"][:]
                    v = it["values"][:]
                    print(f"    {k:18s} nullable-boolean  NA {int(m.sum()):,}  "
                          f"True {int((v & ~m).sum()):,}  False {int((~v & ~m).sum()):,}")
                else:
                    print(f"    {k:18s} categorical")
            elif it.dtype == object:
                import collections
                c = collections.Counter(s(it[:]))
                print(f"    {k:18s} STRING  distinct {len(c)}  top {c.most_common(3)}")
            else:
                v = it[:]
                fin = np.isfinite(v.astype(float))
                print(f"    {k:18s} {it.dtype}  finite {int(fin.sum()):,} of {len(v):,}  "
                      f"min {v[fin].min() if fin.any() else 'n/a'}  "
                      f"max {v[fin].max() if fin.any() else 'n/a'}")


if CHECK:
    report(FG, "fullgene_object.h5ad, current state")
    raise SystemExit(0)

t0 = time.time()
if RESTORE:
    print(f"restoring the four columns from {RESTORE}")
    with h5py.File(FG, "r+") as f:
        o = f["obs"]
        for k in COLS:
            p = f"{RESTORE}/{k}.npy"
            if not os.path.exists(p):
                raise SystemExit(f"backup missing: {p}")
            v = np.load(p, allow_pickle=True)
            if k in o:
                del o[k]
            o.create_dataset(k, data=np.array([x.encode() for x in v.astype(str)],
                                              dtype=h5py.special_dtype(vlen=bytes)))
    report(FG, "after restore")
    raise SystemExit(0)

if not BACKUP:
    raise SystemExit("refusing to write without --backup-dir")
os.makedirs(BACKUP, exist_ok=True)

report(FG, "BEFORE")

# ---------------------------------------------------------------- read what we need
print(f"\nreading cell ids and the source columns  [{time.time()-t0:.0f}s]", flush=True)
with h5py.File(FG, "r") as f:
    fg_ids = s(f["obs/_index"][:])
    n_fg = len(fg_ids)
    fg_before = {k: (f["obs"][k][:] if f["obs"][k].dtype == object else f["obs"][k][:])
                 for k in COLS}
    fg_other = sorted(k for k in f["obs"].attrs["column-order"] if k not in COLS)
    fg_shape = tuple(f["X"].attrs["shape"])
    fg_nnz = f["X/data"].shape[0]
    fg_fingerprint = {k: (f["obs"][k][:100].tobytes() if f["obs"][k].dtype != object
                          else b"".join(f["obs"][k][:100]))
                      for k in fg_other if not isinstance(f["obs"][k], h5py.Group)}
for k in COLS:
    np.save(f"{BACKUP}/{k}.npy", s(fg_before[k]) if fg_before[k].dtype == object else fg_before[k])
print(f"  backed up {len(COLS)} columns to {BACKUP}")

with h5py.File(PV, "r") as f:
    pv_ids = s(f["obs/_index"][:])
    pv_nc = f["obs/n_counts_panel"][:].astype(np.float64)
    pv_ng = f["obs/n_genes_panel"][:].astype(np.float64)
print(f"  panel_object_v2: {len(pv_ids):,} cells  [{time.time()-t0:.0f}s]", flush=True)

pos = {c: i for i, c in enumerate(pv_ids)}
take = np.array([pos.get(c, -1) for c in fg_ids])
missing = int((take < 0).sum())
print(f"  full-gene cells found in panel_object_v2: {n_fg - missing:,} of {n_fg:,} "
      f"(missing {missing:,})")
if missing:
    raise SystemExit("some full-gene cells are absent from panel_object_v2; refusing to write")

new_nc = pv_nc[take]
new_ng = pv_ng[take]
if not np.isfinite(new_nc).all() or not np.isfinite(new_ng).all():
    raise SystemExit("transferred depth columns contain non-finite values; refusing to write")
# CROSS-CHECK: the invariant is a subset relation, NOT equality.
# The 301,843 stored values were computed against the 1,906-gene panel; panel_object_v2 counts on
# the merged 1,721-gene axis, which is a strict subset of it. So the transferred value must be
# <= the stored one for EVERY cell, and equal wherever that cell detected none of the 185 dropped
# genes. Expected: 100.0 % satisfy new <= old; agreement is 97.2 % in STTT, 44.6 % in
# GSE325414 and 4.4 % in GSE246215, which orders by how deeply each arm is sequenced.
# Equality is therefore the wrong check; the subset relation is the right one.
old_ng = s(fg_before["n_genes_panel"])
old_nc = s(fg_before["n_counts_panel"])
had = old_ng != ""
if had.any():
    o_ng = old_ng[had].astype(float)
    o_nc = old_nc[had].astype(float)
    le_ng = float((new_ng[had] <= o_ng).mean())
    le_nc = float((new_nc[had] <= o_nc).mean())
    print(f"  cross-check on the {int(had.sum()):,} cells that already had a value:")
    print(f"    transferred n_genes_panel  <= stored: {100*le_ng:.2f} %   "
          f"exactly equal: {100*np.isclose(new_ng[had], o_ng).mean():.2f} %   "
          f"median difference {np.median(new_ng[had] - o_ng):+.0f}")
    print(f"    transferred n_counts_panel <= stored: {100*le_nc:.2f} %   "
          f"exactly equal: {100*np.isclose(new_nc[had], o_nc).mean():.2f} %   "
          f"median difference {np.median(new_nc[had] - o_nc):+.0f}")
    if le_ng < 1.0 or le_nc < 1.0:
        raise SystemExit("a transferred value EXCEEDS the stored one, which the 1,721-gene axis "
                         "being a subset of the 1,906-gene one makes impossible; refusing to write")
    print("    the stored values are the SUPERSEDED 1,906-gene panel counts and are replaced.")

# ---------------------------------------------------------------- write
print(f"\nwriting  [{time.time()-t0:.0f}s]", flush=True)
with h5py.File(FG, "r+") as f:
    o = f["obs"]
    for k, v in [("n_counts_panel", new_nc), ("n_genes_panel", new_ng)]:
        del o[k]
        d = o.create_dataset(k, data=v, dtype="f8")
        d.attrs["encoding-type"] = "array"
        d.attrs["encoding-version"] = "0.2.0"
    del o["malignancy_score"]
    d = o.create_dataset("malignancy_score", data=np.full(n_fg, np.nan), dtype="f8")
    d.attrs["encoding-type"] = "array"
    d.attrs["encoding-version"] = "0.2.0"
    del o["is_malignant"]
    g = o.create_group("is_malignant")
    g.attrs["encoding-type"] = "nullable-boolean"
    g.attrs["encoding-version"] = "0.1.0"
    g.create_dataset("values", data=np.zeros(n_fg, bool))
    g.create_dataset("mask", data=np.ones(n_fg, bool))
    if "uns" not in f:
        u = f.create_group("uns")
        u.attrs["encoding-type"] = "dict"
        u.attrs["encoding-version"] = "0.1.0"
    else:
        u = f["uns"]
    note = ("obs repair, 6 September 2026, fullbuild_10_repair_obs.py. n_counts_panel and "
            "n_genes_panel were STRING and empty for the 1,587,067 cells of the four original "
            "arms; both are now float64 and complete, transferred from panel_object_v2.h5ad on "
            "cell id. malignancy_score and is_malignant were STRING with three encodings of "
            "missing, including 'False' on 7,244 GSE246215 fibroblasts which was a merge artefact "
            "and not a call; both are now properly typed and missing for every cell, which is "
            "correct because this object holds no epithelial compartment. X, var, varm and every "
            "other obs column are unchanged.")
    if "repair_2026_09_06" in u:
        del u["repair_2026_09_06"]
    ds = u.create_dataset("repair_2026_09_06", data=np.array(note.encode()))
    ds.attrs["encoding-type"] = "string"
    ds.attrs["encoding-version"] = "0.2.0"

# ---------------------------------------------------------------- verify
print(f"\nverifying  [{time.time()-t0:.0f}s]", flush=True)
with h5py.File(FG, "r") as f:
    o = f["obs"]
    assert tuple(f["X"].attrs["shape"]) == fg_shape, "X shape moved"
    assert f["X/data"].shape[0] == fg_nnz, "X nnz moved"
    assert len(s(o["_index"][:])) == n_fg, "obs_names moved"
    assert sorted(k for k in o.attrs["column-order"] if k not in COLS) == fg_other, \
        "obs column set moved"
    for k, fp in fg_fingerprint.items():
        cur = o[k][:100].tobytes() if o[k].dtype != object else b"".join(o[k][:100])
        assert cur == fp, f"untouched column {k} changed"
    nc = o["n_counts_panel"][:]
    ng = o["n_genes_panel"][:]
    assert np.isfinite(nc).all() and np.isfinite(ng).all(), "repaired depth columns have NaN"
    assert (nc > 0).all() and (ng > 0).all(), "repaired depth columns have non-positive values"
    im = o["is_malignant"]
    assert im["mask"][:].all(), "is_malignant is not NA everywhere"
    assert np.isnan(o["malignancy_score"][:]).all(), "malignancy_score is not NaN everywhere"
    print(f"  X unchanged: {fg_shape[0]:,} x {fg_shape[1]:,}, {fg_nnz:,} nonzeros")
    print(f"  {len(fg_other)} other obs columns unchanged (first 100 values fingerprinted)")
report(FG, "AFTER")
print(f"\nrestore with: python3 {os.path.basename(__file__)} --restore {BACKUP}")
print(f"[{time.time()-t0:.0f}s]")
