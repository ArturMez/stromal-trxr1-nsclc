#!/usr/bin/env python
"""
Atlas stage 02, step author_obs: read the author obs of E-MTAB-13526 over HTTP range reads.

Paper items: part of the atlas behind Figure 2, Figure 3, Supplementary Figures 2 to 4 and Supplementary Data 1.
Run: python atlas/steps/stage2_emtab13526_01_read_author_obs.py   (normally through atlas/02_arm_emtab13526.py).
Needs network access; no random numbers are drawn.

Reads the complete `obs` table of both author-annotated h5ad files over HTTP
RANGE READS. The payloads are 97 GiB combined and are NOT downloaded; only the
HDF5 metadata pages and the /obs datasets are fetched. Fetched volume is
reported at the end so the claim is checkable.

Two traps handled here:
  1. The files use the older anndata 0.7 categorical encoding. A categorical
     column is an int8 dataset of codes carrying attrs['categories'], an HDF5
     OBJECT REFERENCE into /obs/__categories/<name>. A naive reader returns the
     raw integer codes silently. The reference is dereferenced explicitly and
     the dereferenced target name is checked against the column name.
  2. fsspec's HTTPFile can stall without a read timeout. Reads go through
     rangefile.py (RangeFile), which sets timeouts and retries.

Each object is checkpointed to its own parquet, so a failure on the second file
does not throw away the first.

Writes:
  output/RNA/meta_dataset/stage2/author_obs_TUMOR.parquet
  output/RNA/meta_dataset/stage2/author_obs_BACKGROUND.parquet
  output/RNA/meta_dataset/stage2/author_obs.parquet      (stacked)
  output/RNA/meta_dataset/stage2/author_obs_columns.csv  (column inventory)
"""

import os
import sys
import time

import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rangefile import RangeFile  # noqa: E402

import os as _os  # noqa: E402
import sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # atlas/
from common import ROOT, SCRNA, SCRATCH, DATA, RSCRIPT  # noqa: E402,F401  paths from config.json
OUT = os.path.join(ROOT, "output/RNA/meta_dataset/stage2")
os.makedirs(OUT, exist_ok=True)

BASE_URL = "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/526/E-MTAB-13526/Files"
FILES = {
    "TUMOUR": "10X_Lung_Tumour_Annotated_v2.h5ad",
    "BACKGROUND": "10X_Lung_Healthy_Background_Annotated_v2.h5ad",
}


def decode(arr):
    """Decode a bytes or object array to str (other dtypes pass through)."""
    arr = np.asarray(arr)
    if arr.dtype.kind in ("S", "O"):
        return np.array([x.decode("utf-8") if isinstance(x, (bytes, np.bytes_)) else x
                         for x in arr], dtype=object)
    return arr


def read_obs_group(h5, g):
    """Read every obs column of an old-encoding h5ad obs group into a data frame, dereferencing the category references."""
    enc = dict(g.attrs)
    col_order = [c.decode() if isinstance(c, bytes) else c
                 for c in np.asarray(enc.get("column-order", [])).ravel()]
    idx_key = enc.get("_index", "_index")
    if isinstance(idx_key, bytes):
        idx_key = idx_key.decode()

    keys = [k for k in g.keys() if k != "__categories"]
    ordered = [idx_key] + [c for c in col_order if c in keys] + \
              [k for k in keys if k != idx_key and k not in col_order]
    seen = set()
    ordered = [k for k in ordered if k in keys and not (k in seen or seen.add(k))]

    out, mismatches = {}, []
    for k in ordered:
        node = g[k]
        if isinstance(node, h5py.Group):
            raise RuntimeError(f"unexpected group column {k!r}")
        t0 = time.time()
        if "categories" in node.attrs:
            catnode = h5[node.attrs["categories"]]      # dereference, never guess by name
            deref = catnode.name.split("/")[-1]
            if deref != k:
                mismatches.append(f"{k} -> {catnode.name}")
            cats = [str(c) for c in decode(catnode[:])]
            codes = np.asarray(node[:]).astype(np.int64)
            if codes.size and (codes.max() >= len(cats) or codes.min() < -1):
                raise RuntimeError(f"{k}: codes out of range "
                                   f"[{codes.min()},{codes.max()}] vs {len(cats)} categories")
            out[k] = pd.Categorical.from_codes(codes, categories=cats)
            kind = f"cat[{len(cats)}]"
        else:
            out[k] = decode(node[:])
            kind = str(np.asarray(out[k]).dtype)
        print(f"    {k:36s} {kind:12s} {time.time()-t0:6.1f}s", flush=True)

    for m in mismatches:
        print("  NOTE categories reference name differs from column:", m, flush=True)

    return pd.DataFrame(out).rename(columns={idx_key: "author_barcode"})


def read_one(obj_key, fn):
    """Read the obs of one author file (or load its checkpoint parquet)."""
    ck = os.path.join(OUT, f"author_obs_{obj_key}.parquet")
    if os.path.exists(ck):
        print(f"=== {obj_key}: checkpoint exists, loading {ck}", flush=True)
        return pd.read_parquet(ck)

    url = f"{BASE_URL}/{fn}"
    print(f"\n=== {obj_key}  {fn} ===", flush=True)
    t0 = time.time()
    rf = RangeFile(url, block_size=1 << 20, verbose=True)
    print(f"  size {rf.size:,} bytes ({rf.size/2**30:.2f} GiB), Range honoured", flush=True)
    with h5py.File(rf, "r") as h5:
        print("  top-level keys:", list(h5.keys()), flush=True)
        if "X" in h5:
            print("  X:", {k: (v.tolist() if hasattr(v, "tolist") else v)
                           for k, v in h5["X"].attrs.items()}, flush=True)
        if "obsm" in h5:
            print("  obsm keys:", list(h5["obsm"].keys()), flush=True)
        if "var" in h5:
            ik = h5["var"].attrs.get("_index", "_index")
            ik = ik.decode() if isinstance(ik, bytes) else ik
            print("  n_var:", h5["var"][ik].shape[0],
                  " var cols:", list(h5["var"].keys()), flush=True)
        df = read_obs_group(h5, h5["obs"])
    print(f"  obs shape {df.shape}  wall {time.time()-t0:.1f}s  net {rf.stats()}", flush=True)
    rf.close()

    df.insert(0, "author_object", obj_key)
    df.to_parquet(ck, index=False)
    print("  checkpoint written:", ck, flush=True)
    return df


def main():
    """Read both author files, stack their obs, and write the parquet and column inventory."""
    frames, colinv = [], []
    for obj_key, fn in FILES.items():
        df = read_one(obj_key, fn)
        for c in df.columns:
            s = df[c]
            nun = int(s.nunique(dropna=True))
            lv = sorted(map(str, pd.unique(s.dropna())))
            colinv.append({"author_object": obj_key, "column": c,
                           "dtype": str(s.dtype), "n_unique": nun,
                           "levels": " | ".join(lv) if nun <= 100 else ""})
        frames.append(df)

    allobs = pd.concat(frames, ignore_index=True, sort=False)
    for c in allobs.columns:
        if str(allobs[c].dtype) == "object":
            allobs[c] = allobs[c].astype("string")
    allobs.to_parquet(os.path.join(OUT, "author_obs.parquet"), index=False)
    pd.DataFrame(colinv).to_csv(os.path.join(OUT, "author_obs_columns.csv"), index=False)

    print("\n=== combined ===")
    print("rows:", len(allobs))
    print(allobs["author_object"].value_counts().to_string())
    print("wrote", os.path.join(OUT, "author_obs.parquet"))


if __name__ == "__main__":
    main()
