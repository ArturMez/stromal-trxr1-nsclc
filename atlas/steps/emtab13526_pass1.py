#!/usr/bin/env python
"""
E-MTAB-13526 scoping pass 1: per-barcode UMI totals only.

This is NOT cell calling. It streams each Matrix Market triplet file once and
records how many barcodes clear a set of UMI floors, plus the wall time and
peak RSS. Those numbers are what a cell-calling budget has to be built on.

Output: output/RNA/meta_dataset/emtab13526_pass1.csv
Input: the 81 Matrix Market triplet files of E-MTAB-13526 (raw_counts/*-matrix.mtx.gz).
Floors: 1, 100, 200, 500, 1000 and 2000 UMI per barcode (counts of barcodes at or above each).
Run by atlas/02_arm_emtab13526.py (step pass1), 4 worker processes; no random numbers are drawn.
Part of the E-MTAB-13526 arm of the atlas behind Figure 2, Figure 3, Supplementary Figures 2 to 4.
"""

import gzip
import os
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

import os as _os  # noqa: E402
import sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # atlas/
from common import ROOT, SCRNA, SCRATCH, DATA, RSCRIPT  # noqa: E402,F401  paths from config.json
RAW = os.path.join(SCRNA, "E-MTAB-13526/raw_counts")
OUT = os.path.join(ROOT, "output/RNA/meta_dataset")

FLOORS = [1, 100, 200, 500, 1000, 2000]


def one(lib):
    """Stream one library's triplet file and return its per-floor barcode counts and UMI summary."""
    path = os.path.join(RAW, f"{lib}-matrix.mtx.gz")
    t0 = time.time()
    nbc = None
    tot = None
    nrows = 0
    with gzip.open(path, "rb") as fh:
        # header: %%MatrixMarket / %comment / dims
        hdr = []
        while len(hdr) < 3:
            line = fh.readline().decode()
            hdr.append(line)
        nfeat, nbc, nent = (int(x) for x in hdr[-1].split())
        tot = np.zeros(nbc + 1, dtype=np.int64)
        ngene = np.zeros(nbc + 1, dtype=np.int32)
        for ch in pd.read_csv(fh, sep=" ", header=None, dtype=np.int64,
                              chunksize=8_000_000, names=["g", "b", "v"],
                              engine="c"):
            b = ch["b"].to_numpy()
            v = ch["v"].to_numpy()
            tot += np.bincount(b, weights=v, minlength=nbc + 1).astype(np.int64)
            ngene += np.bincount(b, minlength=nbc + 1).astype(np.int32)
            nrows += len(ch)
    dt = time.time() - t0
    rec = {"library": lib, "n_features": nfeat, "n_barcodes_whitelist": nbc,
           "nonzero_entries": nent, "rows_read": nrows,
           "total_umi": int(tot.sum()), "pass1_seconds": round(dt, 2)}
    for f in FLOORS:
        rec[f"bc_ge_{f}umi"] = int((tot >= f).sum())
    nz = tot[tot > 0]
    rec["max_umi_barcode"] = int(tot.max())
    rec["median_umi_nonempty"] = float(np.median(nz)) if nz.size else 0.0
    # genes detected in the barcodes that clear 500 UMI, a crude cell-quality read
    sel = tot >= 500
    rec["median_genes_bc_ge500"] = float(np.median(ngene[sel])) if sel.sum() else 0.0
    return rec


if __name__ == "__main__":
    libs = sorted({f.split("-matrix")[0] for f in os.listdir(RAW)
                   if f.endswith("-matrix.mtx.gz")})
    print(f"{len(libs)} libraries")
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=4) as ex:
        recs = list(ex.map(one, libs))
    wall = time.time() - t0
    df = pd.DataFrame(recs).sort_values("library")
    df.to_csv(os.path.join(OUT, "emtab13526_pass1.csv"), index=False)
    print(f"wall {wall:.1f}s for {len(libs)} libraries with 4 workers")
    print(f"sum of per-library cpu seconds = {df['pass1_seconds'].sum():.1f}")
    print(df[["library", "nonzero_entries", "bc_ge_100umi", "bc_ge_500umi",
              "pass1_seconds"]].to_string(index=False))
