#!/usr/bin/env python
"""
Atlas stage 02, step fingerprint: per-barcode UMI fingerprint of all 81 raw libraries.

Paper items: part of the atlas behind Figure 2, Figure 3, Supplementary Figures 2 to 4 and Supplementary Data 1.
Run: python atlas/steps/stage2_emtab13526_06_fingerprint.py   (normally through atlas/02_arm_emtab13526.py).
No random numbers are drawn.

The author obs carries total_counts as an exact integer UMI total per cell.
That makes (barcode -> UMI total) a fingerprint that identifies which raw
library each author sample group came from, with no reliance on any naming
convention, on the author's sangerID free text, or on the SDRF's own site
labels. This is what settles Patient 4 and Donor 2.

Streams each -matrix.mtx.gz once (read-only) and writes, per library, the
barcodes clearing a 300-UMI floor with their UMI total and detected-gene count.
The author minimum total_counts is 400, so a 300 floor cannot lose an author cell.

Output: output/RNA/meta_dataset/stage2/fingerprint/<library>.parquet
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
FP = os.path.join(ROOT, "output/RNA/meta_dataset/stage2/fingerprint")
FLOOR = 300

os.makedirs(FP, exist_ok=True)


def one(lib):
    """Stream one library's matrix and write its barcodes above the 300-UMI floor with UMI total and gene count."""
    out = os.path.join(FP, f"{lib}.parquet")
    if os.path.exists(out):
        return {"library": lib, "status": "cached",
                "n_kept": len(pd.read_parquet(out)), "seconds": 0.0}
    t0 = time.time()
    path = os.path.join(RAW, f"{lib}-matrix.mtx.gz")
    with gzip.open(path, "rb") as fh:
        hdr = []
        while len(hdr) < 3:
            hdr.append(fh.readline().decode())
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

    keep = np.nonzero(tot >= FLOOR)[0]          # 1-based MM column indices
    with gzip.open(os.path.join(RAW, f"{lib}-barcodes.tsv.gz"), "rt") as fh:
        barcodes = np.array([l.rstrip("\n") for l in fh], dtype=object)
    if len(barcodes) != nbc:
        raise RuntimeError(f"{lib}: {len(barcodes)} barcodes vs {nbc} MM columns")

    df = pd.DataFrame({
        "barcode": barcodes[keep - 1],
        "raw_umi": tot[keep].astype(np.int64),
        "raw_ngenes": ngene[keep].astype(np.int32),
    })
    df.to_parquet(out, index=False)
    return {"library": lib, "status": "built", "n_kept": len(df),
            "seconds": round(time.time() - t0, 1)}


if __name__ == "__main__":
    libs = sorted({f.split("-matrix")[0] for f in os.listdir(RAW)
                   if f.endswith("-matrix.mtx.gz")})
    print(f"{len(libs)} libraries, floor {FLOOR} UMI", flush=True)
    # biggest first so the long poles start early; cap workers for memory
    sz = {l: os.path.getsize(os.path.join(RAW, f"{l}-matrix.mtx.gz")) for l in libs}
    libs = sorted(libs, key=lambda l: -sz[l])
    t0 = time.time()
    recs = []
    with ProcessPoolExecutor(max_workers=3) as ex:
        for r in ex.map(one, libs):
            recs.append(r)
            print(f"  {r['library']:8s} {r['status']:7s} kept={r['n_kept']:7d} "
                  f"{r['seconds']:6.1f}s", flush=True)
    df = pd.DataFrame(recs)
    print(f"\nwall {time.time()-t0:.1f}s   total barcodes kept {df['n_kept'].sum():,}")
