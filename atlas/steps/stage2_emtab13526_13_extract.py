#!/usr/bin/env python
"""
Atlas stage 02, step extract: extract panel-gene counts for the author-called barcodes
from the LOCAL Matrix Market triplets of E-MTAB-13526.

Paper items: part of the atlas behind Figure 2, Figure 3, Supplementary Figures 2 to 4 and Supplementary Data 1. No random numbers are drawn.

No remote payload is touched.  Nothing under original_data is modified.

For each of the 81 libraries:
  * the deposited -barcodes.tsv.gz IS the CellRanger whitelist and its line order IS
    the MTX column index, so barcode -> column is read from that library's own file,
    never assumed to be shared between libraries;
  * barcodes are packed 2-bits-per-base into uint32 so a 6.79M-barcode whitelist can be
    searched without a 700 MB python dict;
  * every author barcode must resolve to a whitelist position, asserted per library;
  * the MTX is streamed once.  Full-transcriptome n_counts / n_genes over all 33,538
    features are accumulated BEFORE the panel-gene filter, so the depth control is a
    real library depth and not a panel depth;
  * the panel submatrix is written per library as a CSR npz in author-obs row order.

Output: <scratch>/emtab13526_extract/<library>.npz  plus  _index.parquet
Run:    python stage2_emtab13526_13_extract.py <scratch_dir> [n_workers]
        (normally through atlas/02_arm_emtab13526.py, with 4 workers)
"""
import os
import sys
import gzip
import time

import numpy as np
import pandas as pd
import scipy.sparse as sp
from concurrent.futures import ProcessPoolExecutor

import os as _os  # noqa: E402
import sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # atlas/
from common import ROOT, SCRNA, SCRATCH, DATA, RSCRIPT  # noqa: E402,F401  paths from config.json
RAW = f"{SCRNA}/E-MTAB-13526/raw_counts"
ST2 = f"{ROOT}/output/RNA/meta_dataset/stage2"
PANEL = f"{ROOT}/output/RNA/meta_dataset/stage1a/panel_genes.csv"

SCRATCH = sys.argv[1] if len(sys.argv) > 1 else f"{SCRATCH}/e13526"
NWORK = int(sys.argv[2]) if len(sys.argv) > 2 else 4
OUTDIR = f"{SCRATCH}/emtab13526_extract"

# 2 bits per base.  N never occurs in a CellRanger whitelist; asserted below.
_LUT = np.full(256, 255, dtype=np.uint8)
for b, v in zip(b"ACGT", range(4)):
    _LUT[b] = v


def encode16(byte_matrix):
    """(n, 16) uint8 ASCII -> (n,) uint32, 2 bits per base."""
    code = _LUT[byte_matrix]
    assert code.max() < 4, "barcode contains a base outside ACGT"
    out = np.zeros(code.shape[0], dtype=np.uint32)
    for j in range(16):
        out = (out << np.uint32(2)) | code[:, j].astype(np.uint32)
    return out


def read_whitelist(lib):
    """Decompressed -barcodes.tsv.gz -> encoded uint32 array, index = MTX column - 1."""
    with gzip.open(f"{RAW}/{lib}-barcodes.tsv.gz", "rb") as fh:
        buf = fh.read()
    # every line is '<16bp>-1\n' = 19 bytes
    assert len(buf) % 19 == 0, f"{lib}: barcode file is not a clean 19-byte-per-line file"
    m = np.frombuffer(buf, dtype=np.uint8).reshape(-1, 19)
    assert (m[:, 16] == ord("-")).all() and (m[:, 17] == ord("1")).all() \
        and (m[:, 18] == ord("\n")).all(), f"{lib}: unexpected barcode line format"
    return encode16(m[:, :16])


def one(job):
    """Extract one library: full-transcriptome depth and the panel-gene CSR for its author barcodes."""
    lib, want_bc, gmap, npanel = job
    t0 = time.time()

    wl = read_whitelist(lib)
    nbc = len(wl)

    # author barcodes -> whitelist position, via a sorted search on the packed codes
    q = encode16(np.frombuffer("".join(b[:16] for b in want_bc).encode(),
                               dtype=np.uint8).reshape(-1, 16))
    order = np.argsort(wl, kind="stable")
    pos = np.searchsorted(wl[order], q)
    assert (pos < nbc).all(), f"{lib}: an author barcode sorts past the whitelist"
    col0 = order[pos]                          # 0-based whitelist row
    assert (wl[col0] == q).all(), f"{lib}: author barcode absent from the whitelist"
    assert len(np.unique(col0)) == len(col0), f"{lib}: two author cells share a barcode"

    cmap = np.full(nbc + 1, -1, dtype=np.int32)   # MTX column (1-based) -> output row
    cmap[col0 + 1] = np.arange(len(col0), dtype=np.int32)

    ncell = len(col0)
    n_counts = np.zeros(ncell, dtype=np.int64)
    n_genes = np.zeros(ncell, dtype=np.int32)
    rows, cols, vals = [], [], []
    nent_seen = 0

    with gzip.open(f"{RAW}/{lib}-matrix.mtx.gz", "rb") as fh:
        line = fh.readline().decode()
        assert line.startswith("%%MatrixMarket"), f"{lib}: not a MatrixMarket file"
        while True:
            line = fh.readline().decode()
            if not line.startswith("%"):
                break
        nfeat, nbc_hdr, nent = (int(x) for x in line.split())
        assert nfeat == 33538, f"{lib}: {nfeat} features, expected 33538"
        assert nbc_hdr == nbc, f"{lib}: header says {nbc_hdr} columns, barcode file has {nbc}"

        for ch in pd.read_csv(fh, sep=" ", header=None, dtype=np.int64,
                              chunksize=8_000_000, names=["g", "b", "v"], engine="c"):
            g = ch["g"].to_numpy()
            b = ch["b"].to_numpy()
            v = ch["v"].to_numpy()
            nent_seen += len(g)

            r = cmap[b]
            keep = r >= 0
            if not keep.any():
                continue
            r = r[keep]
            gk = g[keep]
            vk = v[keep]
            # full-transcriptome depth, before any gene subsetting
            n_counts += np.bincount(r, weights=vk, minlength=ncell).astype(np.int64)
            n_genes += np.bincount(r, minlength=ncell).astype(np.int32)

            c = gmap[gk]
            keep2 = c >= 0
            if keep2.any():
                rows.append(r[keep2].astype(np.int32))
                cols.append(c[keep2].astype(np.int32))
                vals.append(vk[keep2].astype(np.int32))

    assert nent_seen == nent, f"{lib}: read {nent_seen} entries, header declares {nent}"
    assert (n_counts > 0).all(), f"{lib}: an author cell has zero counts in the raw matrix"

    rows = np.concatenate(rows) if rows else np.zeros(0, np.int32)
    cols = np.concatenate(cols) if cols else np.zeros(0, np.int32)
    vals = np.concatenate(vals) if vals else np.zeros(0, np.int32)
    X = sp.coo_matrix((vals, (rows, cols)), shape=(ncell, npanel), dtype=np.int32).tocsr()
    X.sum_duplicates()
    assert X.nnz == len(vals), f"{lib}: duplicate (cell, gene) entries in the MTX"

    np.savez(f"{OUTDIR}/{lib}.npz", data=X.data, indices=X.indices, indptr=X.indptr,
             shape=np.array(X.shape), n_counts=n_counts, n_genes=n_genes,
             col0=col0.astype(np.int64), whitelist_size=np.array([nbc]))
    return dict(library=lib, cells=ncell, whitelist=nbc, mtx_entries=nent,
                panel_entries=int(X.nnz), panel_umi=int(X.data.sum()),
                total_umi_full=int(n_counts.sum()), seconds=round(time.time() - t0, 1))


def main():
    """Run the extraction over all 81 libraries with a process pool and write _index.parquet."""
    os.makedirs(OUTDIR, exist_ok=True)

    panel = pd.read_csv(PANEL)
    feat = pd.read_csv(f"{RAW}/P10_B1-features.tsv.gz", sep="\t", header=None,
                       names=["ensembl", "symbol", "type"])
    assert feat.ensembl.is_unique

    # ---- gene resolution: Ensembl id is the join key ---------------------------------
    row_of = {e: i for i, e in enumerate(feat.ensembl)}      # 0-based feature row
    gmap = np.full(33539, -1, dtype=np.int32)                # MTX gene index (1-based)
    present = np.zeros(len(panel), dtype=bool)
    matrix_row = np.full(len(panel), -1, dtype=np.int64)
    for i, e in enumerate(panel.ensembl_id):
        r = row_of.get(e)
        if r is not None:
            gmap[r + 1] = i
            present[i] = True
            matrix_row[i] = r
    print(f"panel genes {len(panel)}; resolved by Ensembl id {present.sum()}; "
          f"absent {int((~present).sum())}: {panel.gene_symbol[~present].tolist()}")

    res = pd.DataFrame({
        "gene_symbol": panel.gene_symbol, "ensembl_id": panel.ensembl_id,
        "present_in_source": present, "matrix_feature_row_0based": matrix_row,
        "source_symbol": [feat.symbol.iloc[r] if r >= 0 else "" for r in matrix_row],
    })
    res["symbol_agrees"] = res.present_in_source & (res.source_symbol == res.gene_symbol)
    res.to_csv(f"{ST2}/arm_emtab13526_gene_resolution.csv", index=False)
    drift = res[res.present_in_source & ~res.symbol_agrees]
    print(f"symbol drift between panel and the CellRanger GRCh38-3.0.0 reference: {len(drift)} genes")

    # ---- author cells per library -----------------------------------------------------
    ob = pd.read_parquet(f"{ST2}/author_obs.parquet",
                         columns=["author_object", "author_barcode", "patient_sample"])
    cw = pd.read_csv(f"{ST2}/crosswalk.csv")
    assert (cw.library_id == cw.library_id_F).all(), "crosswalk library_id != fingerprint route"
    lib_of = dict(zip(cw.patient_sample, cw.library_id))
    ob["library_id"] = ob.patient_sample.map(lib_of)
    assert ob.library_id.notna().all(), "an author cell has no library"
    ob["barcode"] = ob.author_barcode.str.split("-").str[:2].str.join("-")
    assert (ob.barcode.str.len() == 18).all()
    ob = ob.reset_index(drop=True)
    ob["row_in_arm"] = -1

    # deterministic order: library sorted, author order preserved inside a library
    order = np.argsort(ob.library_id.values, kind="stable")
    ob = ob.iloc[order].reset_index(drop=True)
    ob["row_in_arm"] = np.arange(len(ob))
    ob[["library_id", "barcode", "author_object", "author_barcode", "row_in_arm"]] \
        .to_parquet(f"{OUTDIR}/_index.parquet", index=False)

    jobs = []
    for lib, g in ob.groupby("library_id", sort=True):
        jobs.append((lib, g.barcode.tolist(), gmap, len(panel)))
    # biggest first so the tail does not serialise
    size = {f.split("-matrix")[0]: os.path.getsize(f"{RAW}/{f}")
            for f in os.listdir(RAW) if f.endswith("-matrix.mtx.gz")}
    jobs.sort(key=lambda j: -size[j[0]])
    print(f"{len(jobs)} libraries, {len(ob):,} author cells, {NWORK} workers")

    t0 = time.time()
    recs = []
    with ProcessPoolExecutor(max_workers=NWORK) as ex:
        for r in ex.map(one, jobs):
            recs.append(r)
            print(f"  {r['library']:8s} cells {r['cells']:7,} "
                  f"panel_nnz {r['panel_entries']:11,} {r['seconds']:6.1f}s", flush=True)
    wall = time.time() - t0
    df = pd.DataFrame(recs).sort_values("library")
    df.to_csv(f"{ST2}/arm_emtab13526_extract_log.csv", index=False)
    print(f"\nwall {wall:.0f}s, cpu {df.seconds.sum():.0f}s")
    print(f"cells {df.cells.sum():,}  panel entries {df.panel_entries.sum():,}  "
          f"panel UMI {df.panel_umi.sum():,}  full UMI {df.total_umi_full.sum():,}")
    print(f"whitelist sizes: {df.whitelist.value_counts().to_dict()}")


if __name__ == "__main__":
    main()
