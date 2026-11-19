#!/usr/bin/env python3
"""
Atlas stage 03, step gse325414_prepare: the three inputs of wave1_arm_gse325414.py, written from the
GSE325414 deposit.

Paper items: part of the atlas behind Figure 2, Figure 3, Supplementary Figures 2 to 4 and Supplementary Data 1. No random numbers are drawn.

The arm script reads
  g325414_panel.tsv.gz              the deposit Matrix Market entries of the panel genes
  g325414_rowmap.json               matrix row -> gene symbol for those genes
  wave1/gse325414_donor_meta.csv    histology, sex, age and stage per donor
which are built here from the deposit:

  rowmap   {matrix row (1-based, as in the .mtx file): panel gene symbol} for every panel gene of
           stage1a/panel_object.h5ad (var/gene_symbol) found in GSE325414_features.tsv.gz, by its
           own symbol or, failing that, by its previous HGNC symbol from ALIAS below, applied
           only where the alias is in this matrix and is not itself a panel gene. ALIAS is the
           dict of wave1_arm_gse246215.py, copied verbatim; it recovers 23 genes here, which
           gives 1,903 of the 1,906 panel genes (missing EEF1AKNMT, FLACC1, NICN1, as the arm's
           uns merge_safe note says). The deposit has no duplicate symbol.
  panel    tab-separated, no header, one line per nonzero entry of GSE325414_matrix.mtx.gz whose
           row is in the rowmap: row, column (both 1-based, as in the file), value.
  donor    one row per donor of the GEO series matrix (GSE325414_series_matrix.txt.gz), from the
           sample characteristics "donor", "histologic subtype", "Sex", "age" and
           "nsclc clinical_stage", in order of first appearance; the columns are named after the
           characteristic keys, as the arm script reads them.

FROM: raw GSE325414 folder (00_download.sh), stage1a/panel_object.h5ad
TO:   <scratch>/recon/g325414_rowmap.json, <scratch>/recon/g325414_panel.tsv.gz,
      output/RNA/meta_dataset/wave1/gse325414_donor_meta.csv

RUN: python atlas/steps/wave1_gse325414_prepare.py [donor|panel|all]   (normally through
     atlas/03_arms_wave1.py, which passes "all")
"""
import gzip
import json
import os
import sys
import time

import h5py
import numpy as np
import pandas as pd

import os as _os  # noqa: E402
import sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # atlas/
from common import ROOT, SCRNA, SCRATCH  # noqa: E402  paths from config.json

T0 = time.time()
SRC = f"{SCRNA}/GSE325414"
OBJ = f"{ROOT}/output/RNA/meta_dataset/stage1a/panel_object.h5ad"
REC = f"{SCRATCH}/recon"
OUT = f"{ROOT}/output/RNA/meta_dataset/wave1"
KEYS = ["donor", "histologic subtype", "Sex", "age", "nsclc clinical_stage"]
CHUNK = 20_000_000
EXPECT_GENES = 1903       # panel genes found in the deposit; a different deposit or panel stops here

# previous HGNC symbols, applied only where the alias is in this matrix AND is not itself a panel gene
# (verbatim from wave1_arm_gse246215.py)
ALIAS = {"CEP43": "FGFR1OP", "CIBAR2": "FAM92B", "CILK1": "ICK", "HOATZ": "C11orf88",
         "PABIR3": "FAM122C", "PEDS1": "TMEM189", "TAMALIN": "GRASP", "TLCD3A": "FAM57A",
         "SLC66A1": "PQLC2", "SLC49A4": "DIRC2", "POGLUT3": "KDELC2", "CRPPA": "ISPD",
         "ECRG4": "C2orf40", "BNIP5": "C6orf222", "FLACC1": "ALS2CR11", "H1-1": "HIST1H1A",
         "H1-2": "HIST1H1C", "H2AC20": "HIST2H2AA3", "H4-16": "HIST4H4", "H4C3": "HIST1H4C",
         "KRT86": "KRTHB6", "VARS1": "VARS", "CCN1": "CYR61", "CCN2": "CTGF",
         "PALM2AKAP2": "PALM2-AKAP2"}


def log(*a):
    """Print a message with the elapsed time."""
    print(f"[{time.time() - T0:7.1f}s]", *a, flush=True)


def donor_meta():
    """Per-donor characteristics from the series matrix, first appearance order."""
    chars = []
    with gzip.open(f"{SRC}/GSE325414_series_matrix.txt.gz", "rt") as fh:
        for line in fh:
            if line.startswith("!Sample_characteristics_ch1"):
                chars.append([x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]])
    per_sample = []
    for j in range(len(chars[0])):
        d = {}
        for row in chars:
            if ": " in row[j]:
                k, v = row[j].split(": ", 1)
                d[k] = v
        per_sample.append(d)
    S = pd.DataFrame(per_sample)
    missing = [k for k in KEYS if k not in S.columns]
    assert not missing, f"series matrix lacks {missing}"
    S = S[KEYS]
    per_donor = S.drop_duplicates().copy()
    assert per_donor.donor.is_unique, "a donor has two different characteristic sets"
    per_donor["age"] = per_donor["age"].astype(int)
    os.makedirs(OUT, exist_ok=True)
    per_donor.to_csv(f"{OUT}/gse325414_donor_meta.csv", index=False)
    log(f"donor table: {len(per_donor)} donors from {len(S)} samples -> {OUT}/gse325414_donor_meta.csv")
    log(per_donor["histologic subtype"].value_counts().to_string())


def panel_triplets():
    """Row map and panel-gene triplets from the deposit Matrix Market file."""
    with h5py.File(OBJ, "r") as h:
        panel = {x.decode() if isinstance(x, bytes) else str(x) for x in h["var/gene_symbol"][:]}
    feats = [ln.rstrip("\n").split("\t")[0] for ln in gzip.open(f"{SRC}/GSE325414_features.tsv.gz", "rt")]
    assert len(set(feats)) == len(feats), "duplicate gene symbols in the deposit features"
    row_of = {s: i + 1 for i, s in enumerate(feats)}
    rowmap, n_alias = {}, 0
    for g in sorted(panel):
        if g in row_of:
            rowmap[row_of[g]] = g
        elif ALIAS.get(g) in row_of and ALIAS[g] not in panel:
            rowmap[row_of[ALIAS[g]]] = g
            n_alias += 1
    rowmap = dict(sorted(rowmap.items()))
    log(f"features {len(feats):,}; panel genes {len(panel):,}; in this deposit {len(rowmap):,} "
        f"({n_alias} through a previous symbol)")
    assert len(rowmap) == EXPECT_GENES, f"{len(rowmap)} panel genes found, expected {EXPECT_GENES}"
    os.makedirs(REC, exist_ok=True)
    with open(f"{REC}/g325414_rowmap.json", "w") as fh:
        json.dump({str(k): v for k, v in rowmap.items()}, fh)
    keep = np.zeros(len(feats) + 1, dtype=bool)
    keep[list(rowmap)] = True
    n_in = n_out = 0
    with gzip.open(f"{SRC}/GSE325414_matrix.mtx.gz", "rt") as fh, \
            gzip.open(f"{REC}/g325414_panel.tsv.gz", "wt", compresslevel=1) as fo:
        line = fh.readline()
        while line.startswith("%"):
            line = fh.readline()
        nr, nc, nnz = (int(x) for x in line.split())
        assert (nr, nc) == (len(feats), sum(1 for _ in gzip.open(f"{SRC}/GSE325414_barcodes.tsv.gz", "rt")))
        for ch in pd.read_csv(fh, sep=" ", header=None, names=["r", "c", "v"], chunksize=CHUNK,
                              dtype={"r": np.int32, "c": np.int32, "v": np.int64}, engine="c"):
            n_in += len(ch)
            ch = ch[keep[ch.r.values]]
            n_out += len(ch)
            ch.to_csv(fo, sep="\t", header=False, index=False)
            log(f"  read {n_in:,} of {nnz:,} entries, kept {n_out:,}")
    assert n_in == nnz, f"read {n_in} entries, header says {nnz}"
    log(f"wrote {REC}/g325414_panel.tsv.gz ({n_out:,} entries) and g325414_rowmap.json")


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("donor", "all"):
        donor_meta()
    if which in ("panel", "all"):
        panel_triplets()
