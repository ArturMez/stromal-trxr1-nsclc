#!/usr/bin/env python
"""
Atlas stage 02, step verify_arm: independent verification of arm_emtab13526.h5ad (E-MTAB-13526).

Paper items: part of the atlas behind Figure 2, Figure 3, Supplementary Figures 2 to 4 and Supplementary Data 1. No random numbers are drawn.

The counts check does NOT reuse the extractor.  It re-reads four libraries with a plain
python line loop over the gzip stream and a plain dict of barcode -> column, then compares
whole rows of X gene by gene.  A shared bug would have to appear in both code paths.

Run: python stage2_emtab13526_15_verify_arm.py   (normally through atlas/02_arm_emtab13526.py)
"""
import gzip

import anndata as ad
import numpy as np
import pandas as pd

import os as _os  # noqa: E402
import sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # atlas/
from common import ROOT, SCRNA, SCRATCH, DATA, RSCRIPT  # noqa: E402,F401  paths from config.json
ST2 = f"{ROOT}/output/RNA/meta_dataset/stage2"
RAW = f"{SCRNA}/E-MTAB-13526/raw_counts"
PANEL = f"{ROOT}/output/RNA/meta_dataset/stage1a/panel_genes.csv"
SCHEMA = f"{DATA}/build-plan_obs_schema.csv"  # obs schema table, shipped in atlas/data/

CHECK_LIBS = ["P14_B1", "P12_B1", "P7_B1", "P1_T1"]   # small; P1 is the v2-whitelist pair
fails = []


def check(name, ok, detail=""):
    """Print PASS or FAIL for one check and record failures."""
    print(("PASS  " if ok else "FAIL  ") + name + ("  " + detail if detail else ""),
          flush=True)
    if not ok:
        fails.append(name)


def main():
    """Run the schema, cell id, honest-limit, count and metadata checks on the arm."""
    A = ad.read_h5ad(f"{ST2}/arm_emtab13526.h5ad")
    panel = pd.read_csv(PANEL)
    schema = pd.read_csv(SCHEMA)
    o = A.obs
    print(f"object {A.shape}, X {A.X.dtype}, nnz {A.X.nnz:,}\n")

    # ---- schema ---------------------------------------------------------------------
    req = schema.loc[schema.required == "yes", "column"].tolist()
    missing = [c for c in req if c not in o.columns and c != "cell_id"]
    check("every required schema column present", not missing, f"missing {missing}")
    opt = [c for c in schema.column if c not in o.columns and c != "cell_id"]
    print(f"      optional schema columns absent: {opt}")

    # ---- cell_id ---------------------------------------------------------------------
    ix = pd.Index(o.index.astype(str))
    check("cell_id globally unique", ix.is_unique, f"{ix.nunique():,} of {len(ix):,}")
    nfield = np.array([len(s.split("|")) for s in ix])
    check("cell_id has exactly 3 fields", bool((nfield == 3).all()))
    check("cell_id field 1 is the source_key",
          (ix.str.split("|").str[0] == "EMTAB13526").all())
    check("cell_id field 2 is library_id",
          (ix.str.split("|").str[1].values == o.library_id.astype(str).values).all())
    check("cell_id field 3 is the barcode",
          (ix.str.split("|").str[2].values == o.barcode.astype(str).values).all())

    # ---- honest limits ---------------------------------------------------------------
    check("composition_safe False for every cell", (~o.composition_safe).all(),
          f"{int(o.composition_safe.sum())} True")
    check("no library is whole_tissue", (o.enrichment.astype(str) != "whole_tissue").all())
    check("modality single_cell everywhere", (o.modality.astype(str) == "single_cell").all())
    check("celltype_source author everywhere", (o.celltype_source.astype(str) == "author").all())
    check("counts_are_estimates False", (~o.counts_are_estimates).all())
    for k in ["fibroblast_limitation", "structural_compartment_limitation",
              "sorted_fraction_limitation", "immune_contrast", "platform_note"]:
        check(f"uns.{k} is present and non-trivial",
              k in A.uns and len(str(A.uns[k])) > 200)
    check("no 'null' used as a label anywhere in obs",
          not any((o[c].astype(str) == "null").any() for c in o.columns
                  if o[c].dtype == object or isinstance(o[c].dtype, pd.CategoricalDtype)))

    # ---- gene axis --------------------------------------------------------------------
    check("var is the panel in panel_genes.csv order",
          list(A.var_names) == panel.gene_symbol.tolist())
    absent = A.var_names[~A.var.present_in_source.values].tolist()
    Xc = A.X.tocsc()
    colnnz = np.diff(Xc.indptr)
    check("absent genes are exact structural zeros",
          all(colnnz[A.var_names.get_loc(g)] == 0 for g in absent), f"{absent}")
    check("2056 of 2059 panel genes resolved",
          int(A.var.present_in_source.sum()) == 2056)
    check("X is non-negative integer",
          A.X.dtype == np.int32 and (A.X.data > 0).all())

    # ---- depth ----------------------------------------------------------------------
    check("n_counts >= n_counts_panel for every cell",
          (o.n_counts.values >= o.n_counts_panel.values).all())
    check("n_genes >= n_genes_panel for every cell",
          (o.n_genes.values >= o.n_genes_panel.values).all())
    check("n_counts_panel equals the row sums of X",
          (np.asarray(A.X.sum(1)).ravel() == o.n_counts_panel.values).all())
    check("author_total_counts never exceeds n_counts",
          (o.author_total_counts.values <= o.n_counts.values).all())

    # ---- compartment map ---------------------------------------------------------------
    cmap = A.uns["compartment_map"]
    m1 = dict(zip(cmap.author_cell_type, cmap.level1_compartment))
    lab = o.level3_celltype_native.astype(str)
    check("every native label is in the explicit map", set(lab) <= set(m1),
          f"{len(set(lab))} labels")
    check("level1_compartment equals the explicit map for every cell",
          (lab.map(m1).values == o.level1_compartment.astype(str).values).all())
    trail = [x for x in set(lab) if x != x.strip()]
    check("the trailing-space author label is preserved verbatim",
          trail == ["Cycling anti-inflammatory alveolar macrophages "],
          f"{trail}, {int((lab == 'Cycling anti-inflammatory alveolar macrophages ').sum())} cells")
    # the substring trap: 'malignant cell' and 'mast cell' both contain 't cell'
    tsub = lab.str.contains("t cell", case=False, regex=False)
    bad = o.loc[tsub & (o.level1_compartment.astype(str) != "T_NK"),
                "level3_celltype_native"].astype(str).value_counts()
    print(f"      labels containing the substring 't cell' that are NOT T_NK: {bad.to_dict()}")

    # ---- per-library metadata against the crosswalk ------------------------------------
    cw = pd.read_csv(f"{ST2}/crosswalk.csv").set_index("library_id")
    lt = o.groupby("library_id", observed=True).agg(
        cells=("barcode", "size"),
        histology=("histology", lambda s: s.astype(str).unique()[0]),
        tissue_class=("tissue_class", lambda s: s.astype(str).unique()[0]),
        enrichment=("enrichment", lambda s: s.astype(str).unique()[0]))
    check("cells per library equal the crosswalk author cell counts",
          (lt.cells.values == cw.loc[lt.index, "author_n_cells"].values).all())
    for c in ["histology", "tissue_class", "enrichment"]:
        check(f"{c} per library equals the SDRF-derived crosswalk value",
              (lt[c].values == cw.loc[lt.index, c].values).all())
    check("Patient 4 carries histology_confidence disputed",
          (o.loc[o.donor_id_native.astype(str) == "Patient 4",
                 "histology_confidence"].astype(str) == "disputed").all())
    check("only Patient 4 is disputed",
          set(o.loc[o.histology_confidence.astype(str) == "disputed",
                    "donor_id_native"].astype(str)) == {"Patient 4"})
    p1 = o.platform.astype(str) == "10x 3' v2"
    check("only the two P1 libraries are 10x 3' v2",
          set(o.loc[p1, "library_id"].astype(str)) == {"P1_B1", "P1_T1"},
          f"{int(p1.sum())} cells")
    check("block_id splits on platform",
          set(o.block_id.astype(str)) ==
          {"DeZuani_Cvejic_2024_10x 3' v3", "DeZuani_Cvejic_2024_10x 3' v2"})

    # ---- counts, re-read with an independent code path -----------------------------------
    print("\nindependent re-extraction from the Matrix Market files:")
    panel_row = {}          # feature row (1-based) -> panel column
    feat = pd.read_csv(f"{RAW}/P10_B1-features.tsv.gz", sep="\t", header=None,
                       names=["ensembl", "symbol", "type"])
    ens_to_row = {e: i + 1 for i, e in enumerate(feat.ensembl)}
    for j, e in enumerate(panel.ensembl_id):
        if e in ens_to_row:
            panel_row[ens_to_row[e]] = j

    for lib in CHECK_LIBS:
        sel = o.library_id.astype(str) == lib
        sub = o[sel]
        with gzip.open(f"{RAW}/{lib}-barcodes.tsv.gz", "rt") as fh:
            bc_col = {b.strip(): i + 1 for i, b in enumerate(fh)}
        want = {bc_col[b]: k for k, b in enumerate(sub.barcode.astype(str))}
        got = np.zeros((len(sub), len(panel)), dtype=np.int64)
        tot = np.zeros(len(sub), dtype=np.int64)
        ngn = np.zeros(len(sub), dtype=np.int64)
        with gzip.open(f"{RAW}/{lib}-matrix.mtx.gz", "rt") as fh:
            fh.readline()
            line = fh.readline()
            while line.startswith("%"):
                line = fh.readline()
            for line in fh:
                g, b, v = line.split()
                b = int(b)
                k = want.get(b)
                if k is None:
                    continue
                v = int(v)
                tot[k] += v
                ngn[k] += 1
                j = panel_row.get(int(g))
                if j is not None:
                    got[k, j] = v
        arm = np.asarray(A[sel].X.todense())
        check(f"{lib}: X matches an independent re-read exactly, {len(sub)} cells x "
              f"{len(panel)} genes", bool((arm == got).all()),
              f"max abs diff {int(np.abs(arm - got).max())}")
        check(f"{lib}: n_counts matches the independent full-transcriptome total",
              bool((sub.n_counts.values == tot).all()))
        check(f"{lib}: n_genes matches the independent feature count",
              bool((sub.n_genes.values == ngn).all()))

    print("\n" + ("ALL CHECKS PASSED" if not fails else f"FAILURES: {fails}"))
    return 1 if fails else 0


if __name__ == "__main__":
    raise SystemExit(main())
