#!/usr/bin/env python
"""
Atlas stage 01, step arm_emtab6149: the E-MTAB-6149 arm (Lambrechts / Thienpont 2018, discovery
cohort, patients 1 to 5) at the panel gene axis
  ->  output/RNA/meta_dataset/stage1a/arm_emtab6149.h5ad

Paper items: part of the atlas behind Figure 2, Figure 3, Supplementary Figures 2 to 4 and
Supplementary Data 1 (the Lambrechts copies are included once, here and in E-MTAB-6653).
No random numbers are drawn. Asserts 52,698 cells.

Input is the scratch export written by stage1a_export6149.R, which inverts the author
ln(1 + c*1e4/T) normalization back to integer UMI counts and pulls the author 52-subtype
taxonomy out of the six *.Cellview.Rds lineage objects.

Nothing here is taken on trust:
  * the library -> patient -> sampling-site map is DERIVED (barcode length gives the 10x
    chemistry block, an XIST/RPS4Y1 probe gives the sex block, and the published
    Supplementary Table 2 per-sample cell counts give the unique subset partition), then
    checked against the SDRF sample list;
  * integer recoverability is measured per cell, not assumed;
  * the compartment map is an explicit dict over every author label, never a substring
    rule, and it is asserted to be total;
  * Fibro 2 (mural) and the two labels that sit inside the wrong lineage object
    (B_cell 4 = mast, B_cell 8 = pDC) are re-checked against their own marker genes.

Run (normally through atlas/01_arms_luca_lambrechts.py, steps export6149 and arm_emtab6149):
  Rscript stage1a_export6149.R <scratch>
  python  stage1a_arm_emtab6149.py <scratch>
"""
import json
import sys
from itertools import combinations

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp

import os as _os  # noqa: E402
import sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # atlas/
from common import ROOT, SCRNA, SCRATCH, DATA, RSCRIPT  # noqa: E402,F401  paths from config.json
OUT = f"{ROOT}/output/RNA/meta_dataset/stage1a"
SDRF = f"{SCRNA}/E-MTAB-6149/metadata/E-MTAB-6149.sdrf.txt"
PANEL = f"{OUT}/panel_genes.csv"
SCR = sys.argv[1] if len(sys.argv) > 1 else f"{SCRATCH}/e6149"

SOURCE_KEY = "EMTAB6149"
STUDY = "Lambrechts_Thienpont_2018"

# ======================================================================================
# Paper Supplementary Table 1, transcribed from the deposited PDF.  This is the ONLY source of histology: both SDRFs record
# Characteristics[disease] = "lung carcinoma" for every one of the 19 samples.
# ======================================================================================
SUPPL_T1 = pd.DataFrame(
    [
        # patient sex  age copd_gold tnm        stage histology_paper histology   lobe smoking
        (1, "female", 70, "GOLD 2", "pT2bN0M0", "IIA", "Squamous", "LUSC", "left upper", "active"),
        (2, "male", 86, "GOLD 2", "pT2bN0M0", "IB", "Squamous", "LUSC", "right upper", "former"),
        (3, "male", 68, "none", "pT4N2M0", "IIIB", "Adenomatous", "LUAD", "right upper", "former"),
        (4, "female", 64, "GOLD 2", "pT2aN1M0", "IIB", "Adenomatous", "LUAD", "left upper", "former"),
        (5, "male", 60, "GOLD 2", "pT1cN0M0", "IA3", "Large cell", "large_cell", "left upper", "former"),
    ],
    columns=["patient", "sex", "age", "copd_gold", "tnm", "stage", "histology_paper",
             "histology", "lobe", "smoking"],
)
# Suppl. Table 1 footnote, carried into uns verbatim in substance:
P5_NOTE = ("patient 5 is recorded as Large cell with the note 'poorly differentiated "
           "morphology, negative for TTF1 and p63', i.e. a diagnosis of exclusion")

# Paper Supplementary Table 2: cells per sample as published.  Used to derive the
# library -> sample partition; every number below was read off the deposited PDF.
PUB_CELLS = {
    (1, "core"): 424, (1, "middle"): 400, (1, "edge"): 648,
    (2, "core"): 114, (2, "middle"): 806, (2, "edge"): 1651, (2, "normal"): 271,
    (3, "core"): 2149, (3, "middle"): 3266, (3, "edge"): 3554, (3, "normal"): 4235,
    (4, "core"): 5270, (4, "middle"): 5198, (4, "edge"): 3283, (4, "normal"): 2933,
    (5, "core"): 2983, (5, "middle"): 4913, (5, "edge"): 4664, (5, "normal"): 5936,
}
PUB_CHEM = {1: "V1", 2: "V1", 3: "V2", 4: "V2", 5: "V2"}   # Suppl. Table 2, 10X VERSION

# ======================================================================================
# EXPLICIT author-label -> compartment map.  Keys are <lineage>_<dbCluster> exactly as the
# deposited objects give them, plus the two top-level labels that were never subclustered.
# NO SUBSTRING MATCHING: "B_cell_4" is mast, "B_cell_8" is pDC, "Fibro_2" is mural.  A
# substring rule on the lineage prefix would put all three in the wrong compartment.
# Columns: level1_compartment, level2_celltype, author annotation (Suppl. Table 3), low_q
# ======================================================================================
CELLTYPE_MAP = {
    "Alveolar_1":  ("epithelial_malignant", "AT2",         "cuboidal alveolar type 2 cells", False),
    "Alveolar_2":  ("epithelial_malignant", "epithelial",  "lower quality alveolar cells", True),
    "Alveolar_3":  ("epithelial_malignant", "AT1",         "flat alveolar type 1 cells", False),
    "Alveolar_4":  ("epithelial_malignant", "epithelial",  "COPD-injured alveolar cells", False),
    "Alveolar_5":  ("epithelial_malignant", "epithelial",  "respiratory epithelial cells", False),
    "Alveolar_6":  ("epithelial_malignant", "club",        "secretory club cells", False),
    "Alveolar_7":  ("epithelial_malignant", "basal",       "basal cells", False),
    "B_cell_1":    ("B_plasma",             "B",           "follicular B-cells", False),
    "B_cell_2":    ("B_plasma",             "B",           "follicular B-cells", False),
    "B_cell_3":    ("B_plasma",             "plasma",      "plasma B-cells", False),
    "B_cell_4":    ("other",                "mast",        "mast cells", False),
    "B_cell_5":    ("B_plasma",             "B",           "MALT B-cells", False),
    "B_cell_6":    ("B_plasma",             "plasma",      "plasma B-cells", False),
    "B_cell_7":    ("B_plasma",             "B",           "MALT B-cells", False),
    "B_cell_8":    ("myeloid",              "pDC",         "plasmacytoid dendritic cells", False),
    "B_cell_9":    ("other",                "erythroid",   "erythroblasts", False),
    "EC_1":        ("endothelial",          "endothelial", "normal endothelial cells", False),
    "EC_2":        ("endothelial",          "endothelial", "lower quality endothelial cell", True),
    "EC_3":        ("endothelial",          "endothelial", "tumour endothelial cells", False),
    "EC_4":        ("endothelial",          "endothelial", "tumour endothelial cells", False),
    "EC_5":        ("endothelial",          "endothelial", "normal endothelial cells", False),
    "EC_6":        ("endothelial",          "lymphatic_EC", "lymphatic endothelial cells", False),
    "Fibro_1":     ("fibroblast",           "fibroblast",  "fibroblasts", False),
    "Fibro_2":     ("mural",                "pericyte",    "fibroblasts", False),
    "Fibro_3":     ("fibroblast",           "fibroblast",  "lower quality fibroblasts", True),
    "Fibro_4":     ("fibroblast",           "fibroblast",  "fibroblasts", False),
    "Fibro_5":     ("fibroblast",           "fibroblast",  "fibroblasts", False),
    "Fibro_6":     ("fibroblast",           "fibroblast",  "normal lung fibroblasts", False),
    "Fibro_7":     ("fibroblast",           "fibroblast",  "fibroblasts", False),
    "Myeloid_1":   ("myeloid",              "macrophage",  "macrophages", False),
    "Myeloid_2":   ("myeloid",              "macrophage",  "macrophages", False),
    "Myeloid_3":   ("myeloid",              "macrophage",  "macrophages", False),
    "Myeloid_4":   ("myeloid",              "macrophage",  "macrophages", False),
    "Myeloid_5":   ("myeloid",              "DC",          "Langerhans cells", False),
    "Myeloid_6":   ("myeloid",              "macrophage",  "macrophages", False),
    "Myeloid_7":   ("myeloid",              "neutrophil",  "granulocytes", False),
    "Myeloid_8":   ("myeloid",              "macrophage",  "macrophages", False),
    "Myeloid_9":   ("myeloid",              "DC",          "monocyte-derived dendritic cells", False),
    "Myeloid_10":  ("myeloid",              "macrophage",  "macrophages", False),
    "Myeloid_11":  ("myeloid",              "macrophage",  "macrophages", False),
    "Myeloid_12":  ("myeloid",              "DC",          "cross-presenting dendritic cells", False),
    "T_cell_1":    ("T_NK",                 "CD4_T",       "CD4+ T-cells", False),
    "T_cell_2":    ("T_NK",                 "CD8_T",       "CD8+ T-cells", False),
    "T_cell_3":    ("T_NK",                 "CD4_T",       "CD4+ T-cells", False),
    "T_cell_4":    ("T_NK",                 "CD8_T",       "CD8+ T-cells", False),
    "T_cell_5":    ("T_NK",                 "CD8_T",       "CD8+ T-cells", False),
    "T_cell_6":    ("T_NK",                 "NK",          "natural killer cells", False),
    "T_cell_7":    ("T_NK",                 "Treg",        "regulatory T-cells", False),
    "T_cell_8":    ("T_NK",                 "CD8_T",       "CD8+ T-cells", False),
    "T_cell_9":    ("T_NK",                 "CD4_T",       "CD4+ T-cells", False),
    # never subclustered by the authors, so only the top-level label exists
    "Epi":         ("epithelial_malignant", "epithelial",  "epithelial cells, author subcluster not deposited", False),
    "tumor":       ("epithelial_malignant", "malignant",   "tumour cells", False),
}
# Fibro_2 carries the author label "fibroblasts" and is moved to mural on its own markers.
REASSIGNED = {"Fibro_2": "author lineage Fibro (fibroblast); reassigned to mural on "
                         "COX4I2/HIGD1B/RGS5/NOTCH3"}
LOWQ_CONF, EPI_CONF, GOOD_CONF = 0.30, 0.70, 1.00

MARKER_CHECK = {
    "Fibro_2":  ["COX4I2", "HIGD1B", "RGS5", "NOTCH3", "PDGFRB", "ACTA2", "DCN", "COL1A1"],
    "B_cell_4": ["TPSAB1", "CPA3", "KIT", "MS4A2", "MS4A1", "CD79A"],
    "B_cell_8": ["LILRA4", "IRF7", "GZMB", "IL3RA", "MS4A1", "CD79A"],
}

SITE_FROM_SDRF = {
    "tumor core": "core",
    "tumor middle (in between core and edge sample)": "middle",
    "tumor edge": "edge",
    "normal tissue adjacent to tumour": "normal",
}


def log(*a):
    """Print a message immediately (flushed)."""
    print(*a, flush=True)


# ======================================================================================
def derive_library_map(cells):
    """library index -> (patient, site), derived from the data plus Suppl. Table 2."""
    lib = cells.library_index.values
    bclen = cells.barcode.str.len().values
    probe = pd.read_csv(f"{SCR}/sexprobe.tsv", sep="\t").set_index("library_index")

    # (a) 10x chemistry from barcode length: v1 = 14 nt, v2 = 16 nt
    chem = {}
    for L in np.unique(lib):
        u = np.unique(bclen[lib == L])
        assert len(u) == 1, f"library {L} has mixed barcode lengths {u}"
        chem[L] = {14: "V1", 16: "V2"}[int(u[0])]

    counts = pd.Series(lib).value_counts().to_dict()
    sexinf = probe.inferred_sex.to_dict()

    # (b) candidate patients for each library: same chemistry AND same sex as Suppl. T1
    t1 = SUPPL_T1.set_index("patient")
    cand = {L: [p for p in t1.index
                if PUB_CHEM[p] == chem[L] and t1.loc[p, "sex"] == sexinf[L]]
            for L in counts}

    # (c) chemistry plus sex leaves P3 and P5 tied (both male, both v2).  Break that on the
    #     published per-sample cell counts, which are all distinct between the two.
    assign = {}
    for p in t1.index:
        pool = sorted(L for L in counts if p in cand[L])
        if any(len(cand[L]) > 1 for L in pool):
            sizes = {v for (pp, _), v in PUB_CELLS.items() if pp == p}
            pool = sorted(L for L in pool if counts[L] in sizes)
        assign[p] = pool
    allassigned = sorted(x for v in assign.values() for x in v)
    assert allassigned == sorted(counts), "library -> patient assignment is not a partition"

    # (d) within a patient, the unique subset partition matching the published sample sizes
    lib2sample = {}
    for p, libs in assign.items():
        want = {s: PUB_CELLS[(p, s)] for s in ("core", "middle", "edge", "normal")
                if (p, s) in PUB_CELLS}
        assert sum(want.values()) == sum(counts[L] for L in libs), \
            f"patient {p}: local cells {sum(counts[L] for L in libs)} != published {sum(want.values())}"
        sols = []
        for site, target in want.items():
            hits = [c for r in range(1, len(libs) + 1) for c in combinations(libs, r)
                    if sum(counts[L] for L in c) == target]
            sols.append((site, hits))
        # take the assignment that uses every library exactly once; assert it is unique
        good = []

        def rec(i, used, acc):
            """Depth-first search over library-to-sample assignments; collects those that use every library exactly once."""
            if i == len(sols):
                if used == set(libs):
                    good.append(dict(acc))
                return
            site, hits = sols[i]
            for h in hits:
                if not (set(h) & used):
                    acc[site] = h
                    rec(i + 1, used | set(h), acc)
                    acc.pop(site)

        rec(0, set(), {})
        assert len(good) == 1, f"patient {p}: {len(good)} library partitions match, not 1"
        for site, hs in good[0].items():
            for L in hs:
                lib2sample[L] = (p, site)
    return lib2sample, chem, counts, probe


# ======================================================================================
def main():
    """Build the E-MTAB-6149 arm from the R export: recover counts, derive the library map, harmonize obs, write the h5ad and side tables."""
    log("=" * 78)
    cells = pd.read_csv(f"{SCR}/cells.tsv", sep="\t")
    subs = pd.read_csv(f"{SCR}/subclusters.tsv", sep="\t")
    res = pd.read_csv(f"{SCR}/panel_resolution.csv")
    panel = pd.read_csv(PANEL)
    log(f"cells {len(cells)}  sub-object cells {len(subs)}")
    assert len(cells) == 52698, len(cells)
    assert cells.cell.is_unique

    # ---- counts matrix --------------------------------------------------------------
    X = sio.mmread(f"{SCR}/panel_counts.mtx").tocsr()          # genes x cells
    gid = [l.strip() for l in open(f"{SCR}/panel_genes_present.txt")]
    cid = [l.strip() for l in open(f"{SCR}/panel_cells.txt")]
    assert cid == cells.cell.tolist()
    X = X.T.tocsr()                                            # cells x genes
    assert np.allclose(X.data, np.round(X.data)), "recovered counts are not integers"
    log(f"panel matrix {X.shape}, nnz {X.nnz:,}, max {int(X.data.max())}")

    # ---- integer recoverability ------------------------------------------------------
    worst = float(cells.max_integer_deviation.max())
    n_ok = int((cells.max_integer_deviation < 1e-6).sum())
    counts_are_estimates = not (n_ok == len(cells))
    log(f"integer recovery: {n_ok}/{len(cells)} cells within 1e-6, worst deviation {worst:.3e}"
        f"  -> counts_are_estimates = {counts_are_estimates}")

    # ---- library -> patient/site -----------------------------------------------------
    lib2sample, chem, libcounts, probe = derive_library_map(cells)
    log("derived library map (independent of any hard-coded table):")

    sd = pd.read_csv(SDRF, sep="\t", dtype=str).drop_duplicates("Source Name")
    sd["patient"] = sd["Characteristics[individual]"].astype(int)
    sd["site"] = sd["Characteristics[biopsy site]"].map(SITE_FROM_SDRF)
    assert sd.site.notna().all()
    sd = sd.set_index(["patient", "site"])
    run = sd["Comment[RUN]"].to_dict()
    sdrf_sample = sd["Source Name"].to_dict()
    sdrf_sex = sd["Characteristics[sex]"].to_dict()
    sdrf_chem = sd["Comment[library construction]"].to_dict()

    rows = []
    for L in sorted(lib2sample):
        p, s = lib2sample[L]
        rows.append(dict(library_index=L, patient=p, sampling_site=s, cells=libcounts[L],
                         chem_from_barcode=chem[L], sdrf_chem=sdrf_chem[(p, s)],
                         sex_from_probe=probe.loc[L, "inferred_sex"], sdrf_sex=sdrf_sex[(p, s)],
                         sdrf_sample=sdrf_sample[(p, s)], sdrf_run=run[(p, s)]))
    libtab = pd.DataFrame(rows)
    libtab["library_id"] = [f"{r.sdrf_run}_L{r.library_index:02d}" for r in libtab.itertuples()]
    # cross-checks against the SDRF, which knows sample but not library
    assert (libtab.sex_from_probe == libtab.sdrf_sex).all(), "sex probe disagrees with SDRF"
    assert (libtab.chem_from_barcode.map({"V1": "10xV1a", "V2": "10xV2"})
            == libtab.sdrf_chem).all(), "barcode length disagrees with SDRF chemistry"
    log(libtab[["library_index", "library_id", "patient", "sampling_site", "cells",
                "chem_from_barcode", "sex_from_probe"]].to_string(index=False))

    # ---- technical replicate groups: samples that were split over two libraries -------
    per_sample = libtab.groupby(["patient", "sampling_site"], sort=False)
    trg, pairs = {}, []
    for (p, s), g in per_sample:
        if len(g) > 1:
            gid_ = f"{STUDY}_P{p}_{s}"
            for L in g.library_index:
                trg[L] = gid_
            pairs.append(dict(group=gid_, patient=p, sampling_site=s,
                              sdrf_sample=g.sdrf_sample.iloc[0], sdrf_run=g.sdrf_run.iloc[0],
                              libraries="+".join(g.library_id),
                              cells_per_library="+".join(map(str, g.cells)),
                              cells_total=int(g.cells.sum()),
                              cells_published=PUB_CELLS[(p, s)]))
    pairs = pd.DataFrame(pairs)
    log(f"\ntechnical replicate groups: {len(pairs)}")
    log(pairs.to_string(index=False))

    # barcode overlap inside a pair decides what kind of replicate this is
    bc_by_lib = {L: set(g.barcode) for L, g in cells.groupby("library_index")}
    ov = []
    for r in pairs.itertuples():
        a, b = [int(x.split("_L")[1]) for x in r.libraries.split("+")]
        i = len(bc_by_lib[a] & bc_by_lib[b])
        ov.append(dict(group=r.group, lib_a=a, lib_b=b, n_a=len(bc_by_lib[a]),
                       n_b=len(bc_by_lib[b]), shared_barcodes=i,
                       pct_of_smaller=round(100 * i / min(len(bc_by_lib[a]), len(bc_by_lib[b])), 2)))
    ov = pd.DataFrame(ov)
    log("\nbarcode overlap inside each replicate group:")
    log(ov.to_string(index=False))

    # ---- author taxonomy -------------------------------------------------------------
    subs["key"] = subs.sub_lineage + "_" + subs.sub_cluster.astype(str)
    cells = cells.merge(subs[["cell", "key", "sub_lineage", "sub_cluster"]], on="cell", how="left")
    cells["level3_celltype_native"] = cells.key.fillna(cells.lineage_native)
    unmapped = sorted(set(cells.level3_celltype_native) - set(CELLTYPE_MAP))
    assert not unmapped, f"labels with no explicit compartment: {unmapped}"
    unused = sorted(set(CELLTYPE_MAP) - set(cells.level3_celltype_native))
    log(f"\nauthor labels in object: {cells.level3_celltype_native.nunique()}; "
        f"map entries never used: {unused}")

    m = cells.level3_celltype_native.map(CELLTYPE_MAP)
    cells["level1_compartment"] = [x[0] for x in m]
    cells["level2_celltype"] = [x[1] for x in m]
    cells["level3_author_annotation"] = [x[2] for x in m]
    cells["author_low_quality"] = [x[3] for x in m]
    cells["compartment_reassigned"] = cells.level3_celltype_native.isin(REASSIGNED)
    cells["celltype_confidence"] = np.where(
        cells.author_low_quality, LOWQ_CONF,
        np.where(cells.level3_celltype_native == "Epi", EPI_CONF, GOOD_CONF))

    # ---- sample-level metadata -------------------------------------------------------
    l2p = {r.library_index: r for r in libtab.itertuples()}
    cells["patient"] = cells.library_index.map(lambda L: l2p[L].patient)
    cells["sampling_site"] = cells.library_index.map(lambda L: l2p[L].sampling_site)
    cells["library_id"] = cells.library_index.map(lambda L: l2p[L].library_id)
    cells["sample_id"] = cells.library_index.map(lambda L: l2p[L].sdrf_sample)
    cells["sdrf_run"] = cells.library_index.map(lambda L: l2p[L].sdrf_run)
    cells["technical_replicate_group"] = cells.library_index.map(lambda L: trg.get(L, ""))
    t1 = SUPPL_T1.set_index("patient")
    for c in ["sex", "age", "stage", "smoking", "histology", "histology_paper", "copd_gold", "tnm"]:
        cells[c] = cells.patient.map(t1[c])

    cells["cell_id"] = (SOURCE_KEY + "|" + cells.library_id + "|" + cells.barcode)
    assert cells.cell_id.is_unique, "cell_id is not unique"
    cells["platform"] = cells.patient.map({1: "10x 3' v1", 2: "10x 3' v1",
                                           3: "10x 3' v2", 4: "10x 3' v2", 5: "10x 3' v2"})
    cells["tissue_class"] = np.where(cells.sampling_site == "normal",
                                     "adjacent_normal", "primary_tumour")
    cells["donor_id"] = STUDY + "_P" + cells.patient.astype(str)
    cells["block_id"] = STUDY + "_" + cells.platform.str.replace("10x 3' ", "10x_3p_", regex=False)

    # ---- marker re-checks on the reassigned / cross-lineage labels --------------------
    ens2sym = dict(zip(res.loc[res.resolution != "absent", "matrix_ensembl_id"],
                       res.loc[res.resolution != "absent", "gene_symbol"]))
    sym2col = {ens2sym[e]: i for i, e in enumerate(gid)}
    det = {}
    for key, genes in MARKER_CHECK.items():
        sel = (cells.level3_celltype_native == key).values
        ref = {"Fibro_2": cells.level3_celltype_native.str.startswith("Fibro_").values & ~sel,
               "B_cell_4": cells.level3_celltype_native.isin(["B_cell_1", "B_cell_2"]).values,
               "B_cell_8": cells.level3_celltype_native.isin(["B_cell_1", "B_cell_2"]).values}[key]
        for g in genes:
            if g not in sym2col:
                continue
            col = np.asarray((X[:, sym2col[g]] > 0).todense()).ravel()
            det[(key, g)] = (round(100 * col[sel].mean(), 1), round(100 * col[ref].mean(), 1))
    log("\nmarker re-check, pct of cells detecting (label vs its lineage siblings):")
    for (k, g), (a, b) in det.items():
        log(f"   {k:9s} {g:8s} {a:6.1f}  vs {b:6.1f}")

    # ---- var: all 2,059 panel genes, absent ones as explicit structural zeros ---------
    res = res.set_index("gene_symbol").loc[panel.gene_symbol].reset_index()
    pres = (res.resolution != "absent").values
    col_of = {e: i for i, e in enumerate(gid)}
    take = np.array([col_of[e] for e in res.loc[pres, "matrix_ensembl_id"]])
    Xsub = X[:, take].astype(np.int32).tocsc()
    # widen to the full 2,059-gene panel, absent genes as explicit all-zero columns
    indptr = np.zeros(len(res) + 1, dtype=np.int64)
    indptr[1:][pres] = np.diff(Xsub.indptr)
    indptr = np.cumsum(indptr)
    Xfull = sp.csc_matrix((Xsub.data, Xsub.indices, indptr),
                          shape=(len(cells), len(res)), dtype=np.int32).tocsr()
    assert Xfull.nnz == X.nnz
    log(f"\nfull panel matrix {Xfull.shape}, nnz {Xfull.nnz:,}, "
        f"{int((~pres).sum())} genes are structural zeros")

    var = pd.DataFrame({
        "gene_symbol": res.gene_symbol.values,
        "ensembl_id_panel": res.ensembl_id.values,
        "ensembl_id_in_source": res.matrix_ensembl_id.fillna("").values,
        "gene_id_resolution": res.resolution.values,
        "present_in_source": pres,
        "source_set": res.source_set.values,
        "is_hvg": res.is_hvg.values,
    }).set_index("gene_symbol")
    var.index.name = None

    # ---- obs, in the order the schema declares --------------------------------------
    obs = pd.DataFrame(index=pd.Index(cells.cell_id.values, name=None))
    obs["source_key"] = SOURCE_KEY
    obs["study"] = STUDY
    obs["donor_id"] = cells.donor_id.values
    obs["library_id"] = cells.library_id.values
    obs["histology"] = cells.histology.values
    obs["histology_source"] = "paper_supplement"
    obs["histology_confidence"] = "high"
    obs["tissue_class"] = cells.tissue_class.values
    obs["sampling_site"] = cells.sampling_site.values
    obs["platform"] = cells.platform.values
    obs["modality"] = "single_cell"
    obs["enrichment"] = "whole_tissue"
    obs["composition_safe"] = True
    obs["treatment"] = "naive"
    obs["block_id"] = cells.block_id.values
    obs["level1_compartment"] = cells.level1_compartment.values
    obs["level2_celltype"] = cells.level2_celltype.values
    obs["level3_celltype_native"] = cells.level3_celltype_native.values
    obs["celltype_source"] = "author"
    obs["celltype_confidence"] = cells.celltype_confidence.values.astype(float)
    obs["is_malignant"] = pd.array([pd.NA] * len(cells), dtype="boolean")
    obs["malignancy_score"] = np.nan
    obs["n_counts"] = cells.n_counts.values.astype(np.int32)
    obs["n_genes"] = cells.n_genes.values.astype(np.int32)
    obs["counts_are_estimates"] = counts_are_estimates
    obs["doublet_score"] = np.nan
    obs["sex"] = cells.sex.values
    obs["age"] = cells.age.values.astype(float)
    obs["stage"] = cells.stage.values
    obs["smoking"] = cells.smoking.values
    # arm-specific extras
    obs["sample_id"] = cells.sample_id.values
    obs["sdrf_run"] = cells.sdrf_run.values
    obs["patient"] = cells.patient.values.astype(np.int16)
    obs["library_index"] = cells.library_index.values.astype(np.int16)
    obs["technical_replicate_group"] = cells.technical_replicate_group.values
    obs["barcode"] = cells.barcode.values
    obs["lineage_native"] = cells.lineage_native.values
    obs["level3_author_annotation"] = cells.level3_author_annotation.values
    obs["author_low_quality"] = cells.author_low_quality.values
    obs["compartment_reassigned"] = cells.compartment_reassigned.values
    obs["histology_paper"] = cells.histology_paper.values
    obs["copd_gold"] = cells.copd_gold.values
    obs["tnm"] = cells.tnm.values
    obs["n_counts_panel"] = cells.n_counts_panel.values.astype(np.int32)
    obs["n_genes_panel"] = cells.n_genes_panel.values.astype(np.int32)

    for c in ["source_key", "study", "histology", "histology_source", "histology_confidence",
              "tissue_class", "sampling_site", "platform", "modality", "enrichment",
              "treatment", "block_id", "level1_compartment", "level2_celltype",
              "level3_celltype_native", "celltype_source", "sex", "lineage_native",
              "level3_author_annotation", "histology_paper", "copd_gold",
              "technical_replicate_group", "donor_id", "library_id", "sample_id"]:
        obs[c] = obs[c].astype("category")

    assert obs.donor_id.astype(str).str.len().gt(0).all()
    assert (obs.composition_safe == (obs.enrichment == "whole_tissue")).all()

    A = ad.AnnData(X=Xfull, obs=obs, var=var)
    A.layers["counts"] = Xfull.copy()
    A.obsm["X_tsne_author"] = cells[["tsne1", "tsne2"]].values.astype(np.float32)
    A.uns.update({
        "source_key": SOURCE_KEY, "study": STUDY, "accession": "E-MTAB-6149",
        "stage": "1a", "gene_space": "panel_genes.csv (2,059 genes)",
        "X_is": "integer UMI counts, int32 CSR; layers['counts'] is the same matrix",
        "normalisation_note": (
            "the deposited object is named log2cpm but holds ln(1 + c*1e4/T); it was "
            "inverted per cell using s = exp(min positive value) - 1, giving integer UMI "
            f"counts for {n_ok}/{len(cells)} cells at tolerance 1e-6, worst deviation "
            f"{worst:.3e}. counts_are_estimates = {counts_are_estimates}."),
        "histology_note": ("histology is not in the SDRF or IDF, which record "
                           "Characteristics[disease] = 'lung carcinoma' for all 19 samples; "
                           "it comes from paper Supplementary Table 1. " + P5_NOTE),
        "library_map_note": ("the 24 barcode suffixes are libraries, not samples; the SDRF "
                             "lists 19 runs. The library to patient and sampling-site map "
                             "was derived from barcode length (10x chemistry), an "
                             "XIST/RPS4Y1 sex probe and the unique subset partition of "
                             "Supplementary Table 2 per-sample cell counts."),
        "genespace_note": ("26 panel genes carry a different Ensembl id in this 2018 "
                           "featuredata build than in the LuCA reference and were resolved "
                           "by unique gene symbol (TXNIP and MRC1 among them); an Ensembl-id "
                           "join alone silently drops them. 65 panel genes are absent and "
                           "are written as structural zeros; read var['present_in_source'] "
                           "before interpreting any zero."),
        "celltype_note": ("level3_celltype_native is <author lineage>_<author dbCluster>. "
                          "Fibro_2 is labelled fibroblast by the authors and is mapped to "
                          "mural here on its own markers. B_cell_4 is mast and B_cell_8 is "
                          "pDC despite living inside the B_cell lineage object. Epi cells "
                          "carry confidence 0.7 because Epithelial.Cellview.Rds was not "
                          "deposited, so the author 'lower quality epithelial cells' "
                          "subcluster (50 of 214) cannot be separated."),
        "celltype_map": pd.DataFrame(
            [(k, *v) for k, v in CELLTYPE_MAP.items()],
            columns=["level3_celltype_native", "level1_compartment", "level2_celltype",
                     "author_annotation", "author_low_quality"]),
        "reassigned_labels": dict(REASSIGNED),
        "technical_replicate_groups": pairs,
        "technical_replicate_barcode_overlap": ov,
        "library_table": libtab,
        "suppl_table_1": SUPPL_T1,
        "marker_recheck_pct_detected": pd.DataFrame(
            [(k, g, a, b) for (k, g), (a, b) in det.items()],
            columns=["level3_celltype_native", "gene", "pct_in_label", "pct_in_siblings"]),
    })

    fp = f"{OUT}/arm_emtab6149.h5ad"
    A.write_h5ad(fp, compression="gzip")
    log(f"\nwrote {fp}  {A.shape}")

    # ---- side tables and the report --------------------------------------------------
    libtab.to_csv(f"{OUT}/arm_emtab6149_libraries.csv", index=False)
    pairs.merge(ov, on="group").to_csv(f"{OUT}/arm_emtab6149_technical_replicates.csv", index=False)
    ctab = (obs.groupby(["level3_celltype_native", "level1_compartment", "level2_celltype"],
                        observed=True).size().rename("cells").reset_index())
    ctab.to_csv(f"{OUT}/arm_emtab6149_celltype_map.csv", index=False)

    log("\n== donors by histology ==")
    log(obs.groupby("histology", observed=True)
        .agg(donors=("donor_id", "nunique"), cells=("donor_id", "size")).to_string())
    log("\n== cells by compartment ==")
    log(obs.level1_compartment.value_counts().to_string())
    log("\n== cells by sampling_site ==")
    log(obs.sampling_site.value_counts().to_string())
    log("\n== compartment x histology ==")
    log(pd.crosstab(obs.level1_compartment, obs.histology).to_string())
    log("\n== fibroblast detail ==")
    f = obs[obs.lineage_native == "Fibro"]
    log(f.groupby(["level3_celltype_native", "level1_compartment"], observed=True)
        .size().rename("cells").to_string())
    log(f"Fibro lineage cells {len(f)}; kept as fibroblast "
        f"{(f.level1_compartment == 'fibroblast').sum()}; moved to mural "
        f"{(f.level1_compartment == 'mural').sum()}; author low quality "
        f"{int(f.author_low_quality.sum())}")
    log("\n== per-donor median depth ==")
    log(obs.groupby("donor_id", observed=True)
        .agg(cells=("n_counts", "size"), median_umi=("n_counts", "median"),
             median_genes=("n_genes", "median"),
             pct_fibroblast=("level1_compartment", lambda s: round(100 * (s == "fibroblast").mean(), 2)))
        .to_string())
    log("\n== fibroblast yield by sampling_site within donor ==")
    fy = pd.crosstab(obs.donor_id, obs.sampling_site,
                     values=(obs.level1_compartment == "fibroblast").astype(int),
                     aggfunc="mean") * 100
    log(fy.round(2).to_string())

    json.dump({"shape": list(A.shape),
               "counts_are_estimates": bool(counts_are_estimates),
               "integer_recovery": [n_ok, len(cells)],
               "worst_integer_deviation": worst,
               "genes_present_in_source": int(pres.sum()),
               "genes_structural_zero": int((~pres).sum()),
               "donors_by_histology": obs.groupby("histology", observed=True)
                                        .donor_id.nunique().to_dict(),
               "cells_by_compartment": obs.level1_compartment.value_counts().to_dict(),
               "cells_by_sampling_site": obs.sampling_site.value_counts().to_dict(),
               "technical_replicate_groups": pairs.to_dict("records")},
              open(f"{OUT}/arm_emtab6149_build_summary.json", "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
