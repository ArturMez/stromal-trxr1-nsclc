#!/usr/bin/env python
"""
Atlas stage 02, step build_arm: the E-MTAB-13526 arm (De Zuani et al.)
  -> output/RNA/meta_dataset/stage2/arm_emtab13526.h5ad

Paper items: part of the atlas behind Figure 2, Figure 3, Supplementary Figures 2 to 4 and Supplementary Data 1. No random numbers are drawn.
Run: python stage2_emtab13526_14_build_arm.py <scratch_dir>   (normally through atlas/02_arm_emtab13526.py)

Inputs
  * per-library panel-count CSRs written by stage2_emtab13526_13_extract.py from the
    LOCAL Matrix Market triplets (no remote payload, no cell calling);
  * the VERIFIED crosswalk from stage2_emtab13526_09_final_crosswalk.py, whose authority
    route is the (barcode, UMI-total) fingerprint against the deposited matrices;
  * the EXPLICIT 61-label compartment map from stage2_emtab13526_11_compartment_map.py.

Everything the object claims is either measured here or carried from a verified table:
  * the barcode-order alignment between the extraction and author_obs is re-derived and
    asserted row for row, not assumed;
  * platform is taken from the measured barcode whitelist, not from the SDRF, because the
    SDRF claims 10x 3' v3 for all 81 libraries while P1_B1 and P1_T1 were quantified on
    the 737,280-barcode v2 whitelist;
  * compartment mapping is a dict lookup over whole label strings, never a substring;
  * composition_safe is False for every cell: all 81 libraries are sorted fractions.

"""
import json
import os
import sys
import datetime as dt

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

import os as _os  # noqa: E402
import sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # atlas/
from common import ROOT, SCRNA, SCRATCH, DATA, RSCRIPT  # noqa: E402,F401  paths from config.json
ST2 = f"{ROOT}/output/RNA/meta_dataset/stage2"
PANEL = f"{ROOT}/output/RNA/meta_dataset/stage1a/panel_genes.csv"
RAW = f"{SCRNA}/E-MTAB-13526/raw_counts"
SCRATCH = sys.argv[1] if len(sys.argv) > 1 else f"{SCRATCH}/e13526"
EXT = f"{SCRATCH}/emtab13526_extract"

SOURCE_KEY = "EMTAB13526"
STUDY = "DeZuani_Cvejic_2024"
ACCESSION = "E-MTAB-13526"

# level2 vocabulary of the merged panel object.  The stage-2 compartment map uses the
# arm-native words on the left; the merged object uses the words on the right.
LEVEL2_TO_MERGED = {
    "B": "B_cell", "plasma": "plasma_cell", "lymphatic_EC": "endothelial_lymphatic",
    "epithelial_unspecified": "epithelial_other", "fibroblast_unspecified": "fibroblast",
    "T_unspecified": "T_other",          # new level, see uns.level2_new_levels
    "unassigned": "",
    "AT2": "AT2", "CD4_T": "CD4_T", "CD8_T": "CD8_T", "DC": "DC", "NK": "NK",
    "Treg": "Treg", "ciliated": "ciliated", "club": "club", "macrophage": "macrophage",
    "mast": "mast", "monocyte": "monocyte", "pDC": "pDC",
}
CONF = {"high": 1.0, "medium": 0.7, "low": 0.4, "drop": 0.2}

FIBROBLAST_LIMIT = (
    "This is NOT a fibroblast dataset and no LUAD versus LUSC fibroblast contrast can be "
    "run on it at any level of effort. The authors' own annotation gives 2,089 fibroblasts "
    "in the whole atlas under two labels only, 'Activated adventitial fibroblasts' and "
    "'Adventitial fibroblasts'. In primary tumour material there are 84 LUAD fibroblasts "
    "spread over 6 patients (Patient 24 = 67, Patient 22 = 8, Patient 15 = 6, and one cell "
    "each from Patients 5, 13 and 21) against 836 LUSC fibroblasts. There are no CAF "
    "subsets, no pericytes and no mural compartment in the annotation at all, and the only "
    "endothelium is lymphatic. Use this dataset for its immune contrast. If a fibroblast "
    "result is wanted, take it from another arm.")

STRUCTURAL_LIMIT = (
    "Cancer cells, fibroblasts and endothelium sit almost entirely in the CD235a-depleted "
    "libraries (fibroblast 99.5 percent, epithelial 98.6 percent, endothelial 97.7 percent "
    "of each compartment). That fraction was run on 5 LUAD and 5 LUSC patients, not 12 and "
    "8. Any epithelial or stromal contrast here is 5 a side with severe donor imbalance: "
    "89 percent of the LUAD tumour epithelium comes from Patient 24 alone. The 12 versus 8 "
    "figure applies to the IMMUNE compartment only.")

IMMUNE_CONTRAST = (
    "The reason to carry this arm: a genuine within-study, within-platform immune contrast. "
    "12 LUAD patients against 8 LUSC patients, one protocol, one CellRanger version, one "
    "reference, matched tumour and adjacent-normal material per patient.")

SORTED_LIMIT = (
    "Every one of the 81 libraries is a sorted fraction (39 CD45+, 28 CD235a-depleted, "
    "8 CD45+ plus MDSC, 6 MDSC). There is no whole-tissue library. composition_safe is "
    "False for all 896,773 cells and no cell-type proportion from this arm may be read as "
    "tissue composition, not even within one fraction.")


def log(*a):
    """Print a message immediately (flushed)."""
    print(*a, flush=True)


def h5_safe(df):
    """h5py cannot write a mixed object column, so NaN in a text column becomes ''."""
    out = df.copy().reset_index(drop=True)
    for c in out.columns:
        if out[c].dtype == object or pd.api.types.is_string_dtype(out[c]):
            out[c] = np.asarray(out[c].fillna("").astype(str), dtype=object)
    return out


def main():
    """Assemble the arm from the extracts, crosswalk and compartment map, write the h5ad and print the report."""
    panel = pd.read_csv(PANEL)
    cw = pd.read_csv(f"{ST2}/crosswalk.csv")
    cmap = pd.read_csv(f"{ST2}/compartment_map_emtab13526.csv")
    res = pd.read_csv(f"{ST2}/arm_emtab13526_gene_resolution.csv")
    xlog = pd.read_csv(f"{ST2}/arm_emtab13526_extract_log.csv")
    idx = pd.read_parquet(f"{EXT}/_index.parquet")

    # ---- author obs, re-derived in exactly the extraction's row order ----------------
    ob = pd.read_parquet(f"{ST2}/author_obs.parquet")
    lib_of = dict(zip(cw.patient_sample, cw.library_id))
    ob["library_id"] = ob.patient_sample.map(lib_of)
    assert ob.library_id.notna().all()
    ob["barcode"] = ob.author_barcode.str.split("-").str[:2].str.join("-")
    ob = ob.iloc[np.argsort(ob.library_id.values, kind="stable")].reset_index(drop=True)
    assert (ob.library_id.values == idx.library_id.values).all(), "library order drifted"
    assert (ob.barcode.values == idx.barcode.values).all(), "barcode order drifted"
    assert (ob.author_barcode.values == idx.author_barcode.values).all()
    n = len(ob)
    log(f"author cells {n:,} across {ob.library_id.nunique()} libraries")

    # ---- counts ----------------------------------------------------------------------
    libs = sorted(ob.library_id.unique())
    blocks, ncnt, ngen = [], [], []
    for lib in libs:
        z = np.load(f"{EXT}/{lib}.npz")
        X = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"]))
        assert X.shape == (int((ob.library_id == lib).sum()), len(panel)), f"{lib} shape"
        blocks.append(X)
        ncnt.append(z["n_counts"])
        ngen.append(z["n_genes"])
    X = sp.vstack(blocks, format="csr").astype(np.int32)
    del blocks
    n_counts = np.concatenate(ncnt)
    n_genes = np.concatenate(ngen)
    assert X.shape[0] == n
    assert X.nnz == int(xlog.panel_entries.sum())
    assert (X.data > 0).all(), "a stored panel entry is not a positive integer count"
    log(f"panel matrix {X.shape}  nnz {X.nnz:,}  panel UMI {int(X.data.sum()):,}  "
        f"full-transcriptome UMI {int(n_counts.sum()):,}")

    n_counts_panel = np.asarray(X.sum(axis=1)).ravel().astype(np.int64)
    n_genes_panel = np.diff(X.indptr).astype(np.int32)

    # sanity: the author's own total_counts is computed over their 25,718-gene subset, so
    # it must not exceed the full 33,538-feature total measured here.
    ratio = ob.total_counts.values / n_counts
    log(f"author total_counts / measured full total: min {ratio.min():.3f} "
        f"median {np.median(ratio):.3f} max {ratio.max():.3f}; "
        f"cells above 1.0: {int((ratio > 1.0).sum())}")
    assert ratio.max() <= 1.0 + 1e-9, "author UMI total exceeds the raw matrix total"

    # ---- per-library metadata ---------------------------------------------------------
    cwx = cw.set_index("library_id")
    # platform: measured whitelist, not the SDRF claim
    wl = xlog.set_index("library")["whitelist"].to_dict()
    WL2PLAT = {6794880: "10x 3' v3", 737280: "10x 3' v2"}
    plat = {L: WL2PLAT[wl[L]] for L in libs}
    disputed_plat = sorted(L for L in libs if plat[L] != cwx.loc[L, "library_construction_sdrf"])
    log(f"platform from barcode whitelist; disagrees with the SDRF for: {disputed_plat}")

    ind = cwx["individual"].to_dict()
    def donor_short(s):
        """Short donor label (P<n> for patients, D<n> for donors) from an author name such as "Patient 4"."""
        kind, num = s.split()
        return ("P" if kind == "Patient" else "D") + num
    donor = {L: f"{STUDY}_{donor_short(ind[L])}" for L in libs}

    # technical replicates: same individual, same tissue class, same sorted fraction
    key = cwx.apply(lambda r: f"{donor_short(r.individual)}_{r.tissue_class}_{r.enrichment}",
                    axis=1)
    grp_size = key.value_counts()
    trg = {L: (f"{STUDY}_{key[L]}" if grp_size[key[L]] > 1 else "") for L in libs}
    treps = (pd.DataFrame({"library_id": libs, "group": [trg[L] for L in libs]})
             .query("group != ''").groupby("group").library_id
             .apply(lambda s: "+".join(sorted(s))).reset_index())
    log(f"technical replicate groups (same donor, site and fraction): {len(treps)}")

    L = ob.library_id.values
    def perlib(d, dtype=object):
        """Expand a per-library dict d to one value per cell, in the cell order of L."""
        return np.array([d[x] for x in L], dtype=dtype)

    # ---- cell types: EXPLICIT whole-string map, asserted total ------------------------
    ct = ob["Cell types"].astype(str).values
    m1 = dict(zip(cmap.author_cell_type, cmap.level1_compartment))
    m2 = dict(zip(cmap.author_cell_type, cmap.level2_celltype))
    mc = dict(zip(cmap.author_cell_type, cmap.confidence))
    missing = sorted(set(ct) - set(m1))
    assert not missing, f"author labels with no explicit compartment: {missing}"
    unused = sorted(set(m1) - set(ct))
    log(f"author labels {len(set(ct))}, all mapped explicitly; unused map rows: {unused}")

    lvl1 = np.array([m1[c] for c in ct], dtype=object)
    lvl2_native = np.array([m2[c] for c in ct], dtype=object)
    lvl2 = np.array([LEVEL2_TO_MERGED[x] for x in lvl2_native], dtype=object)
    confw = np.array([mc[c] for c in ct], dtype=object)
    conf = np.array([CONF[c] for c in confw], dtype=float)
    author_remove = np.array(["(to remove)" in c for c in ct])
    assert set(ct[author_remove]) == {c for c in set(ct) if "(to remove)" in c}
    # macrophage_strict: unambiguous mononuclear phagocyte labels only.  The three
    # low-confidence myeloid labels (STAB+ anti-inflammatory cells, CAMLs, Collagen+
    # macrophages) and the author's own 'to remove' clusters are excluded.
    mstrict = np.isin(lvl2, ["macrophage", "monocyte"]) & (conf >= 0.7) & ~author_remove

    # ---- obs ---------------------------------------------------------------------------
    cell_id = np.char.add(np.char.add(np.char.add(
        np.full(n, SOURCE_KEY + "|", dtype=object), L.astype(object)),
        np.full(n, "|", dtype=object)), ob.barcode.values.astype(object))
    cell_id = pd.Index(cell_id.astype(str))
    assert cell_id.is_unique, "cell_id is not unique"

    hist = perlib(cwx["histology"].to_dict())
    tclass = perlib(cwx["tissue_class"].to_dict())
    enrich = perlib(cwx["enrichment"].to_dict())
    tdisp = perlib(cwx["tissue_class_disputed"].astype(bool).to_dict()).astype(bool)
    edisp = perlib(cwx["enrichment_disputed"].astype(bool).to_dict()).astype(bool)
    indv = perlib(ind)
    # Patient 4: SDRF says lung squamous cell carcinoma, the author h5ad labels all of its
    # tumor cells 'Squamous dysplasia', and the tumor/adjacent-normal assignment of its
    # six libraries is contradictory between the two sources and was not resolvable.
    hconf = np.where(indv == "Patient 4", "disputed", "high")

    grading = cwx["grading_sdrf"].fillna("").replace("not applicable", "").to_dict()
    smoke = ob["smoking history"].astype(str).replace({"NA": "unknown"}).values

    obs = pd.DataFrame(index=cell_id)
    obs["source_key"] = SOURCE_KEY
    obs["study"] = STUDY
    obs["donor_id"] = perlib(donor)
    obs["library_id"] = L
    obs["histology"] = hist
    obs["histology_source"] = "SDRF"
    obs["histology_confidence"] = hconf
    obs["tissue_class"] = tclass
    obs["sampling_site"] = "NA"
    obs["platform"] = perlib(plat)
    obs["modality"] = "single_cell"
    obs["enrichment"] = enrich
    obs["composition_safe"] = False
    obs["treatment"] = "unknown"
    obs["block_id"] = np.char.add(STUDY + "_", perlib(plat).astype(str))
    obs["level1_compartment"] = lvl1
    obs["level2_celltype"] = lvl2
    obs["level3_celltype_native"] = ct
    obs["celltype_source"] = "author"
    obs["celltype_confidence"] = conf
    obs["is_malignant"] = pd.array([pd.NA] * n, dtype="boolean")
    obs["malignancy_score"] = np.nan
    obs["n_counts"] = n_counts.astype(np.int64)
    obs["n_genes"] = n_genes.astype(np.int32)
    obs["counts_are_estimates"] = False
    obs["doublet_score"] = ob.doublet_scores.values.astype(np.float32)
    obs["sex"] = perlib(cwx["sex_sdrf"].to_dict())
    obs["age"] = perlib(cwx["age_sdrf"].to_dict()).astype(float)
    obs["stage"] = ""
    obs["smoking"] = smoke
    # arm-specific
    obs["donor_id_native"] = indv
    obs["patient_id"] = perlib(donor)
    obs["tnm"] = perlib(grading)
    obs["technical_replicate_group"] = perlib(trg)
    obs["macrophage_strict"] = mstrict
    obs["barcode"] = np.asarray(ob.barcode.astype(str), dtype=object)
    obs["author_object"] = ob.author_object.values
    obs["author_barcode"] = np.asarray(ob.author_barcode.astype(str), dtype=object)
    obs["author_sample"] = ob["sample"].astype(str).values
    obs["author_sorting"] = ob.sorting.astype(str).values
    obs["author_environment"] = ob.environment.astype(str).values
    obs["sangerID"] = ob.sangerID.astype(str).values
    obs["level2_celltype_map_native"] = lvl2_native
    obs["celltype_confidence_word"] = confw
    obs["author_flagged_remove"] = author_remove
    obs["histology_disputed"] = indv == "Patient 4"
    obs["tissue_class_disputed"] = tdisp
    obs["enrichment_disputed"] = edisp
    obs["platform_disputed"] = perlib({x: x in set(disputed_plat) for x in libs}).astype(bool)
    obs["n_counts_panel"] = n_counts_panel
    obs["n_genes_panel"] = n_genes_panel
    obs["author_total_counts"] = ob.total_counts.values.astype(np.int64)
    obs["author_n_genes"] = ob.n_genes.values.astype(np.int32)
    obs["pct_counts_mt"] = ob.pct_counts_mt.values.astype(np.float32)

    for c in ["source_key", "study", "donor_id", "library_id", "histology",
              "histology_source", "histology_confidence", "tissue_class", "sampling_site",
              "platform", "modality", "enrichment", "treatment", "block_id",
              "level1_compartment", "level2_celltype", "level3_celltype_native",
              "celltype_source", "sex", "stage", "smoking", "donor_id_native",
              "patient_id", "tnm", "technical_replicate_group", "author_object",
              "author_sample", "author_sorting", "author_environment",
              "level2_celltype_map_native", "celltype_confidence_word"]:
        obs[c] = obs[c].astype("category")

    # parquet hands back pyarrow-backed strings, which anndata < 0.11 cannot write
    for c in obs.columns:
        if isinstance(obs[c].dtype, pd.CategoricalDtype):
            continue
        if pd.api.types.is_string_dtype(obs[c]) and obs[c].dtype != object:
            obs[c] = np.asarray(obs[c].astype(str), dtype=object)

    assert not obs.composition_safe.any(), "composition_safe must be False everywhere"
    assert (obs.enrichment != "whole_tissue").all()
    assert obs.level1_compartment.astype(str).ne("").all()

    # ---- var -----------------------------------------------------------------------------
    assert (res.gene_symbol.values == panel.gene_symbol.values).all()
    var = pd.DataFrame({
        "gene_symbol": panel.gene_symbol.values,
        "ensembl_id_panel": panel.ensembl_id.values,
        "ensembl_id_in_source": np.where(res.present_in_source, panel.ensembl_id.values, ""),
        "gene_id_resolution": np.where(res.present_in_source, "ensembl_id", "absent"),
        "present_in_source": res.present_in_source.values,
        "source_symbol_in_reference": res.source_symbol.fillna("").values,
        "symbol_agrees_with_reference": res.symbol_agrees.values,
        "source_set": panel.source_set.values,
        "is_hvg": panel.is_hvg.values,
    }).set_index("gene_symbol")
    var.index.name = None
    absent = var.index[~var.present_in_source].tolist()
    zero_cols = np.where(np.diff(X.tocsc().indptr) == 0)[0]
    log(f"var {var.shape}; absent from the deposit and written as structural zeros: {absent}")
    assert set(absent) <= set(var.index[zero_cols]), "an absent gene carries counts"

    # ---- library table -------------------------------------------------------------------
    libtab = cwx.loc[libs, ["individual", "patient_sample", "sangerID", "histology",
                            "tissue_class", "enrichment", "facs_sdrf", "sampling_site_sdrf",
                            "library_construction_sdrf", "author_n_cells", "ena_sample"]].copy()
    libtab["platform_measured"] = [plat[x] for x in libs]
    libtab["barcode_whitelist_size"] = [wl[x] for x in libs]
    libtab["donor_id"] = [donor[x] for x in libs]
    libtab["technical_replicate_group"] = [trg[x] for x in libs]
    libtab["cells_in_arm"] = obs.library_id.value_counts().reindex(libs).values
    libtab = libtab.reset_index()
    assert (libtab.cells_in_arm.values == libtab.author_n_cells.values).all()
    libtab.to_csv(f"{ST2}/arm_emtab13526_libraries.csv", index=False)

    # ---- assemble --------------------------------------------------------------------------
    immc = obs[obs.level1_compartment.astype(str).isin(["T_NK", "myeloid", "B_plasma"])]
    it = pd.crosstab(immc.histology.astype(str), immc.tissue_class.astype(str))
    imm_txt = (" Measured in this object, T_NK plus myeloid plus B_plasma cells: "
               + "; ".join(f"{h} {it.loc[h, 'primary_tumour']:,} tumour and "
                           f"{it.loc[h, 'adjacent_normal']:,} adjacent normal"
                           for h in ["LUAD", "LUSC"]) + ".")

    A = ad.AnnData(X=X, obs=obs, var=var)
    A.uns.update({
        "source_key": SOURCE_KEY, "study": STUDY, "accession": ACCESSION, "stage": "2",
        "study_full": ("De Zuani et al., Cvejic laboratory, Wellcome Sanger Institute, 2024. "
                       "Single cell RNA-seq atlas of human NSCLC lesions and non-involved "
                       "tissue. ArrayExpress E-MTAB-13526, secondary accession ERP153272."),
        "built_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "script": "atlas/steps/stage2_emtab13526_13_extract.py + _14_build_arm.py",
        "X_is": ("integer UMI counts, int32 CSR, panel genes only. Extracted from the "
                 "deposited raw Matrix Market triplets in original_data/external_scRNA/"
                 "E-MTAB-13526/raw_counts, subset to the barcodes the authors called as "
                 "cells. No cell calling was run and no remote payload was downloaded."),
        "gene_space": (f"panel_genes.csv, {len(panel)} genes. {int(res.present_in_source.sum())} "
                       f"resolve against the deposit by Ensembl id; {len(absent)} are absent "
                       f"and are written as explicit all-zero columns: {absent}. Read "
                       "var.present_in_source before interpreting any zero."),
        "gene_id_convention": (
            "The deposit ships one features.tsv.gz per library, all 81 byte identical, "
            "CellRanger GRCh38-3.0.0, 33,538 rows of unversioned Ensembl id, HGNC symbol, "
            "'Gene Expression'. Ensembl id is the join key. Symbols alone would be wrong: "
            f"{int((res.present_in_source & ~res.symbol_agrees).sum())} panel genes carry a "
            "different symbol in that 2018-era reference than in the panel (CCN1 is CYR61, "
            "CCN2 is CTGF, H4C3 is HIST1H4C, SEPTIN8 is SEPT8, VARS1 is VARS, and so on), "
            "and 24 symbols are duplicated in the reference."),
        "crosswalk": (
            "Author cells were assigned to the 81 SDRF libraries by the fingerprint route of "
            "stage2_emtab13526_09_final_crosswalk.py: the author's total_counts is an exact "
            "integer UMI total, so (barcode, UMI total) identifies the source matrix. All 81 "
            "author groups resolve bijectively, minimum score 0.9993 against a maximum "
            "runner-up of 0.1120. The ordinal-name route is WRONG for all six Patient 4 "
            "libraries, where the author's Tumour and Background labels are swapped against "
            "the SDRF."),
        "fibroblast_limitation": FIBROBLAST_LIMIT,
        "structural_compartment_limitation": STRUCTURAL_LIMIT,
        "immune_contrast": IMMUNE_CONTRAST + imm_txt,
        "sorted_fraction_limitation": SORTED_LIMIT,
        "platform_note": (
            "The SDRF records 10x 3' v3, cell barcode 16 nt, UMI 12 nt for all 81 libraries. "
            "The deposited matrices contradict that for Patient 1: P1_B1 and P1_T1 have "
            "exactly 737,280 columns, the 10x v2 737K-august-2016 whitelist, while the other "
            "79 have 6,794,880, the v3 3M-february-2018 whitelist. The two whitelists were "
            "read from the files here: each group is byte identical within itself and the two "
            "share only 77,142 barcodes, 10.46 percent of the v2 list, so this is a different "
            "barcode space and not a truncation. platform is set from the measured whitelist, "
            "P1 gets platform '10x 3' v2' and its own block_id, and obs.platform_disputed "
            f"marks its {int(obs.platform_disputed.sum()):,} cells. The underlying reads were "
            "not re-examined."),
        "histology_note": (
            "histology comes from SDRF Characteristics[disease], never from the author "
            "obs['tumour type'], which has 14 dirty levels including trailing-space duplicates "
            "and an embedded newline. Patient 4 is 'lung squamous cell carcinoma' in the SDRF "
            "and 'Squamous dysplasia' in the author object, so all of its cells carry "
            "histology_confidence 'disputed' and histology_disputed True."),
        "tissue_class_note": (
            "tissue_class comes from SDRF Characteristics[sampling site]. For Patient 4 the "
            "author object reverses tumour and adjacent normal on all six libraries and the "
            "conflict could not be resolved from either source, so those cells carry "
            "tissue_class_disputed True. Author labels are kept verbatim in "
            "obs.author_environment so the alternative reading stays available."),
        "enrichment_note": (
            "enrichment comes from SDRF Characteristics[FACS]. The author 'sorting' field "
            "disagrees for Donor 2 libraries D2_3 and D2_4; the SDRF is followed because the "
            "cell composition supports it (D2_3, SDRF CD45+, is 33.7 percent lymphocyte "
            "against 0.7 percent in D2_4, SDRF MDSC). Those cells carry enrichment_disputed."),
        "celltype_note": (
            "level3_celltype_native is the author 'Cell types' column verbatim, 61 labels, "
            "including 'Cycling anti-inflammatory alveolar macrophages ' with its trailing "
            "space, which is the only form that occurs and must not be stripped. The map to "
            "level1_compartment is an explicit dict over whole label strings and is asserted "
            "total over all 896,773 cells; no substring rule is used anywhere, because "
            "'malignant cell' and 'mast cell' both contain 't cell'. Malignant status is not "
            "annotated by the authors, so is_malignant is left empty for the whole arm."),
        "celltype_low_confidence": (
            "Three myeloid labels are mapped on judgement and carry celltype_confidence 0.4: "
            "'STAB+ anti-inflammatory cells' (14,059 cells; STAB1 is a scavenger-receptor "
            "macrophage marker but the label omits the lineage word), 'CAMLs' (2,520; myeloid "
            "in the literature but a fused giant-cell entity, 95.9 percent in CD235a-), and "
            "'Collagen+ macrophages' (1,131; collagen in a macrophage is a classic "
            "fibroblast-doublet or ambient signature, 93.8 percent in CD235a-). '8,1' (128 "
            "cells) is a bare leiden id with no name and gets level2 empty. Four labels the "
            "authors marked '(to remove)' (839 cells) are mapped but flagged in "
            "obs.author_flagged_remove."),
        "level2_new_levels": (
            "level2_celltype uses the merged panel-object vocabulary. One level is new: "
            "'T_other', for the author's 'Naive T cells' where CD4 and CD8 are not resolved "
            "and for 'Gamma delta T cells' and 'Downregulated T cells', which have no slot in "
            "the existing vocabulary. stage1a_merge_panel_object.LEVEL2_MAP needs the entry "
            "'T_other': 'T_other' added before this arm is merged. The arm-native word is kept "
            "in obs.level2_celltype_map_native."),
        "macrophage_strict_note": (
            "macrophage_strict is True for author labels mapped to macrophage or monocyte at "
            "confidence medium or better, excluding the three low-confidence labels and the "
            "author's 'to remove' clusters."),
        "depth_note": (
            "n_counts and n_genes are full-transcriptome totals over all 33,538 deposited "
            "features, accumulated during extraction before any gene subsetting. "
            "author_total_counts is the authors' own total over their 25,718-gene object and "
            "is always the smaller number; n_counts_panel and n_genes_panel are the totals "
            "over the panel genes in X."),
        "compartment_map": h5_safe(cmap),
        "library_table": h5_safe(libtab),
        "technical_replicate_groups": h5_safe(treps),
        "extract_log": h5_safe(xlog),
        "gene_resolution": h5_safe(res),
    })

    fp = f"{ST2}/arm_emtab13526.h5ad"
    A.write_h5ad(fp, compression="gzip")
    log(f"\nwrote {fp}  {A.shape}  {os.path.getsize(fp)/1e9:.2f} GB")

    # ---- report ---------------------------------------------------------------------------
    o = A.obs
    don = o.groupby("histology", observed=True).donor_id.nunique()
    cells = o.histology.value_counts()
    rep = pd.DataFrame({"donors": don, "cells": cells}).fillna(0).astype(int)
    log("\ncells and donors by histology\n" + rep.to_string())

    ct1 = pd.crosstab(o.level1_compartment, o.histology)
    log("\ncells by compartment by histology\n" + ct1.to_string())

    log("\ncells per sorted fraction\n" +
        pd.DataFrame({"cells": o.enrichment.value_counts(),
                      "libraries": libtab.groupby("enrichment").size()}).to_string())

    fib = o[o.level1_compartment == "fibroblast"]
    log("\nfibroblasts by histology and tissue class\n" +
        pd.crosstab(fib.histology, fib.tissue_class).to_string())
    ftum = fib[fib.tissue_class == "primary_tumour"]
    log("\ntumour fibroblasts per donor\n" +
        ftum.groupby(["histology", "donor_id"], observed=True).size()
        .rename("cells").reset_index().query("cells > 0").to_string(index=False))

    summ = dict(
        shape=list(A.shape), nnz=int(A.X.nnz), libraries=len(libs),
        donors=int(o.donor_id.nunique()),
        donors_by_histology=don.to_dict(), cells_by_histology=cells.to_dict(),
        cells_by_compartment=o.level1_compartment.value_counts().to_dict(),
        cells_by_enrichment=o.enrichment.value_counts().to_dict(),
        fibroblasts_total=int(len(fib)),
        fibroblasts_tumour_by_histology=ftum.histology.value_counts().to_dict(),
        genes_resolved=int(res.present_in_source.sum()), genes_absent=absent,
    )
    with open(f"{ST2}/arm_emtab13526_build_summary.json", "w") as fh:
        json.dump(summ, fh, indent=2, default=str)
    ct1.to_csv(f"{ST2}/arm_emtab13526_compartment_by_histology.csv")
    log("\nOK")


if __name__ == "__main__":
    main()
