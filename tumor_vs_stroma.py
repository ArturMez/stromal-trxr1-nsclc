#!/usr/bin/env python3
"""
TrxR1 in cancer cells and in stromal cells, TrxR1 by stromal cell type, and the genotype and smoking tests.

WHAT: the analyses behind the Results section "TrxR1 in cancer cells and in stromal cells", Figure 1d,
      Supplementary Table 1a and Supplementary Figure 1: the TrxR1-positive fraction of cancer cells
      and of stromal cells per patient, where the positivity cut sits, how the cancer-cell fraction
      relates to survival, to the KEAP1, NFE2L2 and CUL3 genotype and to bulk tumor TXNRD1 RNA, and
      how closely the two fractions follow each other once TMA slide and background staining are
      removed; the TrxR1-positive fraction of each stromal cell type (Figure 1d) and the survival
      block by cell type of Supplementary Table 1a. Also the comparison of current against former
      smokers behind the smoking sentence of the thioredoxin section.
FROM: work_root/intermediate/<cohort>/{cells.tsv, clinical.csv, clinical_full.csv, thresholds.csv};
      the TrxR1 channel of the raw cell export, through redox.staining; under raw_root/Bomi2 the
      mutation summary (genotype of Uppsala II), the Swedish clinical database (smoking on three
      levels, Uppsala II) and the 1995 to 2005 master sheet (its code book, and a check of the
      Uppsala I smoking codes); the bulk RNA table of the Uppsala II tumors (GSE81089), named by the
      optional "gse81089" entry of config.json.
TO:   work_root/output/tumor_vs_stroma/, each CSV with a sibling .md, one line per column:
        percase_BOMI2.csv, percase_BOMI1.csv   s1  one row per patient: the TrxR1-positive
                                                   fractions and cell counts, background indices,
                                                   slide, genotype, smoking, clinical variables
                                                   (patient-level, not shared)
        cut_position.csv                        s1  where the cut sits among cancer and stromal cells
        density.csv                             s1  the intensity densities of Supplementary Figure 1a
        imaging_driver_tests.csv                s2  the genotype tests (KEAP1, NFE2L2 or CUL3 altered)
                                                   of Supplementary Figure 1d, Mann-Whitney
        smoking_detail.csv                      s2  current against former smokers, stromal fraction
                                                   (Mann-Whitney)
        distribution.csv                        s3  per cohort and histology, the spread of each fraction
        cancer_survival.csv                     s3  the cancer-cell fraction against overall survival
        joint_models.csv                        s3  stromal and cancer-cell fractions in one model
        bulk_rna.csv, bulk_per_patient.csv      s3  bulk TXNRD1 against the two fractions, genotype
                                                   and survival (skipped without the gse81089 entry)
        compartment_correlation.csv             s4  cancer-cell against stromal fraction, and T cells
                                                   against fibroblasts, raw, within slide, and within
                                                   slide with background staining removed
        compartment_points.csv                  s4  the per-patient points of Supplementary Figure 1e
        supptable1a.csv, suppfig1_numbers.csv   s5  Supplementary Table 1a, one row per printed cell of
                                                   its blocks other than the survival block by cell type
                                                   and the covariate block; the counts of the
                                                   Supplementary Figure 1 legend and the bulk RNA count
                                                   of Figure 1a
        suppfig1a_density ... suppfig1f_bulk    s6  the six panels of Supplementary Figure 1, PNG and PDF
        lineage_profile_patients.csv            s7  per patient, the TrxR1-positive fraction of each
                                                   cell type (patient-level, not shared)
        lineage_profile_summary.csv             s7  its median and quartiles over both histologies,
                                                   Uppsala II: Figure 1d
        fig1d_lineage_profile.png/pdf           s7  Figure 1d, drawn from lineage_profile_summary.csv
        lineage_survival.csv                    s8  overall survival against the TrxR1-positive fraction
                                                   of all cells and of each stromal cell type
        supptable1a_celltypes.csv               s8  the printed cells of that survival block of
                                                   Supplementary Table 1a, with its patients (deaths) row
        covariate_survival.csv                  s9  the Figure 1b model adjusted for stage and for
                                                   stage, age, sex, WHO performance status and TMA
                                                   slide together (sensitivity analysis); the same
                                                   models of niche connectivity (Figure 5g) are in
                                                   chapter4.py s6 (Supplementary Table 3d)
        supptable1a_covariates.csv              s9  the printed cells of the covariate block of
                                                   Supplementary Table 1a
      Supplementary Table 1a is therefore spread over three tables: supptable1a.csv,
      supptable1a_celltypes.csv (the survival block by cell type) and supptable1a_covariates.csv (the
      covariate block). The table holds only the rows the text or a legend reads and no figure
      shows: the genotype block is Supplementary Figure 1d, the bulk RNA correlations
      Supplementary Figure 1f and the cancer-with-stromal correlations Supplementary Figure 1e. The
      bulk RNA rows for the genotype difference and survival are in the table, because the text
      states them ("four to five times higher", "not associated with survival") and panel f does
      not draw them.
WHY:  one script computes every quantity of the cancer-cell against stromal-cell comparison from
      the package inputs, writes each table before it draws from it, and asserts the checks on the
      drawn values (the density crossing lies at the cut; the adjusted correlations drawn are the
      ones in the table).

Definitions:
  TrxR1-positive   the call of cells.tsv (cytoplasmic intensity above the cohort cut)
  cancer cells     cytokeratin-lineage cells (CK, not CD3, not CD68) in the tumor compartment
  fractions        the mean call over a patient's cells that carry one; a patient enters with at least
                   20 cells of the group (cancer cells, stromal cells, or one stromal cell type in s7
                   and s8), counting every cell of the group
  background index the mean intensity of a group's TrxR1-negative cells, at least 10 of them; for the
                   slide adjustment a missing call counts as negative, in the fraction as well
  slide            the TMA slide that holds most of a patient's cells (redox.staining.patient_slide)
  survival         Cox models with Efron ties, hazard ratio per 1 SD over the analyzed patients,
                   adjusted for stage I to IV as one linear term (every survival model here)
  genotype         altered: "Mutated", "Splicing" or "Splicing:Mutated" in the mutation summary;
                   the pathway column is 1 when KEAP1, NFE2L2 or CUL3 is altered

Statistics: the genotype tests (Supplementary Figure 1d) and the smoking comparison are two-sided
Mann-Whitney tests (redox.ranktests.mw_test, scipy, method "auto"); they use no random numbers. The
survival models are Cox models with Wald p values. The only random numbers are the horizontal
jitter of the dots of panels b and d (drawing only), from numpy default_rng with
  SEED_JITTER = 3
In the tables, NRF2path is the genomic flag of the genotype: 1 when KEAP1, NFE2L2 or CUL3 is altered
in the mutation summary.

RUN: python3 tumor_vs_stroma.py [--only s1 s2 s3 s4 s5 s6 s7 s8 s9]
     s2 to s5 read the tables of the sections before them, s6 draws from all of them, s7 reads
     cells.tsv and percase_<cohort>.csv, s8 the tables of s1 and s7, s9 those of s1 and s7 plus the
     age, sex and WHO performance status columns of clinical_full.csv (redox.clinical.load_covariates).
     s3 checks its stromal models against chapter1/fig1b_models.csv when chapter1.py has written it.
"""
import argparse
import math

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, norm, rankdata, spearmanr

from redox import cells as cellsmod
from redox import clinical as clinicalmod
from redox import paths, plotting, ranktests, staining, survival

FOLDER = "tumor_vs_stroma"

SEED_JITTER = 3

# Least cells for a patient to carry a fraction: stromal cells, cancer cells.
MIN_STROMA, MIN_CANCER = cellsmod.MIN_STROMAL_CELLS, 20

HISTS = ["LUAD", "LUSC"]

# Raw files read here, under raw_root/Bomi2. The clinical database is the one mif/01_clinical.py
# merges; this script reads a different sheet of it, the cleared one, for the three smoking levels.
MUTATION_FILE = "Mutation_summary_BOMI-2_cleaned_original.xlsx"
SWEDISH_FILE = paths.input_file("bomi2_clinical_database")
SWEDISH_SHEET = "Kliniska data cleared"
SWEDISH_ID = "Projektnr:"
SWEDISH_SMOKING = "Smoking: EB, HK-om data ej fanns -->RCC "
MASTER_FILE = "NSCLC 95-05_ 170215.xlsx"

# Genes read from the mutation summary, and the entries that count as altered. Any other entry
# stops the run, so that a new spelling in the file cannot pass as wild type.
GENES = ["KEAP1", "NFE2L2", "CUL3"]
ALTERED = {"Mutated", "Splicing", "Splicing:Mutated"}

# Bulk RNA: the Ensembl identifier of TXNRD1 in the GSE81089 table.
TXNRD1_ID = "ENSG00000198431"

# Colors of Supplementary Figure 1 that belong to this figure alone: a warm neutral for cancer cells
# and the tumor compartment, a blue for stromal cells, and the axis gray.
TUM = "#b07a3b"
STROMA = "#2a78d6"
AXIS = "#c3c2b7"
FS = 5.8                 # the smallest font size of the panels


# ---------------------------------------------------------------- shared helpers

def out(name=None):
    """work_root/output/tumor_vs_stroma/<name>, creating the folder."""
    folder = paths.WORK_ROOT / "output" / FOLDER
    folder.mkdir(parents=True, exist_ok=True)
    return folder if name is None else folder / name


def save(fig, name):
    """A panel as PNG at 300 dpi and as PDF, as redox.plotting.save writes it.

    redox.paths.output names every output folder chapter<N>; this script writes into its own folder,
    so the saver is repeated here with the same crop rule (axis labels passed as extra artists).
    """
    import matplotlib.pyplot as plt
    labels = [lab for ax in fig.axes for lab in (ax.xaxis.label, ax.yaxis.label)
              if lab.get_visible() and lab.get_text()]
    extra = fig.get_default_bbox_extra_artists() + labels
    written = []
    for ext in ("png", "pdf"):
        path = out(f"{name}.{ext}")
        fig.savefig(path, dpi=300, facecolor="white", bbox_inches="tight", bbox_extra_artists=extra)
        written.append(path)
    plt.close(fig)
    for path in written:
        print(f"    -> {path}")


def _sidecar(name, what, columns):
    """Write the sibling .md that describes an output table: a sentence and one line per column."""
    lines = [f"# {name}.csv", "", what, "", "| column | meaning |", "|---|---|"]
    lines += [f"| {c} | {m} |" for c, m in columns] + [""]
    out(f"{name}.md").write_text("\n".join(lines))


def _write(frame, name, what, columns):
    """Write a table and its sidecar, checking that the sidecar names every column."""
    missing = [c for c in frame.columns if c not in dict(columns)]
    if missing:
        raise ValueError(f"{name}: the sidecar does not describe {missing}")
    path = out(f"{name}.csv")
    frame.to_csv(path, index=False)
    _sidecar(name, what, [(c, m) for c, m in columns if c in frame.columns])
    print(f"  -> {path} ({len(frame)} rows)")
    return path


def _read(name, **kw):
    """A table of an earlier section, identifiers as strings."""
    path = out(f"{name}.csv")
    if not path.exists():
        raise FileNotFoundError(f"{path} is missing. Run the section that writes it first.")
    return pd.read_csv(path, dtype={"Case": str, "slide": str}, float_precision="round_trip", **kw)


def _read_percase(cohort):
    """percase_<cohort>.csv, the per-patient table written by s1."""
    return _read(f"percase_{cohort}")


def _p_text(p):
    """p as Supplementary Table 1a prints it: below 0.001 a bound, 3 decimals below 0.1, else 2."""
    p = float(p)
    return "< 0.001" if p < 0.001 else (f"{p:.3f}" if p < 0.1 else f"{p:.2f}")


def _cox(frame, columns):
    """Cox model of overall survival on several exposures, each per 1 SD, adjusted for stage.

    Rows with a missing exposure, follow-up or stage are dropped, and follow-up of zero or less; each
    exposure is z-scored over the rows that remain. One dictionary per exposure, as hr_per_sd.
    """
    x = frame.dropna(subset=list(columns) + ["OS_weeks", "event", "stage"])
    x = x[x.OS_weeks > 0]
    z = [((x[c] - x[c].mean()) / x[c].std()).to_numpy(dtype=float) for c in columns]
    X = np.column_stack(z + [x.stage.to_numpy(dtype=float)])
    b, se = survival.cox_efron(X, x.OS_weeks.to_numpy(dtype=float),
                               x.event.to_numpy(dtype=float).astype(int))
    return [dict(var=c, HR=math.exp(b[i]), lo=math.exp(b[i] - 1.96 * se[i]),
                 hi=math.exp(b[i] + 1.96 * se[i]), p=float(2 * norm.sf(abs(b[i] / se[i]))),
                 logHR=float(b[i]), se=float(se[i]), n=int(len(x)), events=int(x.event.sum()))
            for i, c in enumerate(columns)]


MODEL_COLUMNS = [
    ("var", "the exposure the hazard ratio belongs to"),
    ("HR", "hazard ratio per 1 standard deviation, adjusted for stage"),
    ("lo", "lower 95 % Wald confidence limit"),
    ("hi", "upper 95 % Wald confidence limit"),
    ("p", "Wald p value"),
    ("logHR", "the coefficient on the log hazard scale, which the pooling averages"),
    ("se", "its standard error"),
    ("n", "patients in the model (pooled: summed over the two cohorts)"),
    ("events", "deaths among them (pooled: summed)"),
]


# ---------------------------------------------------------------- s1

def s1_per_patient():
    """Per patient, from the cells and the raw TrxR1 intensity; where the cut sits; panel a densities."""
    cuts, curves = [], []
    for cohort in paths.COHORTS:
        cells = staining.cells_with_intensity(
            cohort, columns=["core", "Case", "tissue", "lineage", "TrxR1"])
        print(f"s1 {cohort}: {len(cells):,} cells, {cells.Case.nunique()} patients")
        fractions = _patient_measures(cells, cohort)
        table = _add_clinical(fractions, cohort)
        _write(table, f"percase_{cohort}",
               f"One row per patient of {cohort} ({len(table)} patients, the order of clinical.csv). "
               "Written by tumor_vs_stroma.py section s1 from cells.tsv, the raw TrxR1 intensity, "
               "clinical.csv, clinical_full.csv and the raw genotype and smoking files. "
               "Patient-level: it stays with the data and is never shared.",
               PERCASE_COLUMNS)
        cuts.append(_cut_position(cells, cohort))
        curves.append(_density(cells, cohort, cuts[-1]))
    _write(pd.DataFrame(cuts), "cut_position",
           "Where the TrxR1 cut sits, per cohort. Written by tumor_vs_stroma.py section s1.",
           CUT_COLUMNS)
    _write(pd.concat(curves, ignore_index=True), "density",
           "The intensity densities of Supplementary Figure 1a. Written by tumor_vs_stroma.py "
           "section s1.",
           [("cohort", "BOMI2 or BOMI1"),
            ("x", "cytoplasmic TrxR1 intensity, the center of one of 1,024 bins from 0 to hi"),
            ("tumor", "density of the intensity of all cells of the tumor compartment"),
            ("stroma", "density of the intensity of all cells of the stromal compartment"),
            ("hi", "the 99th percentile of the intensity of both compartments together; brighter "
                   "cells are left out of both densities")])


PERCASE_COLUMNS = [
    ("Case", "patient identifier, as a string"),
    ("Tumor_type_code", "LUSC or LUAD, empty for other histology (the published column)"),
    ("histology", "LUSC, LUAD or other, from clinical.csv"),
    ("OS_weeks", "overall survival in weeks"),
    ("event", "1 dead, 0 alive at last follow-up"),
    ("stage", "1 to 4 for stage I to IV, the covariate of every survival model"),
    ("smk3", "smoking on three levels: 1 current, 2 former (stopped more than a year before), 3 never; "
             "empty when missing or coded otherwise"),
    ("KEAP1_alt", "Uppsala II: 1 when KEAP1 is altered in the mutation summary, 0 wild type, empty without data"),
    ("NFE2L2_alt", "the same for NFE2L2"),
    ("CUL3_alt", "the same for CUL3"),
    ("NRF2path", "1 when KEAP1, NFE2L2 or CUL3 is altered (the NRF2 pathway), 0 otherwise"),
    ("stromal_all", "TrxR1-positive fraction of all stromal cells with a call; empty below 20 stromal cells"),
    ("n_stroma", "stromal cells of the patient, every one counted"),
    ("cancer", "TrxR1-positive fraction of the cancer cells (cytokeratin lineage in the tumor compartment); "
               "empty below 20 cancer cells"),
    ("n_CK_T", "cancer cells of the patient"),
    ("f_cancer", "TrxR1-positive fraction of the cancer cells with a missing call counted as negative; "
                 "empty below 20 (the slide adjustment)"),
    ("bg_cancer", "background index of the cancer cells: mean intensity of the TrxR1-negative ones; "
                  "empty below 10 of them"),
    ("f_stroma", "the same fraction for all stromal cells"),
    ("bg_stroma", "the background index of all stromal cells"),
    ("f_T", "the same fraction for the stromal T cells"),
    ("bg_T", "the background index of the stromal T cells"),
    ("f_ref", "the same fraction for the reference stromal group: fibroblasts in Uppsala II, stromal "
              "cells other than T cells in Uppsala I"),
    ("bg_ref", "the background index of that group"),
    ("f_mac", "Uppsala II: the same fraction for the stromal macrophages (Uppsala I has no macrophage marker, "
              "empty)"),
    ("bg_mac", "the background index of the stromal macrophages"),
    ("slide", "the TMA slide holding most of the patient's cells"),
    ("first_row", "position of the patient's first cell in cells.tsv; the script reads the Uppsala II "
                  "patients in this order"),
]

CUT_COLUMNS = [
    ("cohort", "BOMI2 or BOMI1"),
    ("cut", "the TrxR1 cut of the cohort, as thresholds.csv records it"),
    ("density_crossing", "where the densities of density.csv cross between their medians; must lie within "
                         "0.1 of the cut"),
    ("pct_cancer_above", "percent of all cancer cells above the cut (cells without an intensity count as below)"),
    ("pct_stroma_above", "percent of all stromal cells above the cut"),
]


def _fraction(sub, floor):
    """Mean call per patient (cells without a call skipped) and cells per patient (all counted)."""
    g = sub.TrxR1.astype(float).groupby(sub.Case)
    n = g.size()
    return g.mean().where(n >= floor), n


def _patient_measures(cells, cohort):
    """Every per-patient measure of the cells: fractions, background, slide, order.

    The fraction of each stromal cell type is computed in s7 (lineage_profile_patients.csv, at least
    20 cells of the type), the definition every printed cell-type value uses.
    """
    stroma = cells[cells.tissue == "Stroma"]
    cancer = cells[(cells.tissue == "Tumor") & (cells.lineage == "CK")]
    t = pd.DataFrame(index=pd.Index(sorted(cells.Case.unique()), name="Case"))
    t["stromal_all"], t["n_stroma"] = _fraction(stroma, MIN_STROMA)
    t["cancer"], t["n_CK_T"] = _fraction(cancer, MIN_CANCER)

    # Fractions and background indices for the slide adjustment.
    reference = (stroma[stroma.lineage == "Fib"] if cohort == "BOMI2"
                 else stroma[stroma.lineage != "T"])
    for key, group in (("cancer", cancer), ("stroma", stroma),
                       ("T", stroma[stroma.lineage == "T"]), ("ref", reference)):
        t[f"f_{key}"], t[f"bg_{key}"] = staining.fraction_and_background(group)
    # The Results say the fractions of fibroblasts, T cells and macrophages correlate with each other
    # (Uppsala II); the macrophage pairs of s4 need the macrophage fraction and background.
    if cohort == "BOMI2":
        t["f_mac"], t["bg_mac"] = staining.fraction_and_background(stroma[stroma.lineage == "Mac"])
    else:
        t["f_mac"], t["bg_mac"] = np.nan, np.nan
    t["slide"] = staining.patient_slide(cells, cohort)
    t["first_row"] = pd.Series(np.arange(len(cells)), index=cells.index).groupby(cells.Case).min()
    return t.reset_index()


def _genotype():
    """Uppsala II genotype per patient from the raw mutation summary (patients as L<Case>)."""
    g = pd.read_excel(paths.raw("Bomi2", MUTATION_FILE), usecols=["Projektnr"] + GENES)
    if g.Projektnr.duplicated().any():
        raise ValueError(f"{MUTATION_FILE}: a patient appears twice")
    out_ = pd.DataFrame({"Case": g.Projektnr.astype(str).str.strip().str.replace("^L", "", regex=True)})
    for gene in GENES:
        v = g[gene].astype(str).str.strip()
        unknown = set(v) - ALTERED - {"WT"}
        if unknown:
            raise ValueError(f"{MUTATION_FILE}: {gene} has entries {sorted(unknown)} that are neither "
                             "wild type nor a known alteration")
        out_[f"{gene}_alt"] = np.where(v == "WT", 0.0, 1.0)
    out_["NRF2path"] = out_[["KEAP1_alt", "NFE2L2_alt", "CUL3_alt"]].max(axis=1)
    return out_


def _smoking_bomi2():
    """Uppsala II smoking on three levels, from the cleared sheet of the Swedish clinical database."""
    k = pd.read_excel(paths.raw("Bomi2", SWEDISH_FILE), sheet_name=SWEDISH_SHEET)
    s = pd.DataFrame({"Case": k[SWEDISH_ID].astype(str).str.replace("L", "", regex=False)
                                           .astype(int).astype(str),
                      "smk3": pd.to_numeric(k[SWEDISH_SMOKING], errors="coerce")})
    if s.Case.duplicated().any():
        raise ValueError(f"{SWEDISH_FILE}: a patient appears twice in {SWEDISH_SHEET}")
    return s


def _check_bomi1_smoking_codes(full):
    """The Uppsala I smoking codes follow the code book of the 1995 to 2005 master sheet.

    The Uppsala I clinical file holds the codes without a legend. The master sheet of the same
    cohort carries the legend (1 current, 2 former, 3 never) and the same codes for the same
    patients; both are checked, so that the three levels are read as the legend says.
    """
    legend = pd.read_excel(paths.raw("Bomi2", MASTER_FILE), sheet_name="Legende")
    if not {"1=current", "2=ex >1year", "3=aldrig"} <= set(legend.Smoking.dropna()):
        raise ValueError(f"{MASTER_FILE}: the smoking legend is not 1 current, 2 former, 3 never")
    master = pd.read_excel(paths.raw("Bomi2", MASTER_FILE), sheet_name="Data",
                           usecols=["Project ID", "Smoking"])
    master["Case"] = master["Project ID"].astype("Int64").astype(str)
    both = full[["Case", "Smoking"]].merge(master[["Case", "Smoking"]], on="Case",
                                           suffixes=("", "_master"))
    if len(both) != len(full) or not (both.Smoking.fillna(-1) == both.Smoking_master.fillna(-1)).all():
        raise ValueError("the Uppsala I smoking codes differ from the master sheet")


def _add_clinical(fractions, cohort):
    """The clinical, genotype and smoking columns, on the patients of clinical.csv in its order."""
    clin = clinicalmod.load_clinical(cohort)
    # Smoking is read for smk3 (Uppsala I) and for the check of the three levels (Uppsala II); the
    # output keeps smk3 only (PERCASE_COLUMNS).
    keep = ["Case", "Tumor_type_code", "Smoking"]
    full = pd.read_csv(paths.intermediate(cohort, "clinical_full.csv"), dtype={"Case": str},
                       usecols=keep)
    t = clin[["Case", "histology", "OS_weeks", "event", "stage"]].merge(full, on="Case", how="left")
    code = t.Tumor_type_code.where(t.Tumor_type_code.isin(["LUSC", "LUAD"]), "other")
    if not (code == t.histology).all():
        raise ValueError(f"{cohort}: clinical.csv and clinical_full.csv disagree on histology")
    if cohort == "BOMI2":
        t = t.merge(_smoking_bomi2(), on="Case", how="left")
        if t.smk3.isna().any() or not (t.Smoking == t.smk3.map({1: "Yes", 2: "Yes", 3: "No"})).all():
            raise ValueError("the three smoking levels do not reproduce the Yes/No column of Uppsala II")
        t = t.merge(_genotype(), on="Case", how="left")
    else:
        _check_bomi1_smoking_codes(full[full.Case.isin(clin.Case)])
        t["smk3"] = t.Smoking.where(t.Smoking.isin([1, 2, 3]))
        t = t.drop(columns=["Smoking"])
    t = t.merge(fractions, on="Case", how="left")
    if t.n_stroma.isna().all():
        raise ValueError(f"merging the cells of {cohort} with clinical.csv kept no patient")
    order = [c for c, _ in PERCASE_COLUMNS if c in t.columns]
    return t[order]


def _cut_position(cells, cohort):
    """Where the cut sits among the cancer cells and the stromal cells of a cohort."""
    cut = staining.trxr1_cut(cohort)
    cancer = cells.loc[(cells.tissue == "Tumor") & (cells.lineage == "CK"), "TrxR1_intensity"]
    stroma = cells.loc[cells.tissue == "Stroma", "TrxR1_intensity"]
    row = dict(cohort=cohort, cut=cut,
               pct_cancer_above=100 * float((cancer > cut).mean()),
               pct_stroma_above=100 * float((stroma > cut).mean()))
    print(f"  {cohort}: cut {cut:.4f}; cancer cells above {row['pct_cancer_above']:.1f} %, "
          f"stromal cells {row['pct_stroma_above']:.1f} %")
    return row


def _nrd0(x):
    """R's default bandwidth (bw.nrd0)."""
    s = min(np.std(x, ddof=1), (np.quantile(x, .75) - np.quantile(x, .25)) / 1.34)
    return 0.9 * s * len(x) ** -0.2


def _kde(x, grid):
    """Binned Gaussian density on a grid, the way R's density bins before it smooths."""
    bw = _nrd0(x)
    h, e = np.histogram(x, bins=len(grid) - 1, range=(grid[0], grid[-1]))
    c = (e[:-1] + e[1:]) / 2
    step = c[1] - c[0]
    k = np.arange(-int(4 * bw / step) - 1, int(4 * bw / step) + 2) * step
    y = np.convolve(h, np.exp(-0.5 * (k / bw) ** 2), mode="same")
    return c, y / (y.sum() * step)


def _density(cells, cohort, cut_row):
    """The two compartment densities of panel a, and the check that they cross at the cut.

    As the cut was set (mif/02_cells.R, crossover_cut): intensities above the 99th percentile of both
    compartments left out, a Gaussian density with R's default bandwidth on a grid from 0 to that
    percentile, the crossing searched between the two medians.
    """
    a = cells.loc[cells.tissue == "Tumor", "TrxR1_intensity"].dropna().to_numpy()
    b = cells.loc[cells.tissue == "Stroma", "TrxR1_intensity"].dropna().to_numpy()
    hi = float(np.quantile(np.concatenate([a, b]), 0.99))
    a, b = a[a <= hi], b[b <= hi]
    x, ya = _kde(a, np.linspace(0, hi, 1025))
    _, yb = _kde(b, np.linspace(0, hi, 1025))
    m1, m2 = sorted([np.median(a), np.median(b)])
    window = np.where((x >= m1) & (x <= m2))[0]
    cross = float(x[window][np.where(np.diff(np.sign(ya[window] - yb[window])) != 0)[0][0]])
    if abs(cross - cut_row["cut"]) >= 0.1:
        raise ValueError(f"{cohort}: the densities cross at {cross:.3f}, not at the cut {cut_row['cut']:.3f}")
    cut_row["density_crossing"] = cross
    return pd.DataFrame({"cohort": cohort, "x": x, "tumor": ya, "stroma": yb, "hi": hi})


# ---------------------------------------------------------------- s2

TEST_COLUMNS = ["family", "cohort", "histology", "variable", "outcome", "test", "n1", "n0", "group1", "group0",
                "median1", "median0", "stat", "p"]
GENOTYPE = ("NRF2path", "NRF2path (KEAP1/NFE2L2/CUL3)")   # tumors with a KEAP1, NFE2L2 or CUL3 alteration


def s2_driver_tests():
    """The genotype tests of Supplementary Figure 1d, and the smoking comparison.

    Uppsala II tumors with a KEAP1, NFE2L2 or CUL3 alteration against the others, within each histology:
    the TrxR1-positive fraction of cancer cells (the positive control of the assay) and of stromal cells,
    two-sided Mann-Whitney test. No random numbers.
    """
    b2 = _read_percase("BOMI2").sort_values("first_row").reset_index(drop=True)
    b1 = _read_percase("BOMI1")
    rows = []

    def two_group(df, grp, outcome, hist, family, label, g1="1", g0="0"):
        """Mann-Whitney test of one fraction between tumors with (1) and without (0) the alteration."""
        x = df.dropna(subset=[grp, outcome])
        a = x.loc[x[grp] == 1, outcome].to_numpy()
        b = x.loc[x[grp] == 0, outcome].to_numpy()
        U, p = ranktests.mw_test(a, b)
        return dict(family=family, cohort="BOMI2", histology=hist, variable=label,
                    outcome=outcome, test="Mann-Whitney", n1=len(a), n0=len(b), group1=g1, group0=g0,
                    median1=np.median(a) if len(a) else np.nan,
                    median0=np.median(b) if len(b) else np.nan, stat=U, p=p)

    grp, label = GENOTYPE
    for outcome, family in (("cancer", "control"), ("stromal_all", "B")):
        for h in HISTS:
            rows.append(two_group(b2[b2.Tumor_type_code == h], grp, outcome, h, family, label))
    table = pd.DataFrame(rows)[TEST_COLUMNS]
    _write(table, "imaging_driver_tests",
           "Uppsala II tumors with a KEAP1, NFE2L2 or CUL3 alteration against the others: the TrxR1-positive "
           "fraction of cancer cells (positive control) and of stromal cells, Mann-Whitney p within each "
           "histology. The p values of Supplementary Figure 1d. Written by tumor_vs_stroma.py section s2.",
           [("family", "control for the positive control (the cancer-cell fraction), B for the stromal fraction"),
            ("cohort", "BOMI2"),
            ("histology", "LUAD or LUSC"),
            ("variable", "the genotype: NRF2path (KEAP1, NFE2L2 or CUL3 altered)"),
            ("outcome", "the TrxR1-positive fraction tested: cancer or stromal_all (columns of "
                        "percase_BOMI2.csv)"),
            ("test", "Mann-Whitney"),
            ("n1", "tumors with the alteration"),
            ("n0", "tumors without it"),
            ("group1", "what group 1 is (1, altered)"),
            ("group0", "what group 0 is (0, not altered)"),
            ("median1", "median fraction in group 1"),
            ("median0", "median fraction in group 0"),
            ("stat", "Mann-Whitney U"),
            ("p", "two-sided Mann-Whitney p")])
    _smoking_detail(b1, _read_percase("BOMI2"))


def _smoking_detail(b1, b2):
    """Current against former smokers, the TrxR1-positive fraction of stromal cells, per cohort and histology.

    The raw clinical sheets hold smoking on three levels in both cohorts. Current and former smokers are
    compared within each cohort and histology with the two-sided Mann-Whitney test, the comparison behind
    the smoking sentence of the Results. Both cohorts are read in the order of clinical.csv. No random numbers.
    """
    cl = {}
    for cohort, c in (("BOMI1", b1), ("BOMI2", b2)):
        c = c[c.Tumor_type_code.isin(HISTS) & c.smk3.notna()].copy()
        c["current"] = c.smk3.map({1: 1, 2: 0})
        print(f"  {cohort}: " + ", ".join(
            f"{h} current {int((g.smk3 == 1).sum())} former {int((g.smk3 == 2).sum())} "
            f"never {int((g.smk3 == 3).sum())}" for h, g in c.groupby("Tumor_type_code")))
        cl[cohort] = c
    outcome = "stromal_all"
    rows = []
    for cohort, c in cl.items():
        for h in HISTS:
            x = c[c.Tumor_type_code == h].dropna(subset=["current", outcome])
            a = x.loc[x.current == 1, outcome].to_numpy()
            b = x.loc[x.current == 0, outcome].to_numpy()
            U, p = ranktests.mw_test(a, b)
            rows.append(dict(test_id="B", cohort=cohort, histology=h, comparison="current vs former",
                             outcome=outcome, test="Mann-Whitney", n1=len(a), n0=len(b),
                             median1=np.median(a), median0=np.median(b), p=p))
    _write(pd.DataFrame(rows), "smoking_detail",
           "Current against former smokers, the TrxR1-positive fraction of stromal cells, per cohort and "
           "histology. Backs the sentence that smoking was not associated with the stromal fraction in "
           "Uppsala I and II. Written by tumor_vs_stroma.py section s2.",
           [("test_id", "B, current against former smokers within one cohort and histology"),
            ("cohort", "BOMI2 or BOMI1"),
            ("histology", "LUAD or LUSC"),
            ("comparison", "group 1 against group 0"),
            ("outcome", "the TrxR1-positive fraction compared: stromal_all (column of percase_<cohort>.csv)"),
            ("test", "Mann-Whitney"),
            ("n1", "current smokers"), ("n0", "former smokers (stopped more than a year before)"),
            ("median1", "median fraction of the current smokers"),
            ("median0", "median fraction of the former smokers"),
            ("p", "two-sided Mann-Whitney p")])


# ---------------------------------------------------------------- s3

def _patients():
    """Both cohorts, LUSC and LUAD, with the two fractions and the clinical columns."""
    frames = []
    for cohort in paths.COHORTS:
        d = _read_percase(cohort)
        keep = ["Case", "histology", "OS_weeks", "event", "stage", "stromal_all", "n_stroma",
                "cancer", "n_CK_T"] + (["NRF2path"] if cohort == "BOMI2" else [])
        frames.append(d[keep].assign(cohort=cohort))
    D = pd.concat(frames, ignore_index=True)
    return D[D.histology.isin(["LUSC", "LUAD"])].reset_index(drop=True)


def s3_cancer_and_stroma():
    """Distribution, survival, joint models and bulk RNA of the cancer-cell and stromal fractions."""
    D = _patients()

    # Distribution per cohort and histology.
    rows = []
    for (cohort, h), g in D.groupby(["cohort", "histology"]):
        for m in ("cancer", "stromal_all"):
            v = g[m].dropna()
            rows.append(dict(cohort=cohort, histology=h, measure=m, n=len(v), median=v.median(),
                             q1=v.quantile(.25), q3=v.quantile(.75)))
    T1 = pd.DataFrame(rows)
    _write(T1, "distribution",
           "The TrxR1-positive fractions of cancer cells and of all stromal cells, per cohort and "
           "histology, patients with at least 20 cells of the kind. Written by tumor_vs_stroma.py "
           "section s3.",
           [("cohort", "BOMI2 or BOMI1"), ("histology", "LUSC or LUAD"),
            ("measure", "cancer or stromal_all"), ("n", "patients"),
            ("median", "median fraction"), ("q1", "25th percentile"), ("q3", "75th percentile")])

    # The cancer-cell fraction against overall survival.
    rows, per = [], {}
    for h in ("LUSC", "LUAD"):
        for cohort in paths.COHORTS:
            x = D[(D.cohort == cohort) & (D.histology == h)]
            r = survival.hr_per_sd(x, "cancer")
            per[cohort] = r
            rows.append(dict(histology=h, cohort=cohort, model="stage", var="cancer", **r))
        rows.append(dict(histology=h, cohort="pooled", model="stage, fixed effect", var="cancer",
                         **survival.meta_fixed([per["BOMI2"], per["BOMI1"]])))
    _write(pd.DataFrame(rows)[["histology", "cohort", "model"] + [c for c, _ in MODEL_COLUMNS]],
           "cancer_survival",
           "The TrxR1-positive fraction of cancer cells against overall survival, Cox model adjusted "
           "for stage, per cohort and histology and pooled by fixed effect. Written by "
           "tumor_vs_stroma.py section s3.",
           [("histology", "LUSC or LUAD"), ("cohort", "BOMI2, BOMI1 or pooled"),
            ("model", "stage, or stage, fixed effect for the pooled row")] + MODEL_COLUMNS)

    # Stromal and cancer-cell fractions in one model, on the patients who carry both.
    rows, per = [], {}
    for h in ("LUSC", "LUAD"):
        for cohort in paths.COHORTS:
            x = D[(D.cohort == cohort) & (D.histology == h)].dropna(subset=["cancer", "stromal_all"])
            alone = _cox(x, ["stromal_all"])[0]
            joint = _cox(x, ["stromal_all", "cancer"])
            per[(cohort, "alone")], per[(cohort, "joint")], per[(cohort, "cancer")] = alone, joint[0], joint[1]
            rows.append(dict(histology=h, cohort=cohort, model="stromal alone (same patients)", **alone))
            for r in joint:
                rows.append(dict(histology=h, cohort=cohort, model="stromal + cancer", **r))
        for key, label, var in (("alone", "stromal alone, fixed effect", "stromal_all"),
                                ("joint", "stromal joint, fixed effect", "stromal_all"),
                                ("cancer", "cancer joint, fixed effect", "cancer")):
            rows.append(dict(histology=h, cohort="pooled", model=label, var=var,
                             **survival.meta_fixed([per[("BOMI2", key)], per[("BOMI1", key)]])))
    T3 = pd.DataFrame(rows)[["histology", "cohort", "model"] + [c for c, _ in MODEL_COLUMNS]]
    _write(T3, "joint_models",
           "The stromal TrxR1-positive fraction alone and together with the cancer-cell fraction, Cox "
           "models adjusted for stage, on the patients who carry both fractions; each histology pooled "
           "over the two cohorts by fixed effect. Written by tumor_vs_stroma.py section s3.",
           [("histology", "LUSC or LUAD"), ("cohort", "BOMI2, BOMI1 or pooled"),
            ("model", "stromal alone (same patients); stromal + cancer, one row per term; the pooled rows "
                      "of the stromal term alone, the stromal term of the joint model and the "
                      "cancer-cell term of the joint model")] + MODEL_COLUMNS)
    _check_against_figure1(T3)
    _bulk(D)


def _check_against_figure1(T3):
    """The stromal fraction alone reproduces Figure 1b, when chapter1.py has written it.

    Every patient with at least 20 stromal cells in these histologies also has 20 cancer cells, so
    the two patient sets are the same and the hazard ratios must agree. A difference stops the run.
    """
    path = paths.WORK_ROOT / "output" / "chapter1" / "fig1b_models.csv"
    if not path.exists():
        print("  (chapter1/fig1b_models.csv not written yet; the Figure 1b check is skipped)")
        return
    f1 = pd.read_csv(path)
    for r in T3[T3.model.isin(["stromal alone (same patients)", "stromal alone, fixed effect"])].itertuples():
        ref = f1[(f1.histology == r.histology) & (f1.cohort == r.cohort)].iloc[0]
        if abs(ref.HR - r.HR) > 1e-9 or int(ref.n) != int(r.n):
            raise ValueError(f"{r.histology} {r.cohort}: the stromal fraction alone gives HR {r.HR} on "
                             f"{r.n} patients, Figure 1b {ref.HR} on {ref.n}")
    print("  the stromal fraction alone reproduces the six Figure 1b models")


def _bulk(D):
    """Bulk tumor TXNRD1 (GSE81089) against the two fractions, the genotype and survival."""
    named = paths.CFG.get("gse81089")
    if not named:
        print("  no gse81089 entry in the configuration; bulk_rna.csv is not written")
        return
    B = pd.read_csv(named, usecols=["Case", TXNRD1_ID])
    B["Case"] = B.Case.astype(str)
    B["bulk_TXNRD1"] = np.log2(B[TXNRD1_ID] + 1)
    X = D[D.cohort == "BOMI2"].merge(B[["Case", "bulk_TXNRD1"]], on="Case")
    if X.empty:
        raise ValueError(f"merging {named} with percase_BOMI2.csv kept no patient")
    _write(X[["Case", "histology", "cancer", "stromal_all", "NRF2path", "bulk_TXNRD1"]],
           "bulk_per_patient",
           "Uppsala II patients with bulk RNA sequencing (GSE81089), the points of Supplementary "
           "Figure 1f. Written by tumor_vs_stroma.py section s3. Patient-level: never shared.",
           [("Case", "patient identifier, as a string"), ("histology", "LUSC or LUAD"),
            ("cancer", "TrxR1-positive fraction of cancer cells"),
            ("stromal_all", "TrxR1-positive fraction of stromal cells"),
            ("NRF2path", "1 when KEAP1, NFE2L2 or CUL3 is altered, 0 otherwise"),
            ("bulk_TXNRD1", "log2(TXNRD1 expression + 1) in the deposited table")])
    rows = []
    for h, g in X.groupby("histology"):
        for m in ("cancer", "stromal_all"):
            gg = g.dropna(subset=[m])
            r, p = spearmanr(gg.bulk_TXNRD1, gg[m])
            rows.append(dict(histology=h, test=f"Spearman bulk TXNRD1 vs {m}", estimate=r, p=p,
                             n=len(gg), n1=np.nan, n0=np.nan))
        gm = g.dropna(subset=["NRF2path"])
        a = gm[gm.NRF2path == 1].bulk_TXNRD1
        b = gm[gm.NRF2path == 0].bulk_TXNRD1
        rows.append(dict(histology=h,
                         test="bulk TXNRD1, NRF2-pathway mutant minus wild type (median, log2)",
                         estimate=a.median() - b.median(),
                         p=mannwhitneyu(a, b).pvalue if len(a) and len(b) else np.nan,
                         n=f"{len(a)} vs {len(b)}", n1=len(a), n0=len(b)))
        r = survival.hr_per_sd(g, "bulk_TXNRD1")
        rows.append(dict(histology=h, test="bulk TXNRD1 and OS, HR per 1 SD, stage-adjusted",
                         estimate=r["HR"], lo=r["lo"], hi=r["hi"], p=r["p"], n=r["n"],
                         events=r["events"], n1=np.nan, n0=np.nan))
    _write(pd.DataFrame(rows)[["histology", "test", "estimate", "p", "n", "lo", "hi", "events", "n1", "n0"]],
           "bulk_rna",
           "Bulk tumor TXNRD1 (GSE81089, Uppsala II) against the TrxR1-positive fractions, the genotype "
           "and overall survival. Written by tumor_vs_stroma.py section s3.",
           [("histology", "LUSC or LUAD"), ("test", "what is estimated"),
            ("estimate", "Spearman rho, the difference in medians (log2), or the hazard ratio per 1 SD "
                         "adjusted for stage"),
            ("p", "Spearman, Mann-Whitney or Wald p"),
            ("n", "patients; for the genotype, altered vs other"),
            ("lo", "lower 95 % confidence limit of the hazard ratio"), ("hi", "upper limit"),
            ("events", "deaths in the survival model"),
            ("n1", "tumors with an alteration"), ("n0", "tumors without one")])


# ---------------------------------------------------------------- s4

def s4_compartments():
    """Cancer-cell against stromal fraction, and T cells against the other stromal group, across patients.

    Three estimates per cohort and histology: the raw Spearman rho; within slide (the ranks of both
    fractions regressed on slide indicators, Pearson correlation of the residuals); and within slide
    with the ranks of both groups' background indices as further covariates. Patients enter with
    both fractions (20 cells), both background indices (10 negative cells) and a slide.
    """
    rows, points = [], []
    for cohort in paths.COHORTS:
        X = _read_percase(cohort)
        X = X[X.histology.isin(["LUSC", "LUAD"])]
        for h in ("LUSC", "LUAD"):
            pairs = (("cancer", "stroma", "cancer cells with stromal cells"),
                     ("T", "ref", "T cells with fibroblasts" if cohort == "BOMI2"
                      else "T cells with other stromal cells"))
            if cohort == "BOMI2":
                pairs += (("mac", "T", "macrophages with T cells"), ("mac", "ref", "macrophages with fibroblasts"))
            for a, b, label in pairs:
                x = X[X.histology == h].dropna(subset=[f"f_{a}", f"f_{b}", f"bg_{a}", f"bg_{b}", "slide"])
                S = staining.slide_indicators(x.slide)
                fa, fb = x[f"f_{a}"].to_numpy(), x[f"f_{b}"].to_numpy()
                r0, p0 = spearmanr(fa, fb)
                r1, p1 = staining.partial_rank_correlation(fa, fb, S)
                r2, p2, ea, eb = staining.partial_rank_correlation(
                    fa, fb, S + [rankdata(x[f"bg_{a}"]), rankdata(x[f"bg_{b}"])], return_residuals=True)
                rows.append(dict(cohort=cohort, histology=h, pair=label, n=len(x),
                                 n_slides=x.slide.nunique(), rho_raw=r0, p_raw=p0,
                                 rho_within_slide=r1, p_within_slide=p1,
                                 rho_slide_background=r2, p_slide_background=p2))
                if a == "cancer":
                    n = len(x)
                    points.append(pd.DataFrame({
                        "cohort": cohort, "Case": x.Case.to_numpy(), "histology": h,
                        "slide": x.slide.to_numpy(), "f_cancer": fa, "f_stroma": fb,
                        "bg_cancer": x.bg_cancer.to_numpy(), "bg_stroma": x.bg_stroma.to_numpy(),
                        "cancer_adjusted": ea / n, "stroma_adjusted": eb / n}))
    _write(pd.DataFrame(rows), "compartment_correlation",
           "How closely the TrxR1-positive fractions of two cell groups follow each other across "
           "patients, per cohort and histology. Written by tumor_vs_stroma.py section s4.",
           [("cohort", "BOMI2 or BOMI1"), ("histology", "LUSC or LUAD"),
            ("pair", "the two groups: cancer cells with all stromal cells; stromal T cells with "
                     "fibroblasts (Uppsala I: with the other stromal cells); Uppsala II also macrophages "
                     "with T cells and macrophages with fibroblasts"),
            ("n", "patients with both fractions, both background indices and a slide"),
            ("n_slides", "TMA slides among them"),
            ("rho_raw", "Spearman rho"), ("p_raw", "its p"),
            ("rho_within_slide", "correlation of the ranks after slide indicators are regressed out"),
            ("p_within_slide", "t-test p, n - k - 2 degrees of freedom"),
            ("rho_slide_background", "the same with the ranks of both background indices added"),
            ("p_slide_background", "its p")])
    _write(pd.concat(points, ignore_index=True), "compartment_points",
           "The points of Supplementary Figure 1e: per patient the two fractions as measured and the "
           "rank residuals after slide and background. Written by tumor_vs_stroma.py section s4. "
           "Patient-level: never shared.",
           [("cohort", "BOMI2 or BOMI1"), ("Case", "patient identifier, as a string"),
            ("histology", "LUSC or LUAD"), ("slide", "the patient's TMA slide"),
            ("f_cancer", "TrxR1-positive fraction of cancer cells, missing calls as negative"),
            ("f_stroma", "the same for stromal cells"),
            ("bg_cancer", "background index of the cancer cells"),
            ("bg_stroma", "background index of the stromal cells"),
            ("cancer_adjusted", "rank of f_cancer within the histology, slide and both background "
                                "ranks regressed out, divided by the patients"),
            ("stroma_adjusted", "the same for f_stroma")])


# ---------------------------------------------------------------- s5

SUPP1A_COLUMNS = [
    ("section", "the block of the table"),
    ("row", "short name of the row, used by manifest/check.py"),
    ("label", "the row as the table names it"),
    ("histology", "LUSC or LUAD"),
    ("cohort", "BOMI2, BOMI1 or pooled"),
    ("value", "the value the cell prints: a fraction, a hazard ratio, rho, a difference in medians, "
              "a ratio of medians, or a count (corr_patients: the fewest patients, the lower end of a "
              "printed range)"),
    ("lo", "lower limit: 25th percentile or lower 95 % confidence limit; for corr_patients the fewest "
           "patients over the rows of the correlation block"),
    ("hi", "upper limit: 75th percentile or upper 95 % confidence limit; for corr_patients the most "
           "patients over the rows of the correlation block"),
    ("p", "p value of the cell"),
    ("n", "patients or tumors (corr_patients: the most patients over the rows of its block)"),
    ("n0", "bulk RNA genotype row: tumors without an alteration"),
    ("events", "deaths"),
    ("shown", "the cell as the table prints it"),
]


def s5_supptable1a():
    """Supplementary Table 1a, one row per printed cell, from the tables of s3 and s4; legend counts.

    supptable1a.csv holds every block of the table except the survival block by cell type with its
    patients (deaths) row (s8, supptable1a_celltypes.csv) and the covariate block (s9,
    supptable1a_covariates.csv). The table does not repeat what a figure panel shows: the genotype
    block is Supplementary Figure 1d (imaging_driver_tests.csv), the bulk RNA correlations
    Supplementary Figure 1f and the correlations of cancer cells with stromal cells Supplementary
    Figure 1e. The patients row of the correlation block gives the patients of the printed rows
    (the correlations between stromal cell types), one number when they agree, else "a to b".
    suppfig1_numbers.csv holds the counts of the Supplementary Figure 1 legend and the bulk RNA count
    drawn in Figure 1a.
    """
    T1 = _read("distribution")
    T2 = _read("cancer_survival")
    T3 = _read("joint_models")
    corr = _read("compartment_correlation")
    rows = []

    def add(section, row, label, h, cohort, shown, **values):
        """Append one printed cell of the table, with the numbers behind it."""
        rows.append(dict(section=section, row=row, label=label, histology=h, cohort=cohort,
                         shown=shown, **values))

    def hr_text(r):
        """A hazard ratio as the table prints it: HR (lower to upper), p."""
        return f"{r.HR:.2f} ({r.lo:.2f} to {r.hi:.2f}), p {_p_text(r.p)}"

    dist = "TrxR1-positive cells per patient, median (IQR)"
    for h in ("LUSC", "LUAD"):
        for cohort in paths.COHORTS:
            for m, row, label in (("cancer", "dist_cancer", "cancer cells"),
                                  ("stromal_all", "dist_stroma", "stromal cells")):
                r = T1[(T1.cohort == cohort) & (T1.histology == h) & (T1.measure == m)].iloc[0]
                add(dist, row, label, h, cohort,
                    f"{100 * r['median']:.0f} % ({100 * r.q1:.0f} to {100 * r.q3:.0f} %)",
                    value=r["median"], lo=r.q1, hi=r.q3, n=r.n)
            r = T1[(T1.cohort == cohort) & (T1.histology == h) & (T1.measure == "cancer")].iloc[0]
            add(dist, "dist_patients", "patients", h, cohort, f"{int(r.n)}", value=r.n, n=r.n)

    os_ = "Overall survival, hazard ratio per 1 SD (95 % CI), adjusted for stage"
    for h in ("LUSC", "LUAD"):
        for cohort in ("BOMI2", "BOMI1", "pooled"):
            sel = [("os_cancer", "cancer cells", T2[(T2.histology == h) & (T2.cohort == cohort)]),
                   ("os_stroma", "stromal cells (Figure 1b)",
                    T3[(T3.histology == h) & (T3.cohort == cohort)
                       & T3.model.isin(["stromal alone (same patients)", "stromal alone, fixed effect"])]),
                   ("os_stroma_joint", "stromal cells, with cancer cells in the model",
                    T3[(T3.histology == h) & (T3.cohort == cohort) & (T3["var"] == "stromal_all")
                       & T3.model.isin(["stromal + cancer", "stromal joint, fixed effect"])]),
                   ("os_cancer_joint", "cancer cells, with stromal cells in the model",
                    T3[(T3.histology == h) & (T3.cohort == cohort) & (T3["var"] == "cancer")
                       & T3.model.isin(["stromal + cancer", "cancer joint, fixed effect"])])]
            for row, label, frame in sel:
                if len(frame) != 1:
                    raise ValueError(f"{row} {h} {cohort}: {len(frame)} model rows, expected one")
                r = frame.iloc[0]
                add(os_, row, label, h, cohort, hr_text(r), value=r.HR, lo=r.lo, hi=r.hi, p=r.p,
                    n=r.n, events=r.events)
            r = T2[(T2.histology == h) & (T2.cohort == cohort)].iloc[0]
            add(os_, "os_patients", "patients (deaths)", h, cohort, f"{int(r.n)} ({int(r.events)})",
                value=r.n, n=r.n, events=r.events)

    if out("bulk_rna.csv").exists():
        T4 = _read("bulk_rna")
        bulk = "Bulk tumor TXNRD1 (GSE81089, log2)"
        for h in ("LUSC", "LUAD"):
            g = T4[T4.histology == h].set_index("test")
            r = g.loc["bulk TXNRD1, NRF2-pathway mutant minus wild type (median, log2)"]
            add(bulk, "bulk_geno_diff", "altered against other tumors, difference of medians", h, "BOMI2",
                f"{r.estimate:+.1f}, p {_p_text(r.p)}", value=r.estimate, p=r.p, n=r.n1, n0=r.n0)
            add(bulk, "bulk_geno_fold", "the same as a ratio of medians, 2 to the power of the difference "
                "(the text: four to five times higher)", h, "BOMI2", f"{2 ** r.estimate:.1f}",
                value=2 ** r.estimate)
            r = g.loc["bulk TXNRD1 and OS, HR per 1 SD, stage-adjusted"]
            add(bulk, "bulk_os", "overall survival, hazard ratio per 1 SD", h, "BOMI2",
                f"{r.estimate:.2f} ({r.lo:.2f} to {r.hi:.2f}), p {_p_text(r.p)}",
                value=r.estimate, lo=r.lo, hi=r.hi, p=r.p, n=float(r.n), events=r.events)
            n = g.loc["Spearman bulk TXNRD1 vs cancer"].n
            add(bulk, "bulk_tumors", "tumors", h, "BOMI2", f"{int(float(n))}", value=float(n), n=float(n))

    cc = "TrxR1-positive fractions correlated across patients, Spearman rho"
    for h in ("LUSC", "LUAD"):
        for cohort in paths.COHORTS:
            g = corr[(corr.cohort == cohort) & (corr.histology == h)]
            kept = []                                   # patients of the printed rows, for the patients row
            t = g[g.pair.isin(["T cells with fibroblasts", "T cells with other stromal cells"])].iloc[0]
            add(cc, "corr_T_ref", "T cells with fibroblasts (Uppsala I: with other stromal cells), within "
                "slide and adjusted for background staining", h, cohort,
                f"{t.rho_slide_background:.2f}, p {_p_text(t.p_slide_background)}",
                value=t.rho_slide_background, p=t.p_slide_background, n=t.n)
            kept.append(int(t.n))
            for pair, row in (("macrophages with T cells", "corr_mac_T"),
                              ("macrophages with fibroblasts", "corr_mac_ref")):
                m = g[g.pair == pair]
                if len(m):
                    m = m.iloc[0]
                    add(cc, row, f"{pair[0].upper()}{pair[1:]}, within slide and adjusted for background staining",
                        h, cohort, f"{m.rho_slide_background:.2f}, p {_p_text(m.p_slide_background)}",
                        value=m.rho_slide_background, p=m.p_slide_background, n=m.n)
                    kept.append(int(m.n))
            lo, hi = min(kept), max(kept)
            add(cc, "corr_patients", "patients", h, cohort, f"{lo}" if lo == hi else f"{lo} to {hi}",
                value=lo, lo=lo, hi=hi, n=hi)
    table = pd.DataFrame(rows)
    table = table[[c for c, _ in SUPP1A_COLUMNS if c in table.columns]]
    _write(table, "supptable1a",
           "Supplementary Table 1a, TrxR1 in cancer cells and in stromal cells: one row per printed "
           "cell, with the numbers behind it (bulk_geno_fold is not printed; the text states it as "
           "four to five times higher). Written by tumor_vs_stroma.py section s5 from distribution.csv, "
           "cancer_survival.csv, joint_models.csv, bulk_rna.csv and compartment_correlation.csv. Two "
           "parts of the printed table are in other files of this folder: the survival block by cell "
           "type, with its patients (deaths) row, in supptable1a_celltypes.csv (section s8), and the "
           "covariate block in supptable1a_covariates.csv (section s9). The patients row of the correlation block (corr_patients) gives "
           "the patients of its printed rows, the correlations between stromal cell types: one number "
           "when they agree, else the fewest to the most.",
           SUPP1A_COLUMNS)
    counts = []
    for cohort in paths.COHORTS:
        t = _read_percase(cohort)
        t = t[t.histology.isin(["LUSC", "LUAD"])]
        for m in ("cancer", "stromal_all"):
            counts.append(dict(key="suppfig1b_patients", cohort=cohort, measure=m,
                               value=int(t[m].notna().sum())))
    t = _read_percase("BOMI2")                    # the key of panel d: altered tumors per gene, cancer-cell plot
    t = t[t.histology.isin(["LUSC", "LUAD"]) & t.NRF2path.notna() & (t.n_CK_T >= 20)]
    for g in ("KEAP1", "NFE2L2", "CUL3"):
        counts.append(dict(key="suppfig1d_altered", cohort="BOMI2", measure=g, value=int((t[f"{g}_alt"] == 1).sum())))
    if out("bulk_rna.csv").exists():
        # Figure 1a, "bulk RNA (147)": the tumors of the bulk RNA rows of Supplementary Table 1a, both histologies.
        T4 = _read("bulk_rna")
        n = T4[T4.test == "Spearman bulk TXNRD1 vs cancer"].n.astype(float)
        counts.append(dict(key="fig1a_bulk_tumors", cohort="BOMI2", measure="LUSC and LUAD", value=int(n.sum())))
    _write(pd.DataFrame(counts), "suppfig1_numbers",
           "Numbers of the Supplementary Figure 1 legend, and the bulk RNA count drawn in Figure 1a, that no "
           "other table holds. Written by tumor_vs_stroma.py section s5.",
           [("key", "suppfig1b_patients: patients drawn in panel b (LUSC and LUAD, at least 20 cells "
                    "of the kind); suppfig1d_altered: tumors with that gene altered in the cancer-cell "
                    "plot of panel d (LUSC and LUAD, at least 20 cancer cells); fig1a_bulk_tumors: Uppsala "
                    "II tumors with bulk RNA and a cancer-cell fraction, LUSC and LUAD together (the "
                    "tumors row of the bulk block of Supplementary Table 1a, summed)"),
            ("cohort", "BOMI2 or BOMI1"),
            ("measure", "cancer or stromal_all (panel b); the gene (panel d); LUSC and LUAD (Figure 1a)"),
            ("value", "the number")])


# ---------------------------------------------------------------- s6

def s6_panels():
    """The six panels of Supplementary Figure 1, each drawn from the tables of s1 to s5."""
    plotting.style()
    rng = np.random.default_rng(SEED_JITTER)      # consumed by panel b, then panel d
    _panel_a()
    _panel_b(rng)
    _panel_c()
    _panel_d(rng)
    _panel_e()
    _panel_f()


def _ax_style(ax):
    """The shared axis look of the panels: no top or right spine, gray left and bottom, small ticks."""
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(AXIS)
    ax.tick_params(labelsize=FS, length=1.6)


def _letter(fig, letter):
    """The bold panel letter at the top left of the figure."""
    fig.text(0.0, 1.0, letter, fontsize=9.5, fontweight="bold", va="bottom", ha="left")


def _panel_a():
    """Density of the TrxR1 intensity in the tumor and the stromal compartment, with the cut."""
    import matplotlib.pyplot as plt
    dens, cuts = _read("density"), _read("cut_position").set_index("cohort")
    fig = plt.figure(figsize=(3.6, 1.75))
    axes = [fig.add_axes([0.02, 0.22, 0.46, 0.60]), fig.add_axes([0.53, 0.22, 0.46, 0.60])]
    for ax, cohort in zip(axes, paths.COHORTS):
        d = dens[dens.cohort == cohort]
        x, ya, yb, hi = d.x.to_numpy(), d.tumor.to_numpy(), d.stroma.to_numpy(), float(d.hi.iloc[0])
        cut = float(cuts.loc[cohort, "cut"])
        top = max(ya.max(), yb.max())
        ax.fill_between(x, yb, color=STROMA, alpha=0.18, lw=0)
        ax.plot(x, yb, color=STROMA, lw=1.0, label="stromal compartment")
        ax.fill_between(x, ya, color=TUM, alpha=0.18, lw=0)
        ax.plot(x, ya, color=TUM, lw=1.0, label="tumor compartment")
        ax.axvline(cut, color=plotting.INK, lw=0.8, ls=(0, (3, 2)))
        ax.set_ylim(0, top * 1.55)
        ax.set_xlim(0, hi)
        ax.text(cut - hi * 0.02, top * 1.5, f"cut {cut:.2f}", fontsize=FS, va="top", ha="right")
        ax.text(1.0, 0.97, f"cancer cells {cuts.loc[cohort, 'pct_cancer_above']:.0f} %\n"
                           f"stromal cells {cuts.loc[cohort, 'pct_stroma_above']:.0f} %",
                transform=ax.transAxes, ha="right", va="top", fontsize=FS, color=plotting.INK2,
                linespacing=1.2)
        ax.set_title(plotting.COHORT_LABEL[cohort], fontsize=6.5, pad=2)
        ax.set_yticks([])
        ax.spines[["top", "right", "left"]].set_visible(False)
        ax.spines["bottom"].set_color(AXIS)
        ax.tick_params(labelsize=FS, length=1.6)
        ax.set_xlabel("cytoplasmic TrxR1 intensity", fontsize=FS, labelpad=1)
    axes[0].legend(loc="upper left", bbox_to_anchor=(0.0, 1.30), ncol=2, frameon=False, fontsize=FS,
                   handlelength=1.1, columnspacing=0.7, handletextpad=0.35)
    _letter(fig, "a")
    save(fig, "suppfig1a_density")


def _dots(ax, xp, v, color, rng, alpha=0.45, spread=0.28, half=0.36):
    """One column of jittered dots with the median (horizontal) and the interquartile range."""
    ax.scatter(xp + rng.uniform(-spread, spread, len(v)), v, s=3, color=color, alpha=alpha, lw=0)
    q1, med, q3 = np.quantile(v, [.25, .5, .75])
    ax.plot([xp - half, xp + half], [med, med], color=plotting.INK, lw=1.4)
    ax.plot([xp, xp], [q1, q3], color=plotting.INK, lw=0.8)


def _percent_axis(ax, label):
    """A y axis for a fraction, 0 to 100 percent, with the given label."""
    ax.set_ylim(0, 1.12)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_yticklabels(["0", "25", "50", "75", "100"])
    ax.set_ylabel(label, fontsize=FS, labelpad=1)


def _panel_b(rng):
    """TrxR1-positive fraction of cancer cells and of stromal cells, one dot per patient."""
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(1.9, 1.75))
    ax = fig.add_axes([0.25, 0.22, 0.72, 0.60])
    positions = {("BOMI2", "cancer"): 0, ("BOMI2", "stromal_all"): 1,
                 ("BOMI1", "cancer"): 2.6, ("BOMI1", "stromal_all"): 3.6}
    tables = {c: _read_percase(c) for c in paths.COHORTS}
    for (cohort, m), xp in positions.items():
        t = tables[cohort]
        v = t.loc[t.histology.isin(["LUSC", "LUAD"]), m].dropna().to_numpy()
        _dots(ax, xp, v, TUM if m == "cancer" else STROMA, rng)
    ax.set_xticks(list(positions.values()))
    ax.set_xticklabels(["cancer\ncells", "stromal\ncells"] * 2, fontsize=FS)
    ax.text(0.5, 1.07, plotting.COHORT_LABEL["BOMI2"], ha="center", fontsize=6.3)
    ax.text(3.1, 1.07, plotting.COHORT_LABEL["BOMI1"], ha="center", fontsize=6.3)
    _percent_axis(ax, "TrxR1-positive cells,\n% per patient")
    ax.set_xlim(-0.6, 4.2)
    _ax_style(ax)
    _letter(fig, "b")
    save(fig, "suppfig1b_fractions")


def _panel_c():
    """Hazard ratio of the cancer-cell fraction per cohort, pooled, and pooled beside the stromal one."""
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    t = _read("supptable1a")
    t = t[t.row.isin(["os_cancer", "os_cancer_joint", "os_stroma"])]
    check = t[(t.row == "os_stroma") & (t.histology == "LUSC") & (t.cohort == "pooled")].iloc[0]
    print(f"    (Figure 1b, squamous both cohorts, for reference: HR {check.value:.3f})")
    rows = [(0.0, plotting.COHORT_LABEL["BOMI2"], "os_cancer", "BOMI2"),
            (1.0, plotting.COHORT_LABEL["BOMI1"], "os_cancer", "BOMI1"),
            (2.0, "combined", "os_cancer", "pooled"),
            (3.6, "combined,\nadjusted for\nstromal cells", "os_cancer_joint", "pooled")]
    fig = plt.figure(figsize=(2.3, 1.75))
    ax = fig.add_axes([0.42, 0.22, 0.55, 0.60])
    every = []
    for y, _, row, cohort in rows:
        for dy, h in ((-0.16, "LUSC"), (0.16, "LUAD")):
            r = t[(t.row == row) & (t.cohort == cohort) & (t.histology == h)].iloc[0]
            every += [r.lo, r.hi]
            ax.plot([r.lo, r.hi], [y + dy] * 2, color=TUM, lw=1.0)
            ax.scatter([r.value], [y + dy], s=14, marker=plotting.HM[h],
                       color=TUM if h == "LUSC" else "white", edgecolor=TUM, lw=0.8, zorder=3)
    ax.axvline(1, color=AXIS, lw=0.7, ls=":")
    ax.set_yticks([r[0] for r in rows])
    ax.set_yticklabels([r[1] for r in rows], fontsize=FS)
    ax.set_ylim(4.1, -0.6)
    ax.set_xscale("log")
    ticks = [v for v in (0.5, 0.7, 1, 1.4, 2) if min(every) * 0.85 <= v <= max(every) * 1.15]
    ax.set_xlim(min(every) * 0.85, max(every) * 1.15)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{v:g}" for v in ticks])
    ax.minorticks_off()
    ax.set_xlabel("hazard ratio per 1 SD\n(stage-adjusted)", fontsize=FS, labelpad=1)
    _ax_style(ax)
    ax.legend(handles=[Line2D([], [], marker=plotting.HM["LUSC"], color=plotting.INK2, lw=0,
                              markersize=3.6, label="LUSC"),
                       Line2D([], [], marker=plotting.HM["LUAD"], color=plotting.INK2,
                              markerfacecolor="white", lw=0, markersize=3.6, label="LUAD")],
              loc="lower left", bbox_to_anchor=(0.0, 1.0), ncol=2, frameon=False, fontsize=FS,
              handletextpad=0.3, columnspacing=0.8, borderaxespad=0.1)
    _letter(fig, "c")
    save(fig, "suppfig1c_survival")


# Altered gene in panel d: Okabe-Ito colors (safe for color-blind readers), gray for tumors without an alteration.
GENE_COLORS = {"none": "#c9c8c2", "KEAP1": "#CC79A7", "NFE2L2": "#009E73", "CUL3": "#E69F00"}


def _panel_d(rng):
    """Uppsala II, the two fractions in tumors with and without a KEAP1, NFE2L2 or CUL3 alteration.

    Cancer cells over stromal cells; within each, LUSC and LUAD side by side, other against altered, the
    altered tumors colored by gene and every tumor shaped by histology (plotting.HM). Bars: median and
    interquartile range. p over each pair: the two-sided Mann-Whitney test within that histology (s2)."""
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    B = _read_percase("BOMI2")
    B = B[B.histology.isin(["LUSC", "LUAD"]) & B.NRF2path.notna()].copy()
    genes = ["KEAP1", "NFE2L2", "CUL3"]
    assert (B[[f"{g}_alt" for g in genes]].sum(axis=1) <= 1).all(), "a tumor with two of the three alterations"
    B["gene"] = np.select([B[f"{g}_alt"] == 1 for g in genes], genes, "none")
    tests = _read("imaging_driver_tests")
    tests = tests[tests.variable.str.startswith("NRF2path")]
    W, H = 2.1, 3.0
    fig = plt.figure(figsize=(W, H))
    left, width, bottom, gap, top = 0.42, 1.6, 0.36, 0.30, 0.78     # inches; top holds the title and the key
    h = (H - bottom - gap - top) / 2
    axes = {"cancer": fig.add_axes([left / W, (bottom + h + gap) / H, width / W, h / H]),
            "stromal_all": fig.add_axes([left / W, bottom / H, width / W, h / H])}
    xpos = {("LUSC", 0): 0, ("LUSC", 1): 1, ("LUAD", 0): 2.3, ("LUAD", 1): 3.3}
    n_gene = {}
    for m, floor, label in (("cancer", "n_CK_T", "cancer cells"), ("stromal_all", "n_stroma", "stromal cells")):
        ax = axes[m]
        D = B[B[floor] >= 20]
        for (hist, altered), x0 in xpos.items():
            sub = D[(D.histology == hist) & (D.NRF2path == altered)]
            jit = rng.uniform(-0.27, 0.27, len(sub))
            for g in (["none"] if altered == 0 else genes):
                q = (sub.gene == g).to_numpy()
                ax.scatter(x0 + jit[q], 100 * sub[m].to_numpy()[q], s=3.2 if altered == 0 else 6,
                           marker=plotting.HM[hist], color=GENE_COLORS[g], alpha=0.7 if altered == 0 else 0.95,
                           lw=0, zorder=2 if altered == 0 else 3)
                if m == "cancer" and altered == 1:
                    n_gene[g] = n_gene.get(g, 0) + int(q.sum())
            q1, med, q3 = np.quantile(100 * sub[m], [.25, .5, .75])
            ax.plot([x0 - 0.33, x0 + 0.33], [med, med], color=plotting.INK, lw=1.3, zorder=4)
            ax.plot([x0, x0], [q1, q3], color=plotting.INK, lw=0.8, zorder=4)
        for hist in ("LUSC", "LUAD"):
            r = tests[(tests.histology == hist) & (tests.outcome == m)]
            assert len(r) == 1
            p = float(r.p.iloc[0])
            x0, x1 = xpos[(hist, 0)], xpos[(hist, 1)]
            ax.plot([x0, x0, x1, x1], [103, 106, 106, 103], color=plotting.INK2, lw=0.6, clip_on=False)
            ax.text((x0 + x1) / 2, 107.5, "p < 0.001" if p < 0.001 else (f"p = {p:.3f}" if p < 0.1 else f"p = {p:.2f}"),
                    ha="center", va="bottom", fontsize=FS, color=plotting.INK2)
        ax.text(0.0, 1.19, label, transform=ax.transAxes, fontsize=6.0, color=plotting.INK, va="bottom")
        ax.set_ylim(0, 104); ax.set_yticks([0, 50, 100]); ax.set_xlim(-0.6, 3.9)
        ax.set_xticks(list(xpos.values()))
        ax.set_xticklabels(["other", "altered"] * 2 if m == "stromal_all" else [], fontsize=FS)
        ax.set_ylabel("TrxR1+ cells, %", fontsize=FS, labelpad=1)
        _ax_style(ax)
    for x, hist in ((0.5, "LUSC"), (2.8, "LUAD")):
        axes["stromal_all"].text(x, -0.34, hist, ha="center", va="top", fontsize=FS,
                                 transform=axes["stromal_all"].get_xaxis_transform())
    fig.legend(handles=[Line2D([], [], marker="o", ls="", color=GENE_COLORS[g], ms=3.4,
                               label="other" if g == "none" else f"{g} ({n_gene[g]})")
                        for g in ["none"] + genes],
               loc="upper left", bbox_to_anchor=((left - 0.05) / W, (H - top + 0.50) / H), ncol=2, frameon=False,
               fontsize=FS, handletextpad=0.05, columnspacing=0.8, labelspacing=0.25, borderaxespad=0)
    fig.text((left + width / 2) / W, (H - top + 0.55) / H,
             f"KEAP1, NFE2L2 or CUL3 alteration\n({plotting.COHORT_LABEL['BOMI2']})",
             ha="center", va="bottom", fontsize=6.3)
    _letter(fig, "d")
    save(fig, "suppfig1d_genotype")


def _scatter_by_histology(ax, x, y, hist, color, s=4, alpha=0.6):
    """Points by histology (LUSC filled squares, LUAD open circles) with a least-squares line each."""
    x, y, hist = np.asarray(x, float), np.asarray(y, float), np.asarray(hist)
    for h in ("LUSC", "LUAD"):
        m = hist == h
        if h == "LUSC":
            ax.scatter(x[m], y[m], s=s, marker=plotting.HM[h], color=color, alpha=alpha, lw=0)
        else:
            ax.scatter(x[m], y[m], s=s, marker=plotting.HM[h], facecolors="none", edgecolors=color,
                       alpha=alpha, lw=0.5)
        ok = m & np.isfinite(x) & np.isfinite(y)
        b1, b0 = np.polyfit(x[ok], y[ok], 1)
        xs = np.array([x[ok].min(), x[ok].max()])
        ax.plot(xs, b0 + b1 * xs, color=plotting.INK if color == plotting.INK2 else color, lw=1.0,
                zorder=4, ls="-" if h == "LUSC" else (0, (3, 1.5)))


# Geometry of the two-row panels e and f, inches: axes width and height, and the bottom of the
# upper and the lower row (room below for the x label and the histology key, above each row for
# its two-line subtitle).
GRID_W, GRID_H = 1.10, 0.95
GRID_BOTTOM = (1.97, 0.62)


def _subtitle(ax, name, text, height, color=plotting.INK):
    """A two-line subtitle above an axis: the descriptive text, and above it the group name."""
    ax.text(0.5, 1.04, text, transform=ax.transAxes, ha="center", va="bottom", fontsize=FS,
            color=plotting.INK2)
    ax.text(0.5, 1.04 + 0.14 / height, name, transform=ax.transAxes, ha="center", va="bottom",
            fontsize=6.2, color=color)


def _histology_legend(fig, anchor):
    """The key for the histology markers (LUSC square, LUAD circle), placed at anchor in figure units."""
    from matplotlib.lines import Line2D
    fig.legend(handles=[Line2D([], [], marker=plotting.HM["LUSC"], color=plotting.INK2, lw=0.9,
                               markersize=3.6, label="LUSC"),
                        Line2D([], [], marker=plotting.HM["LUAD"], color=plotting.INK2,
                               markerfacecolor="none", lw=0.9, ls=(0, (3, 1.5)), markersize=3.6,
                               label="LUAD")],
               loc="lower center", bbox_to_anchor=anchor, ncol=2, frameon=False, fontsize=FS,
               handletextpad=0.3, columnspacing=0.8)


def _panel_e():
    """Cancer-cell against stromal fraction, as measured and adjusted for slide and background staining."""
    import matplotlib.pyplot as plt
    points = _read("compartment_points")
    corr = _read("compartment_correlation")
    cc = corr[corr.pair == "cancer cells with stromal cells"].set_index(["cohort", "histology"])
    for (cohort, h), g in points.groupby(["cohort", "histology"]):
        r = float(np.corrcoef(g.cancer_adjusted, g.stroma_adjusted)[0, 1])
        if abs(r - cc.loc[(cohort, h), "rho_slide_background"]) > 1e-6:
            raise ValueError(f"panel e, {cohort} {h}: the drawn residuals give {r}, the table "
                             f"{cc.loc[(cohort, h), 'rho_slide_background']}")
    lim = 1.15 * float(np.abs(points[["cancer_adjusted", "stroma_adjusted"]].to_numpy()).max())
    fig = plt.figure(figsize=(3.3, 3.4))
    fw, fh = fig.get_size_inches()
    for j, (kind, head) in enumerate((("meas", "as measured"), ("adj", "adjusted for slide and background staining"))):
        left = 0.45 + j * (GRID_W + 0.55)
        fig.text((left + GRID_W / 2) / fw, 3.33 / fh, head, ha="center", va="top", fontsize=6.2)
        for i, cohort in enumerate(paths.COHORTS):
            ax = fig.add_axes([left / fw, GRID_BOTTOM[i] / fh, GRID_W / fw, GRID_H / fh])
            g = points[points.cohort == cohort]
            if kind == "meas":
                _scatter_by_histology(ax, g.f_cancer * 100, g.f_stroma * 100, g.histology, plotting.INK2)
                ax.set_xlim(0, 100)
                ax.set_ylim(0, 100)
                ax.set_xticks([0, 50, 100])
                ax.set_yticks([0, 50, 100])
                ax.set_ylabel("stromal cells, %", fontsize=FS, labelpad=1)
                column = "rho_raw"
            else:
                _scatter_by_histology(ax, g.cancer_adjusted, g.stroma_adjusted, g.histology, plotting.INK2)
                ax.set_xlim(-lim, lim)
                ax.set_ylim(-lim, lim)
                ax.axhline(0, color=AXIS, lw=0.4, ls=":", zorder=0)
                ax.axvline(0, color=AXIS, lw=0.4, ls=":", zorder=0)
                ax.set_ylabel("stromal, adjusted rank", fontsize=FS, labelpad=1)
                column = "rho_slide_background"
            _ax_style(ax)
            _subtitle(ax, plotting.COHORT_LABEL[cohort],
                      f"rho: LUSC {cc.loc[(cohort, 'LUSC'), column]:.2f}, "
                      f"LUAD {cc.loc[(cohort, 'LUAD'), column]:.2f}", GRID_H)
            if i == 0:
                ax.tick_params(labelbottom=False)
            else:
                ax.set_xlabel("cancer cells, %" if kind == "meas" else "cancer, adjusted rank",
                              fontsize=FS, labelpad=1)
    _histology_legend(fig, (0.5, 0.0))
    _letter(fig, "e")
    save(fig, "suppfig1e_compartments")


def _panel_f():
    """Bulk tumor TXNRD1 against the cancer-cell and the stromal fraction, Uppsala II."""
    import matplotlib.pyplot as plt
    if not out("bulk_rna.csv").exists():
        print("  no bulk_rna.csv; panel f is skipped")
        return
    T4 = _read("bulk_rna")
    X = _read("bulk_per_patient")
    W = 1.25
    fig = plt.figure(figsize=(1.9, 3.4))
    fw, fh = fig.get_size_inches()
    axes = {}
    for i, (m, color, title) in enumerate((("cancer", TUM, "cancer cells"),
                                           ("stromal_all", STROMA, "stromal cells"))):
        ax = axes[m] = fig.add_axes([0.5 / fw, GRID_BOTTOM[i] / fh, W / fw, GRID_H / fh])
        g = X.dropna(subset=[m])
        _scatter_by_histology(ax, g[m] * 100, g.bulk_TXNRD1, g.histology, color)
        text = []
        for h in ("LUSC", "LUAD"):
            gg = g[g.histology == h]
            r, _ = spearmanr(gg[m], gg.bulk_TXNRD1)
            ref = T4[(T4.histology == h) & (T4.test == f"Spearman bulk TXNRD1 vs {m}")].iloc[0]
            if abs(r - ref.estimate) > 1e-6 or len(gg) != int(float(ref.n)):
                raise ValueError(f"panel f, {m} {h}: the drawn points give rho {r}, the table {ref.estimate}")
            text.append(f"{h} {r:.2f}")
        _subtitle(ax, title, "rho: " + ", ".join(text), GRID_H, color=color)
        ax.set_xlim(0, 100)
        ax.set_xticks([0, 50, 100])
        _ax_style(ax)
        ax.set_ylabel("bulk TXNRD1 (log2)", fontsize=FS, labelpad=1)
        if i == 0:
            ax.tick_params(labelbottom=False)
        else:
            ax.set_xlabel("TrxR1-positive cells, %", fontsize=FS, labelpad=1)
    lo = min(a.get_ylim()[0] for a in axes.values())
    hi = max(a.get_ylim()[1] for a in axes.values())
    for a in axes.values():
        a.set_ylim(lo, hi)
    _histology_legend(fig, ((0.5 + W / 2) / fw, 0.0))
    _letter(fig, "f")
    save(fig, "suppfig1f_bulk")


# ---------------------------------------------------------------- s7: TrxR1 by cell type

# The groups of the profile: cancer cells, then the stromal cell types in the order of Figure 1 (T cells,
# fibroblasts, macrophages) and the stromal cells without a lineage marker. Uppsala I has no CD68 and no
# alpha-SMA, so its cells without a lineage marker hold the fibroblasts and macrophages that Uppsala II separates
# (redox.cells.LINEAGES); the panel shows Uppsala II only. All stromal cells together are computed only to check
# the stromal fraction against percase (they are Figure 1b's measure, not written here). The summary keeps the
# rows the manuscript reads: every group over both histologies for Figure 1d (Uppsala II). The rows of each
# stromal cell type per histology are not written; lineage_profile_patients.csv holds every group of
# both cohorts for s8.
PROFILE_GROUPS = {"BOMI2": ["cancer", "T", "Fib", "Mac", "Other"],
                  "BOMI1": ["cancer", "T", "Other"]}
MIN_PROFILE = 20          # cells of the group per patient, as in panel b


def _profile_histologies(cohort, group):
    """The summary rows of one group: both histologies together in Uppsala II (Figure 1d)."""
    return ["both"] if cohort == "BOMI2" else []


def s7_lineage_profile():
    """TrxR1-positive fraction of cancer cells and of each stromal cell type per patient, summary, panel."""
    long, summary = [], []
    for cohort in paths.COHORTS:
        cells = cellsmod.load_cells(cohort, columns=["Case", "tissue", "lineage", "TrxR1"])
        cells["Case"] = cells.Case.astype(str)
        stroma = cells[cells.tissue == "Stroma"]
        per = _read_percase(cohort).set_index("Case")
        wide = pd.DataFrame(index=pd.Index(sorted(cells.Case.unique()), name="Case"))
        for g in dict.fromkeys(PROFILE_GROUPS[cohort] + ["stromal_all"]):
            sub = (cells[(cells.tissue == "Tumor") & (cells.lineage == "CK")] if g == "cancer"
                   else stroma if g == "stromal_all" else stroma[stroma.lineage == g])
            wide[g], wide[f"n_{g}"] = _fraction(sub, MIN_PROFILE)
        # the cancer-cell and stromal fractions are those of percase (same cells, same floor)
        for g in ("cancer", "stromal_all"):
            check = per[g].reindex(wide.index)
            if not np.allclose(check.fillna(-1), wide[g].fillna(-1), atol=1e-12):
                raise SystemExit(f"s7 {cohort}: the {g} fraction differs from percase_{cohort}.csv")
        wide["histology"] = per.histology.reindex(wide.index)
        wide = wide[wide.histology.isin(HISTS)]
        for g in PROFILE_GROUPS[cohort]:
            v = wide[[g, f"n_{g}", "histology"]].dropna(subset=[g])
            long.append(pd.DataFrame(dict(cohort=cohort, Case=v.index, histology=v.histology.to_numpy(),
                                          group=g, fraction=v[g].to_numpy(), n_cells=v[f"n_{g}"].to_numpy())))
            for h in _profile_histologies(cohort, g):
                x = v.loc[v.histology.isin(HISTS if h == "both" else [h]), g].to_numpy()
                q1, med, q3 = np.quantile(x, [.25, .5, .75])
                summary.append(dict(cohort=cohort, histology=h, group=g, n=len(x),
                                    median=med, q1=q1, q3=q3))
        print(f"s7 {cohort}: {len(wide)} LUSC and LUAD patients")
    groups = ("cancer (cytokeratin-positive cells of the tumor compartment), T (CD3), Fib (alpha-SMA, BOMI2), "
              "Mac (CD68, BOMI2), Other (stromal cells without a lineage marker; in BOMI1 this includes "
              "fibroblasts and macrophages)")
    _write(pd.concat(long, ignore_index=True), "lineage_profile_patients",
           "Per patient, the TrxR1-positive fraction of cancer cells and of each stromal cell type (LUSC and "
           "LUAD, at least 20 cells of the group); section s8 fits its survival block by cell type on it. Written by "
           "tumor_vs_stroma.py section s7 from cells.tsv and percase_<cohort>.csv. Patient-level: it stays with "
           "the data and is never shared.",
           [("cohort", "BOMI2 (Uppsala II) or BOMI1 (Uppsala I)"), ("Case", "patient identifier"),
            ("histology", "LUSC or LUAD"), ("group", groups),
            ("fraction", "TrxR1-positive fraction of the patient's cells of the group (mean call)"),
            ("n_cells", "cells of the group of the patient")])
    _write(pd.DataFrame(summary), "lineage_profile_summary",
           "Per group: patients and the median and quartiles of the TrxR1-positive fraction, BOMI2, LUSC and "
           "LUAD together, the rows Figure 1d draws (cancer cells and the four stromal cell types). Written by "
           "tumor_vs_stroma.py section s7.",
           [("cohort", "BOMI2"), ("histology", "both (LUSC and LUAD together)"),
            ("group", "as in lineage_profile_patients"),
            ("n", "patients with at least 20 cells of the group"), ("median", "median fraction"),
            ("q1", "25th percentile"), ("q3", "75th percentile")])
    _panel_profile()


def _panel_profile():
    """Figure 1d: Uppsala II, TrxR1-positive fraction of cancer cells, T cells, fibroblasts, macrophages and stromal cells
    without a lineage marker ("no marker"), median per patient (bar, value at right) with the interquartile range; five
    compact rows in the height of the three rows of Figure 1c and 1e; lineage colors of redox.plotting, cancer cells as
    in Supplementary Figure 1. Letter, title, subtitle, labels and font sizes as in the article panel, which the separate
    layout code composes from lineage_profile_summary.csv as well."""
    import matplotlib.pyplot as plt
    plotting.style()
    s = _read("lineage_profile_summary")
    s = s[(s.cohort == "BOMI2") & (s.histology == "both")].set_index("group")
    rows = [("cancer", "cancer cells", 0.30, TUM), ("T", "T cells", -0.42, plotting.LC["T"]),
            ("Fib", "fibroblasts", -1.05, plotting.LC["Fib"]), ("Mac", "macrophages", -1.68, plotting.LC["Mac"]),
            ("Other", "no marker", -2.31, plotting.LC["Other"])]
    fig = plt.figure(figsize=(1.75, 2.05))
    ax = fig.add_axes([0.42, 0.20, 0.42, 0.60])
    for g, lab, y, col in rows:
        r = s.loc[g]
        ax.barh(y, r["median"], height=0.40, color=col, lw=0, zorder=2)
        ax.plot([r.q1, r.q3], [y, y], color=plotting.INK, lw=0.8, zorder=3)
        ax.text(-0.04, y, lab, ha="right", va="center", fontsize=6.3)
        ax.text(r.q3 + 0.03, y, f"{100 * r['median']:.0f} %", ha="left", va="center", fontsize=FS,
                color=plotting.INK2)
    ax.set_xlim(0, 1); ax.set_ylim(-2.5, 0.55)
    ax.set_yticks([]); ax.set_xticks([0, 0.5, 1]); ax.set_xticklabels(["0", "50", "100"])
    ax.set_xlabel("TrxR1+ cells, %", fontsize=6.3, color=plotting.INK2, labelpad=2)
    _ax_style(ax)
    ax.spines["left"].set_visible(False)
    fig.text(0.01, 0.94, "d", fontsize=9.5, fontweight="bold")
    fig.text(0.11, 0.94, "TrxR1+ cells", fontsize=7.5, fontweight="bold")
    fig.text(0.11, 0.875, f"{plotting.COHORT_LABEL['BOMI2']}, per patient", fontsize=6.0, color=plotting.INK2)
    save(fig, "fig1d_lineage_profile")


# ---------------------------------------------------------------- s8: TrxR1 in all cells and per stromal cell type, survival

# The measures of the survival block of Supplementary Table 1a: the TrxR1-positive fraction
# of all cells of the patient (cancer and stromal, every lineage; patients who carry both the cancer-cell and the
# stromal fraction, as in the joint models), the stromal fraction for reference, and the fraction of each stromal
# cell type (at least 20 cells of the type, lineage_profile_patients.csv). Cox models adjusted for stage, per 1 SD,
# pooled by fixed effect where both cohorts have the marker.
SURVIVAL_MEASURES = [("all_cells", "all cells (cancer and stromal)"), ("stromal_all", "stromal cells (Figure 1b)"),
                     ("T", "T cells"), ("Fib", "fibroblasts"), ("Mac", "macrophages"),
                     ("Other", "stromal cells without a lineage marker")]


def _survival_frame():
    """Both cohorts, LUSC and LUAD: the clinical columns, the stromal and cancer-cell fractions, and the
    TrxR1-positive fraction of all cells and of each stromal cell type (the frame of s8 and s9)."""
    D = _patients()
    D["Case"] = D.Case.astype(str)
    P = _read("lineage_profile_patients")
    P["Case"] = P.Case.astype(str)
    frames = []
    for cohort in paths.COHORTS:
        cells = cellsmod.load_cells(cohort, columns=["Case", "tissue", "lineage", "TrxR1"])
        cells["Case"] = cells.Case.astype(str)
        f, _ = _fraction(cells, 1)
        w = P[P.cohort == cohort].pivot(index="Case", columns="group", values="fraction")
        w = w.drop(columns=["cancer"])
        w["all_cells"] = f.reindex(w.index.union(f.index))
        frames.append(w.reset_index().assign(cohort=cohort))
    D = D.merge(pd.concat(frames, ignore_index=True), on=["cohort", "Case"], how="left")
    D.loc[D.cancer.isna() | D.stromal_all.isna(), "all_cells"] = np.nan
    return D



def _figure1b_lusc_pooled():
    """The pooled LUSC hazard ratio of Figure 1b, from chapter1/fig1b_models.csv (the s8 and s9 checks)."""
    f = paths.output(1, "fig1b_models.csv")
    if not f.exists():
        raise SystemExit(f"{f} is missing: run chapter1.py first")
    m = pd.read_csv(f)
    return float(m[(m.histology == "LUSC") & (m.cohort == "pooled")].HR.iloc[0])

def s8_lineage_survival():
    """Overall survival against the TrxR1-positive fraction of all cells and of each stromal cell type."""
    D = _survival_frame()
    rows = []
    for m, label in SURVIVAL_MEASURES:
        for h in HISTS[::-1]:                                   # LUSC first, as the table prints it
            per = {}
            for cohort in paths.COHORTS:
                x = D[(D.cohort == cohort) & (D.histology == h)]
                if x[m].notna().sum() == 0:
                    continue
                per[cohort] = survival.hr_per_sd(x, m)
                rows.append(dict(measure=m, label=label, histology=h, cohort=cohort, **per[cohort]))
            if len(per) == 2:
                rows.append(dict(measure=m, label=label, histology=h, cohort="pooled",
                                 **survival.meta_fixed([per["BOMI2"], per["BOMI1"]])))
    T8 = pd.DataFrame(rows)
    ref = T8[(T8.measure == "stromal_all") & (T8.histology == "LUSC") & (T8.cohort == "pooled")].iloc[0]
    if abs(ref.HR - _figure1b_lusc_pooled()) > 1e-9:
        raise SystemExit(f"s8: the stromal reference row does not reproduce Figure 1b ({ref.HR:.4f})")
    _write(T8, "lineage_survival",
           "Overall survival against the TrxR1-positive fraction of all cells and of each stromal cell type, "
           "Cox model adjusted for stage, per cohort and histology and pooled by fixed effect where both cohorts "
           "carry the measure. Written by tumor_vs_stroma.py section s8.",
           [("measure", "all_cells, stromal_all (the stromal fraction of percase_<cohort>.csv, the reference), "
                        "T, Fib, Mac or Other (as in lineage_profile_patients)"),
            ("label", "the row label of Supplementary Table 1a"), ("histology", "LUSC or LUAD"),
            ("cohort", "BOMI2, BOMI1 or pooled")] + MODEL_COLUMNS)
    # The printed block: every measure but the stromal reference (printed in the block above it), then the
    # patients (deaths) row, the range over the rows of the block for each histology and cohort.
    printed = T8[T8.measure != "stromal_all"]
    shown = []
    for _, r in printed.iterrows():
        shown.append(dict(measure=r.measure, label=r.label, histology=r.histology, cohort=r.cohort,
                          text=f"{r.HR:.2f} ({r.lo:.2f} to {r.hi:.2f}), p {_p_text(r.p)}",
                          patients=int(r.n), deaths=int(r.events)))
    for h in HISTS[::-1]:
        for cohort in ("BOMI2", "BOMI1", "pooled"):
            g = printed[(printed.histology == h) & (printed.cohort == cohort)]
            lo, hi, dlo, dhi = int(g.n.min()), int(g.n.max()), int(g.events.min()), int(g.events.max())
            shown.append(dict(measure="range", label="patients (deaths)", histology=h, cohort=cohort,
                              text=f"{lo} to {hi} ({dlo} to {dhi})", patients_min=lo, patients_max=hi,
                              deaths_min=dlo, deaths_max=dhi))
    counts = ["patients", "deaths", "patients_min", "patients_max", "deaths_min", "deaths_max"]
    _write(pd.DataFrame(shown).astype({c: "Int64" for c in counts}), "supptable1a_celltypes",
           "The printed cells of the survival block by cell type of Supplementary Table 1a: one row per "
           "measure, histology and cohort, then the patients (deaths) row of the block. Written by "
           "tumor_vs_stroma.py section s8 from lineage_survival.csv (its stromal_all rows, the stromal "
           "reference printed in the block above, are left out).",
           [("measure", "as in lineage_survival: all_cells, T, Fib, Mac or Other; range for the patients "
                        "(deaths) row"),
            ("label", "row label"), ("histology", "LUSC or LUAD"),
            ("cohort", "BOMI2, BOMI1 or pooled"),
            ("text", "the printed cell: HR (95 % CI), p; for the range row, patients (deaths) as the "
                     "smallest to the largest over the rows of the block"),
            ("patients", "patients in the model (empty on the range row)"),
            ("deaths", "deaths among them (empty on the range row)"),
            ("patients_min", "range row: fewest patients over the rows of the block, same histology and cohort"),
            ("patients_max", "range row: most patients"),
            ("deaths_min", "range row: fewest deaths"),
            ("deaths_max", "range row: most deaths")])


# ---------------------------------------------------------------- s9: further covariates (Supplementary Table 1a)

# The sensitivity models of the covariate block of Supplementary Table 1a: the model of Figure 1b with age (years), sex and
# WHO performance status (one linear term, 0 to 4; redox.clinical.load_covariates) and the tissue microarray slide
# (the slide holding most of the patient's cells, the slide column of percase_<cohort>.csv; one indicator per slide
# among the analyzed patients, the first slide in sorted order as reference) added to stage together. No patient of
# these analysis sets lacks a covariate (checked below); one who did would leave that model. The stage model is the
# reference: it reproduces Figure 1b (checked below). The block fits the stage model and the fully
# adjusted model only.
COVARIATE_MODELS = [("stage", "stage", []),
                    ("all", "stage, age, sex, WHO performance status and slide", ["age", "male", "ps", "slide"])]
# The measure: the stromal fraction of Figure 1b. The same models for the niche connectivity of Figure 5g are
# Supplementary Table 3d, fitted by chapter4.py s6 with the same covariates and the same code path.
FULL_MEASURES = [("stromal_all", "stromal cells (Figure 1b)")]


def _covariate_fit(x, measure, extra):
    """hr_per_sd with stage and the extra covariates; "slide" stands for one indicator per slide."""
    plain = [c for c in extra if c != "slide"]
    need = [measure, "OS_weeks", "event", "stage"] + plain + (["slide"] if "slide" in extra else [])
    x = x.dropna(subset=need)
    x = x[x.OS_weeks > 0]
    covars = ["stage"] + plain
    if "slide" in extra:
        S = pd.get_dummies(x.slide.astype(str), prefix="slide", drop_first=True).astype(float)
        x = pd.concat([x, S], axis=1)
        covars += list(S.columns)
    return survival.hr_per_sd(x, measure, covars=covars)


def s9_covariates():
    """Supplementary Table 1a, covariate block: the Figure 1b model with age, sex, performance status and slide
    added to stage."""
    D = _survival_frame()
    cov = pd.concat([clinicalmod.load_covariates(c).assign(cohort=c) for c in paths.COHORTS], ignore_index=True)
    sl = pd.concat([_read_percase(c)[["Case", "slide"]].assign(cohort=c) for c in paths.COHORTS], ignore_index=True)
    D = D.merge(cov, on=["cohort", "Case"], how="left").merge(sl, on=["cohort", "Case"], how="left")
    plan = [(m, lab, COVARIATE_MODELS) for m, lab in FULL_MEASURES]
    rows = []
    for m, label, models in plan:
        for key, model_label, extra in models:
            for h in HISTS[::-1]:                               # LUSC first, as the table prints it
                per = {}
                for cohort in paths.COHORTS:
                    x = D[(D.cohort == cohort) & (D.histology == h)]
                    if x[m].notna().sum() == 0:
                        continue
                    per[cohort] = _covariate_fit(x, m, extra)
                    rows.append(dict(measure=m, label=label, model=key, model_label=model_label,
                                     histology=h, cohort=cohort, **per[cohort]))
                if len(per) == 2:
                    rows.append(dict(measure=m, label=label, model=key, model_label=model_label, histology=h,
                                     cohort="pooled", **survival.meta_fixed([per["BOMI2"], per["BOMI1"]])))
    T9 = pd.DataFrame(rows)
    # Checks: every model of a measure uses the same patients; the stage model reproduces Figure 1b.
    n = T9.groupby(["measure", "histology", "cohort"]).n.nunique()
    if (n > 1).any():
        raise SystemExit(f"s9: a covariate is missing for some patients: {n[n > 1].index.tolist()}")
    st = T9[T9.model == "stage"].set_index(["measure", "histology", "cohort"])
    if abs(st.loc[("stromal_all", "LUSC", "pooled"), "HR"] - _figure1b_lusc_pooled()) > 1e-9:
        raise SystemExit("s9: the stromal stage model does not reproduce Figure 1b")
    _write(T9, "covariate_survival",
           "Overall survival against the stromal TrxR1 fraction (Figure 1b), adjusted for stage (the reference, "
           "equal to Figure 1b) and for stage, age, sex, WHO performance status and TMA slide together. Cox model "
           "per cohort and histology, per 1 SD, pooled by fixed effect. Written by tumor_vs_stroma.py section s9. "
           "(The same models of niche connectivity, Figure 5g, are in chapter4/ch4_niche_survival.csv.)",
           [("measure", "stromal_all (Figure 1b)"),
            ("label", "the row label of the measure in Supplementary Table 1a"),
            ("model", "stage or all"),
            ("model_label", "the covariates of the model, as the table prints them"),
            ("histology", "LUSC or LUAD"), ("cohort", "BOMI2, BOMI1 or pooled")]
           + [(c, m.replace("adjusted for stage", "adjusted for the covariates of the model")) for c, m in MODEL_COLUMNS])
    # The printed block: the stromal fraction fully adjusted (its stage model is printed above). The connectivity
    # rows and their patients (deaths) row are Supplementary Table 3d (chapter4.py s6).
    shown = []

    def cell(r):
        """A hazard ratio as the table prints it: HR (lower to upper), p."""
        return f"{r.HR:.2f} ({r.lo:.2f} to {r.hi:.2f}), p {_p_text(r.p)}"
    for _, r in T9.iterrows():
        if r.measure == "stromal_all" and r.model != "stage":
            shown.append(dict(block="stromal", measure=r.measure, model=r.model, label=r.model_label,
                              histology=r.histology, cohort=r.cohort, text=cell(r)))
    _write(pd.DataFrame(shown), "supptable1a_covariates",
           "The printed cells of the covariate block of Supplementary Table 1a. Written by tumor_vs_stroma.py "
           "section s9 from covariate_survival.csv.",
           [("block", "stromal (Figure 1b, fully adjusted)"),
            ("measure", "as in covariate_survival.csv"), ("model", "as in covariate_survival.csv"),
            ("label", "the printed row label"), ("histology", "LUSC or LUAD"),
            ("cohort", "BOMI2, BOMI1 or pooled"),
            ("text", "the printed cell: HR (95 % CI), p")])


# ---------------------------------------------------------------- run

SECTIONS = {"s1": s1_per_patient,
            "s2": s2_driver_tests,
            "s3": s3_cancer_and_stroma,
            "s4": s4_compartments,
            "s5": s5_supptable1a,
            "s6": s6_panels,
            "s7": s7_lineage_profile,
            "s8": s8_lineage_survival,
            "s9": s9_covariates}


def main():
    """Command line: run the sections named by --only, or all of them in order."""
    ap = argparse.ArgumentParser(description=__doc__.split("RUN:")[0].strip(),
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", nargs="*", choices=sorted(SECTIONS),
                    help="run these sections only; the default runs all of them in order")
    a = ap.parse_args()
    for name in (a.only or sorted(SECTIONS)):
        print(f"[{name}] {SECTIONS[name].__doc__.strip().splitlines()[0]}")
        SECTIONS[name]()


if __name__ == "__main__":
    main()
