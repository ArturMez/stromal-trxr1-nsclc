"""
The TrxR1 staining behind the positivity call: intensity, tissue microarray slide and background.

WHAT: the cytoplasmic TrxR1 intensity of every cell of cells.tsv, the tissue microarray (TMA) slide
      of every core and patient, the background index of a group of cells, and the rank correlation
      of two per-patient measures with slide and background removed.
FROM: raw_root/Raw_original_<cohort>/<cell export> (the TrxR1 channel only, together with the four
      columns that identify a cell) and work_root/intermediate/<cohort>/cells.tsv.
TO:   DataFrames and floats; nothing is written here.
WHY:  cells.tsv carries the TrxR1 call, not the intensity it was made from. The tumor-versus-stroma
      analyses (tumor_vs_stroma.py), chapter1.py (the slide and background checks) and chapter4.py
      need the intensity itself: where the cut sits among cancer cells, the mean intensity of the
      TrxR1-negative cells of a patient as an index of background staining, and the slide a core
      was stained on. They read it through this module, so that the alignment of the raw export
      with cells.tsv is checked in one place.

The raw export is read with the missing-value spellings mif/02_cells.R uses, restricted to the
Tumor and Stroma categories and to the cores cells.tsv keeps. mif/02_cells.R writes cells.tsv from
the same rows in the same order, so the two tables are aligned row by row; the alignment is
checked on the core, the tissue category and both raw positions of every cell, and the call in
cells.tsv is checked against the intensity and the cut recorded in thresholds.csv.
"""
import numpy as np
import pandas as pd
from scipy.stats import rankdata, t as t_dist

from . import cells as cellsmod
from . import paths

# The raw exports mif/02_cells.R reads, per cohort: folder under raw_root and file name (the Uppsala I
# file name comes from config.json "input_files", key bomi1_cell_export).
RAW_EXPORT = {"BOMI2": ("Raw_original_BOMI2", "BOMI2_REDOX_merged_cell_seg_data.txt"),
              "BOMI1": ("Raw_original_BOMI1", paths.input_file("bomi1_cell_export"))}

# The TrxR1 channel (Opal 570, cytoplasm, mean), the channel mif/02_cells.R calls TrxR1.
TRXR1_CHANNEL = "Cytoplasm Opal 570 Mean (Normalized Counts, Total Weighting)"

# The slide number inside a core name. A core is named after its TMA slide, so the slide is the
# number that follows the cohort prefix of the name.
SLIDE_PATTERN = {"BOMI2": r"BOMI2_Redox_(\d+)_Core", "BOMI1": r"RedOx_BOMI_1_(\d+)_Core"}

# Least cells for a background index: the TrxR1-negative cells of the group must number this many.
MIN_BACKGROUND_CELLS = 10


def cells_with_intensity(cohort, columns=None):
    """cells.tsv of a cohort with one more column, TrxR1_intensity, read from the raw export.

    columns  the cells.tsv columns to return besides TrxR1_intensity; none returns every column of
             cells.tsv, positions as x_raw and y_raw (the alignment is checked on the core, the
             tissue category and the raw positions, so those are always read).
    """
    every = (["core", "Case", "tissue", "x_raw", "y_raw"]
             + [marker for marker, _ in cellsmod.MARKERS[cohort]] + ["lineage", "TrxR1"])
    if columns is None:
        columns = every
    want = list(dict.fromkeys(["core", "tissue", "x_raw", "y_raw", "TrxR1"] + list(columns)))
    cells = cellsmod.load_cells(cohort, columns=want)
    folder, name = RAW_EXPORT[cohort]
    source = paths.raw(folder, name)
    raw = pd.read_csv(source, sep="\t",
                      usecols=["Sample Name", "Tissue Category", "Cell X Position",
                               "Cell Y Position", TRXR1_CHANNEL],
                      na_values=["NA", "#N/A", ""], keep_default_na=False,
                      float_precision="round_trip", low_memory=False)
    raw = raw[raw["Tissue Category"].isin(["Tumor", "Stroma"])
              & raw["Sample Name"].isin(set(cells.core))].reset_index(drop=True)
    if len(raw) != len(cells):
        raise ValueError(f"{source} holds {len(raw)} Tumor and Stroma cells of the kept cores, "
                         f"cells.tsv {len(cells)}; rerun mif/02_cells.R")
    same = ((raw["Sample Name"].to_numpy() == cells.core.to_numpy()).all()
            and (raw["Tissue Category"].to_numpy() == cells.tissue.to_numpy()).all()
            and np.array_equal(raw["Cell X Position"].to_numpy(dtype=float),
                               cells.x_raw.to_numpy(dtype=float))
            and np.array_equal(raw["Cell Y Position"].to_numpy(dtype=float),
                               cells.y_raw.to_numpy(dtype=float)))
    if not same:
        raise ValueError(f"{source} and cells.tsv are not aligned cell by cell; rerun mif/02_cells.R")
    cells["TrxR1_intensity"] = raw[TRXR1_CHANNEL].to_numpy(dtype=float)
    cut = trxr1_cut(cohort)
    known = cells.TrxR1_intensity.notna()
    call = (cells.TrxR1_intensity[known] > cut).astype(int).to_numpy()
    if not (cells.TrxR1[known].astype(int).to_numpy() == call).all() or cells.TrxR1[~known].notna().any():
        raise ValueError(f"the TrxR1 call of cells.tsv does not follow the intensity and the cut "
                         f"{cut} of thresholds.csv ({cohort})")
    return cells[list(columns) + ["TrxR1_intensity"]]


def trxr1_cut(cohort):
    """The TrxR1 positivity cut of a cohort, as mif/02_cells.R recorded it in thresholds.csv."""
    t = pd.read_csv(paths.intermediate(cohort, "thresholds.csv"), float_precision="round_trip")
    return float(t.loc[t.marker == "TrxR1", "cut"].iloc[0])


def slide_of(core, cohort):
    """The TMA slide of each core name, as a string (the digits of the name, leading zeros kept)."""
    return pd.Series(core).str.extract(SLIDE_PATTERN[cohort], expand=False)


def patient_slide(cells, cohort):
    """Each patient's slide: the slide of most of the patient's cells, the lowest name on a tie.

    Most patients have all their cores on one slide. For the others the slide that holds most of
    their cells is taken (pandas mode, whose first value is the smallest of the tied names).
    """
    s = cells.assign(slide=slide_of(cells.core, cohort).to_numpy()).dropna(subset=["slide"])
    return s.groupby("Case").slide.agg(lambda v: v.mode().iloc[0])


def fraction_and_background(group, min_cells=20, min_background=MIN_BACKGROUND_CELLS):
    """Per patient: TrxR1-positive fraction and background index of a group of cells.

    The call is read with a missing call as negative, so that the fraction and the background use
    the same cells: the fraction is the mean call over every cell of the group (patients with fewer
    than min_cells cells left empty), and the background index is the mean intensity of the
    TrxR1-negative cells (patients with fewer than min_background such cells left empty).
    """
    pos = group.TrxR1.astype(float).fillna(0.0)
    g = pos.groupby(group.Case)
    frac = g.mean().where(g.size() >= min_cells)
    neg = group[pos.to_numpy() == 0].groupby("Case").TrxR1_intensity
    background = neg.mean().where(neg.size() >= min_background)
    return frac, background


def partial_rank_correlation(a, b, covariates, return_residuals=False):
    """Correlation of the ranks of a and b after the covariates are regressed out of both.

    Ranks with ties averaged; each rank vector is regressed by least squares on an intercept and
    the covariates (slide indicators, ranks of background indices), and the Pearson correlation of
    the two residual vectors is returned with a t-test p on n - k - 2 degrees of freedom, k being
    the number of covariates. With no covariates this is Spearman's rho.
    """
    ra, rb = rankdata(a), rankdata(b)
    X = np.column_stack([np.ones(len(ra))] + list(covariates))
    ea = ra - X @ np.linalg.lstsq(X, ra, rcond=None)[0]
    eb = rb - X @ np.linalg.lstsq(X, rb, rcond=None)[0]
    r = float(np.corrcoef(ea, eb)[0, 1])
    df = len(ra) - X.shape[1] - 1
    p = float(2 * t_dist.sf(abs(r) * np.sqrt(df / (1 - r * r)), df))
    return (r, p, ea, eb) if return_residuals else (r, p)


def slide_indicators(slides):
    """Slide indicator columns for partial_rank_correlation, the first slide as reference."""
    S = pd.get_dummies(pd.Series(slides).reset_index(drop=True), drop_first=True).astype(float)
    return [S[c].to_numpy() for c in S.columns]
