#!/usr/bin/env python3
"""
Chapter 4 analysis: TrxR1 niches, their composition, their connectivity and overall survival.

WHAT: six steps that make Figure 5a-h, Supplementary Figure 6 and Supplementary Tables 3a to 3d, plus
      ch4_niche_survival.csv (every Cox model behind Supplementary Table 3d): per-core niche measures,
      their patient means, every statistic the chapter prints, the individual panels (each drawn from the
      table it belongs to), the niche-area densities behind the density block of Supplementary Table 3b,
      and the niche measures against overall survival. The chapter and file names (chapter4, fig4*,
      supptable4*) are historical; the items they hold are the ones named here.
      Naming: "confluence" in the code (columns biggest_30, ch4_confluence_*) is the niche connectivity of
      the paper: the largest group of linked TrxR1-positive cells (links at <= 30 um) over all
      TrxR1-positive stromal cells of the core.
FROM: work_root/intermediate/<cohort>/cells.tsv (positions, lineage, TrxR1 call), percase.csv (the
      stromal TrxR1 fraction, the T-cell proportion and the stromal cell density of Chapter 1),
      clinical.csv, clinical_full.csv (s6: age, sex and WHO performance status, through
      redox.clinical.load_covariates), and the example-core file named by the "examples" entry of
      config.json.
TO:   work_root/output/chapter4/, every table with a sibling .md describing its columns:
        s1  per core, both cohorts, cores passing RULE_CH4 (>= 60 stromal cells, >= 25 of them
            TrxR1-positive and >= 25 TrxR1-negative)
              ch4_confluence_percore.csv   connectivity at link distances of 10 to 100 um, the proportion
                                           of stromal cells in niche areas and the mean niche size
              ch4_makeup_percore.csv       lineage of the niche cells and of all stromal cells
              ch4_nichearea_percore.csv    niche area at 15 um against the rest of the stroma:
                                           observed, same local density, label shuffle; all
                                           stromal cells and TrxR1-negative cells; the lineages
                                           of Supplementary Table 3b (TABLE4C_LINEAGES)
              ch4_cancer_percore.csv       cancer-cell TrxR1 (cytokeratin-positive cells of the
                                           Tumor category)
        s2  per patient, the case means
              ch4_confluence_percase.csv, ch4_makeup_percase.csv, ch4_nichearea_percase.csv,
              ch4_cancer_percase.csv
        s3  statistics
              supptable4a.csv      Supplementary Table 3a: niche-cell and stromal lineage means,
                                   Spearman rho
              supptable4a_nichearea.csv
                                   Supplementary Table 3a: the percentage of stromal cells in niche
                                   areas per patient, median (IQR) per cohort with LUSC and LUAD
                                   together (also the medians the Results paragraph on niche areas states)
              supptable4c.csv      Supplementary Table 3b: niche areas against the rest of the stroma,
                                   patient means and Wilcoxon p, LUSC and LUAD together (Figure 5b)
              supptable4e.csv      Supplementary Table 3c: the fraction, connectivity, T-cell and
                                   fibroblast proportions and density correlated with each other, the
                                   fraction against the niche-area proportion and the mean niche
                                   size, and connectivity against cancer-cell TrxR1
              fig4_stats.csv       every statistic of Figure 5a, b, e, f, g and h
              fig4ef_patients.csv  the patients of Figure 5e and 5f, one row each
              fig4ac_window.csv    the window of Figure 5a and 5c and the connectivity printed on c
              fig4d_examples.csv   the two example cores of Figure 5d with core and patient values
              suppfig6_stats.csv   the four cores of Supplementary Figure 6, checked against the
                                   selection rule its legend states
              ch4_cancer_trxr1.csv connectivity against cancer-cell TrxR1, LUSC (the Results paragraph
                                   on niche connectivity)
              ch4_numbers.csv      the text values that are a count, a maximum or a percentage of
                                   the values above
        s4  panels, PNG and PDF: fig4a_niches (Figure 5a), fig4b_bars (5b), fig4c_confluence (5c),
            fig4d_cores (5d), fig4e_histology (5e), fig4f_fraction (5f), fig4g_forest (5g),
            fig4h_sweep (5h), suppfig6_matrix (Supplementary Figure 6)
        s5  niche areas as densities (cells per 1,000 um2 of stroma), BOMI2
              ch4_density_percore.csv  per core, T-cell and fibroblast counts inside and outside the
                                       niche areas, and the stromal areas
              ch4_density_summary.csv  patient medians of the densities of Supplementary Table 3b
                                       and the text, and of their ratios, with Wilcoxon p
              supptable4c_density.csv  the printed cells of the density block of Supplementary
                                       Table 3b, cells per mm2 of stroma
        s6  niche measures and overall survival (Supplementary Table 3d), on the patients of Figure 5g
              ch4_niche_survival.csv   every Cox model: the fraction, connectivity, the niche-area
                                       proportion and the mean niche size alone; each niche measure
                                       beside the fraction; connectivity with the fraction and
                                       stromal density or the T-cell proportion; connectivity
                                       fully adjusted; per cohort, histology and term, and combined
              supptable3d.csv          the printed cells of Supplementary Table 3d
WHY:  one script with the shared niche and survival code; every table is written before the panels
      that draw it, and each panel reads its table back from disk.

Definitions:
  positions      raw export positions x redox.coords.MPP, computed on load by redox.cells.load_cells
                 (BOMI2 integer pixels x 0.4964 um per pixel, so that equal distances stay exactly
                 equal; BOMI1 as exported, in um)
  niche          a connected group of TrxR1-positive stromal cells linked at <= 30 um, any size
  niche area     every stromal cell within 15 um (inclusive) of a niche cell
  confluence     niche connectivity of the paper: cells of the largest niche / all TrxR1-positive
                 stromal cells of the core
  lineage        the lineage column of cells.tsv (priority CD3 > CD68 > CK > alpha-SMA > other; a cell
                 positive for more than one lineage marker is assigned to the first lineage in this
                 order)
  same density   the inside minus outside difference within value-binned fifths of local stromal
                 density (stromal cells within 30 um), weighted by the cells compared, fifths with
                 >= 5 cells on each side
  label shuffle  TrxR1 labels permuted within lineage, 20 times per core, the niche area taken as
                 the same proportion of cells nearest to the shuffled positive cells; net = observed
                 minus the shuffled mean. One generator, seeded SEED, consumed over the BOMI2 cores
                 and then the BOMI1 cores in the order of cells.tsv.
  case value     mean over the patient's qualifying cores weighted by stromal cells, each measure
                 over the cores that carry it; cancer-cell TrxR1 pooled over the cores with >= 20
                 cytokeratin-positive Tumor cells
  niche summaries
                 per core, the proportion of stromal cells in niche areas (niche cells and every
                 stromal cell within 15 um of one, over the stromal cells with a TrxR1 call) and
                 the mean niche size (TrxR1-positive stromal cells / niches, an isolated
                 TrxR1-positive cell a niche of one); case values as above
  models         Cox with Efron ties, HR per 1 SD within the analyzed set, stage as one linear term;
                 cohorts pooled by fixed effect; the histology interaction standardizes once over
                 LUSC and LUAD and pools the interaction terms of the two cohorts. A model with
                 several exposures (s6) standardizes each over the analyzed set and pools each
                 term over the two cohorts by fixed effect
  densities      BOMI2 only; every stromal cell of a core, a missing TrxR1 call read as negative (so
                 the core rule counts the few stromal cells without a call, which s1 leaves out); a
                 5 um grid over the core; the stroma is mask A, the grid points within 10 um of a
                 stromal cell; the niche area of a grid point is defined as for a cell. Per patient
                 the mean over cores, then medians over LUSC and LUAD patients of the inside and
                 outside values and of their ratio; units are cells per 1,000 um2 in the core tables
                 and cells per mm2 in the printed block

Seeds: SEED = 20260921 (the label shuffle of the niche-area control); SEED_BOOTSTRAP = 20260914 (the
bootstrap intervals of Figure 5e, one generator for the four cohort and histology groups in table
order); SEED_JITTER = 4 (horizontal jitter of the Figure 5e dots, drawing only).

RUN: python3 chapter4.py [--only s1 s2 s3 s4 s5 s6]   (about two minutes; s1, the label shuffle, and s5,
     the density grids, are most of it; s6 reads fig4_stats.csv of s3 and takes seconds)
"""
import argparse
import json
import time

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from redox import cells as cellsmod
from redox import clinical as clinicalmod
from redox import niches, paths, plotting, staining, stats, survival

CHAPTER = 4

SEED = 20260921              # the label shuffle of the niche-area control
SEED_BOOTSTRAP = 20260914    # bootstrap intervals of the Figure 5e medians
SEED_JITTER = 4              # horizontal jitter of the Figure 5e dots (drawing only)
# The paper shows "confluence" as "niche connectivity" (or "connectivity") and a "confluent" example as
# "connected". File, column and role names inside the package keep "confluence"/"confluent".
SHOWN_ROLE = {"confluent": "connected"}

RADII = (10, 15, 20, 30, 50, 75, 100)   # link distances of the confluence sweep, um
LINK_UM = 30                            # the niche link distance of the manuscript
ZONE_UM = 15.0                          # the niche-area radius, half the link distance
DENS_UM = 30.0                          # the radius of the local stromal density
MIN_SIDE = 20                           # cells needed inside and outside for a difference
MIN_BAND = 5                            # cells needed on each side within one density fifth
N_SHUFFLE = 20                          # label shuffles per core
MIN_CK = 20                             # cytokeratin-positive Tumor cells for cancer-cell TrxR1
BOOT_REPS = 2000
DENSITY_GRID_UM = 5.0                   # grid step of the stromal masks of s5
DENSITY_COHORT = "BOMI2"                # s5 covers the cohort the density statements are about
DENSITY_LINEAGES = ["T", "Fib", "Mac"]  # the lineages the density statements name 
# The densities of the density block of Supplementary Table 3b and of the text, as (feature, lineage, TrxR1
# states counted); lineage None counts every stromal cell (n_in and n_out of ch4_density_percore.csv), the
# first row of the block.
DENSITY_FEATURES = [("Fib_any", "Fib", ("pos", "neg")), ("T_neg", "T", ("neg",)), ("Fib_neg", "Fib", ("neg",)),
                    ("T_any", "T", ("pos", "neg")), ("Mac_any", "Mac", ("pos", "neg")), ("Mac_neg", "Mac", ("neg",)),
                    ("all_cells", None, ("pos", "neg"))]
# The density block of Supplementary Table 3b: its rows in order, as (feature, printed label), and
# the stromal mask: mask A, the grid points within 10 um of a stromal cell.
TABLE4C_DENSITY_ROWS = [("all_cells", "all stromal cells"), ("T_any", "T cells"),
                        ("T_neg", "TrxR1-negative T cells"), ("Fib_any", "fibroblasts"),
                        ("Fib_neg", "TrxR1-negative fibroblasts"), ("Mac_any", "macrophages"),
                        ("Mac_neg", "TrxR1-negative macrophages")]
DENSITY_PRINTED_MASK = "A"
PER_MM2 = 1000.0                        # cells per 1,000 um2 to cells per mm2 (1 mm2 = 1,000,000 um2)
MIN_CASES = 10                          # patients needed for a summary row

HISTOLOGIES = ("LUSC", "LUAD")

# The two simple per-core niche summaries: the proportion of stromal
# cells in niche areas (printed in Supplementary Tables 3c and 3d as "percentage of stromal cells in niche areas"; the
# code keeps "share" in its column and term names) and the mean niche size, columns of ch4_confluence_percore.csv and
# ch4_confluence_percase.csv.
NICHE_SUMMARIES = ["niche_area_share", "niche_mean_size"]
# Supplementary Table 3d (s6): every model as (exposures, extra covariates). Each exposure is per 1 SD of the
# analyzed set and every model is adjusted for stage; FULL_COVARIATES adds age, sex, WHO performance status and the
# TMA slide (one indicator per slide), the fully adjusted model of connectivity. The two cohorts are combined per term by fixed effect (survival.meta_fixed). In the
# tables a model is named by its exposures joined with "+" and by its adjustment, "stage" or "all".
FULL_COVARIATES = ["age", "male", "ps", "slide"]
FULL_LABEL = "stage, age, sex, WHO performance status and slide"
SURVIVAL_MODELS = [(["frac"], []), (["conf"], []), (["share"], []), (["size"], []),
                   (["frac", "conf"], []), (["frac", "share"], []), (["frac", "size"], []),
                   (["frac", "conf", "dens"], []), (["frac", "conf", "prop_T"], []),
                   (["conf"], FULL_COVARIATES)]
# The printed rows of Supplementary Table 3d, in its order, as (block, label, model, adjustment, term): the niche
# measures alone and beside the fraction, for both histologies and both cohorts; then the joint models, printed for
# LUSC with the two cohorts combined only (the LUAD and per-cohort values stay in
# ch4_niche_survival.csv). The fully adjusted connectivity row sits next to the stage-adjusted one. The labels are
# the row labels the supplement prints; the joint models sit under the block header
# "Joint models, squamous cell carcinoma, both cohorts combined".
TABLE3D_ROWS = [
    ("single", "TrxR1-positive fraction of stromal cells", "frac", "stage", "frac"),
    ("single", "niche connectivity (Figure 5g)", "conf", "stage", "conf"),
    ("full", f"niche connectivity, adjusted for {FULL_LABEL}", "conf", "all", "conf"),
    ("beside", "niche connectivity, beside the TrxR1-positive fraction", "frac+conf", "stage", "conf"),
    ("single", "percentage of stromal cells in niche areas", "share", "stage", "share"),
    ("beside", "percentage of stromal cells in niche areas, beside the TrxR1-positive fraction", "frac+share", "stage",
     "share"),
    ("single", "mean niche size", "size", "stage", "size"),
    ("beside", "mean niche size, beside the TrxR1-positive fraction", "frac+size", "stage", "size"),
]
TABLE3D_JOINT = [("fraction + connectivity", ["frac", "conf"]),
                 ("fraction + connectivity + stromal cell density", ["frac", "conf", "dens"]),
                 ("fraction + connectivity + T-cell proportion", ["frac", "conf", "prop_T"])]
TABLE3D_JOINT_GROUPS = [("LUSC", "both")]
TERM_LABEL = {"frac": "TrxR1-positive fraction", "conf": "connectivity", "share": "percentage of stromal cells in "
              "niche areas", "size": "mean niche size", "dens": "stromal cell density",
              "prop_T": "T-cell proportion of stromal cells"}

# The readouts of the niche-area comparison, in the order of Supplementary Table 3b.
AREA_KEYS = ["inall", "outall", "all", "dens", "shuf_net", "neg", "negdens", "negshuf_net"]
# The label of the niche measure in ch4_cancer_trxr1.csv (niche connectivity, biggest_30, patient mean).
CANCER_METRIC = "connectivity"
# The rows of Supplementary Table 3b (LUSC and LUAD together), in its order, by cohort and lineage. Uppsala I
# (BOMI1) has T cells and the CD3-negative, cytokeratin-negative cells (Other); Uppsala II (BOMI2) also has
# fibroblasts and macrophages.
TABLE4C_LINEAGES = {"BOMI2": ["T", "Fib", "Mac", "Other"], "BOMI1": ["T", "Other"]}
# The selection rule of the Supplementary Figure 6 legend: at least this many stromal cells in the core
# shown, and at least this TrxR1-positive fraction of them.
SUPP_MIN_CELLS = 1000
SUPP_MIN_CORE_FRACTION = 0.05

# Panel colors of this figure alone; lineage and histology colors come from redox.plotting.
BIG, OTHER, NEGC = "#e34948", "#f4b9b3", "#d6d5cf"
HALO, NEG_IN, NEG_OUT = "#fbe6e2", "#8f8e89", "#dcdbd5"
D_NEG, D_BIG, D_BIG_EDGE, D_OTH_EDGE = "#bdbcb6", "#e8645e", "#a12b2a", "#dd6a61"
D_TUM, D_TUM_EDGE = "#efe4d6", "#b08d68"
LUSC_LIGHT = "#9cc0ee"       # BOMI1 LUSC in panel f

# Supplementary Figure 6 cells as examples.json names them, and how the table labels them.
SUPP_CELLS = {"LL": ("low", "low"), "LH": ("low", "high"), "HL": ("high", "low"), "HH": ("high", "high")}


# ---------------------------------------------------------------- shared helpers

def _out(name):
    return paths.output(CHAPTER, name)


def _read(name):
    """A table this chapter wrote, re-read from disk with identifiers as strings."""
    path = _out(name)
    if not path.exists():
        raise FileNotFoundError(f"{path} is missing. Run the section that writes it first.")
    return pd.read_csv(path, dtype={"Case": str, "core": str, "cohort": str, "patient": str},
                       float_precision="round_trip")


def _write(frame, name, what, columns):
    """Write a table and its sibling .md description of every column. Returns the path."""
    path = _out(name)
    frame.to_csv(path, index=False)
    lines = [f"# {name}", "", what, "", "| column | meaning |", "|---|---|"]
    lines += [f"| {c} | {m} |" for c, m in columns] + [""]
    _out(name.rsplit(".", 1)[0] + ".md").write_text("\n".join(lines))
    print(f"   wrote {path.name} ({len(frame)} rows)")
    return path


def _merge(left, right, left_name, right_name, **kw):
    """pandas merge that stops the run when it keeps no row, naming the two tables."""
    m = left.merge(right, **kw)
    if len(m) == 0:
        raise ValueError(f"merging {left_name} with {right_name} kept no patient")
    return m


def _names(cohort):
    """Lineages of a cohort in the priority order of cells.tsv."""
    return list(cellsmod.LINEAGES[cohort])


def _clinical(cohort):
    return clinicalmod.load_clinical(cohort)


def _fmt_p(p):
    return "p < 0.001" if float(p) < 0.001 else f"p = {float(p):.3f}"


def _comp(lin, a, b, names):
    """Lineage proportions of set a minus those of set b, one value per lineage."""
    return np.array([(lin[a] == nm).mean() - (lin[b] == nm).mean() for nm in names])


def _banded(lin, a, b, q, names):
    """The a minus b difference within density fifths, weighted by the cells compared."""
    acc, w = np.zeros(len(names)), 0
    for v in np.unique(q):
        i, o = a & (q == v), b & (q == v)
        if i.sum() >= MIN_BAND and o.sum() >= MIN_BAND:
            n = int(i.sum() + o.sum())
            acc += n * _comp(lin, i, o, names)
            w += n
    return acc / w if w else np.full(len(names), np.nan)


def _density_fifths(xy):
    """Value-binned fifths of local stromal density, cut on all stromal cells of the core."""
    dens = cKDTree(xy).query_ball_point(xy, DENS_UM, return_length=True) - 1
    return np.digitize(dens, np.unique(np.quantile(dens, [.2, .4, .6, .8])), right=True)


def _stroma_and_tumor(cohort):
    """Stromal cells with a TrxR1 call and a position, and Tumor cells, from cells.tsv."""
    cols = ["core", "Case", "tissue", "x_um", "y_um", "lineage", "TrxR1", "CK"]
    c = cellsmod.load_cells(cohort, columns=cols)
    stroma = c[(c.tissue == "Stroma") & c.TrxR1.notna()].dropna(subset=["x_um", "y_um"])
    tumor = c[c.tissue == "Tumor"]
    return stroma, tumor


def _examples():
    """The example-core file, or none with a message when it is absent."""
    try:
        path = paths.examples()
    except KeyError as exc:
        print(f"   examples: {exc}; example panels and rows skipped")
        return None
    if not path.exists():
        print(f"   examples: {path} not found; example panels and rows skipped")
        return None
    return json.loads(path.read_text())


# ---------------------------------------------------------------- s1 per-core tables

def _core_area(xy, trx, lin, names, labels, q, rng):
    """Niche area at ZONE_UM against the rest of the stroma, with its two controls.

    Returns the per-core values and the niche-area membership of every cell. The shuffle consumes
    rng in a fixed order (lineage order within a shuffle, N_SHUFFLE shuffles, only on cores that
    pass the observed gate), so that its values reproduce. With niches
    of any size every TrxR1-positive cell is a niche cell, so the shuffled niche cells are
    the shuffled positive cells themselves and need no linking, and no positive cell lies outside
    the niche area (checked).
    """
    N = len(trx)
    memb = niches.area_membership(xy, trx, labels, ZONE_UM)
    inside = memb >= 0
    if (trx & ~inside).any():
        raise ValueError("a TrxR1-positive stromal cell lies outside the niche area")
    r = {}
    if inside.sum() >= MIN_SIDE and (~inside).sum() >= MIN_SIDE:
        v = {"all": _comp(lin, inside, ~inside, names), "dens": _banded(lin, inside, ~inside, q, names),
             "inall": np.array([(lin[inside] == nm).mean() for nm in names]),
             "outall": np.array([(lin[~inside] == nm).mean() for nm in names])}
        ni, no = inside & ~trx, ~inside & ~trx
        okneg = ni.sum() >= MIN_SIDE and no.sum() >= MIN_SIDE
        v["neg"] = _comp(lin, ni, no, names) if okneg else np.full(len(names), np.nan)
        v["negdens"] = _banded(lin, ni, no, q, names) if okneg else np.full(len(names), np.nan)
        p_in = inside.mean()
        p_in_neg = ni.sum() / (~trx).sum()
        sh, shn = [], []
        for _ in range(N_SHUFFLE):
            t2 = trx.copy()
            for nm in names:
                m = lin == nm
                t2[m] = rng.permutation(trx[m])
            if not t2.any():
                continue
            dd = cKDTree(xy[t2]).query(xy)[0]
            order = np.argsort(dd, kind="stable")
            a = np.zeros(N, bool)
            a[order[:int(round(p_in * N))]] = True
            sh.append(_comp(lin, a, ~a, names))
            if okneg:
                pool = np.flatnonzero(~t2)
                oo = pool[np.argsort(dd[pool], kind="stable")]
                b = np.zeros(N, bool)
                b[oo[:int(round(p_in_neg * len(pool)))]] = True
                shn.append(_comp(lin, b, ~t2 & ~b, names))
        v["shuf_net"] = v["all"] - np.mean(sh, axis=0) if sh else np.full(len(names), np.nan)
        v["negshuf_net"] = v["neg"] - np.mean(shn, axis=0) if shn else np.full(len(names), np.nan)
        for key, arr in v.items():
            for nm, x in zip(names, arr):
                r[f"{key}_{nm}"] = float(x)
    return r, memb


def _printed(values, lineages):
    """The readouts <key>_<lineage> of one core that a printed table reads: every readout of the listed
    lineages. The values are computed over every lineage (the label shuffle draws for each), so leaving the
    others out changes no kept value."""
    keep = {}
    for col, x in values.items():
        key, nm = col.rsplit("_", 1)
        if nm in lineages:
            keep[col] = x
    return keep


def s1_percore():
    """Every per-core niche measure, both cohorts, the core rule applied once.

    One pass over the qualifying cores computes all four tables. The label-shuffle generator is
    created once, before the BOMI2 cores, and consumed by the niche-area section alone.
    """
    rng = np.random.default_rng(SEED)
    rows = {k: [] for k in ("confluence", "makeup", "nichearea", "cancer")}
    for cohort in paths.COHORTS:
        names = _names(cohort)
        stroma, tumor = _stroma_and_tumor(cohort)
        ok = cellsmod.qualifying_cores(stroma, *cellsmod.RULE_CH4)
        ck = tumor[tumor.TrxR1.notna() & (tumor.CK.fillna(0) == 1)]
        n_ck = ck.groupby("core").size()
        pos_ck = ck.groupby("core").TrxR1.sum()
        t0, done = time.time(), 0
        for core, s in stroma.groupby("core", sort=False):
            if core not in ok:
                continue
            trx = s.TrxR1.to_numpy(dtype=int).astype(bool)
            lin = s.lineage.to_numpy(dtype=object)
            xy = s[["x_um", "y_um"]].to_numpy(dtype=float)
            base = dict(cohort=cohort, core=core, Case=str(s.Case.iloc[0]), n_stroma=len(s),
                        n_pos=int(trx.sum()))
            # confluence at every link distance
            conf = dict(base)
            labels = None
            for rad in RADII:
                lab = niches.link(xy, trx, float(rad))
                conf[f"biggest_{rad}"] = niches.confluence(lab, trx)
                if rad == LINK_UM:
                    labels = lab
            # lineage of niche cells (every positive cell) and of all stromal cells
            r = dict(base)
            for nm in names:
                r[f"mk_pos_{nm}"] = float((lin[trx] == nm).mean())
                r[f"mk_all_{nm}"] = float((lin == nm).mean())
            rows["makeup"].append(r)
            # niche area against the rest
            q = _density_fifths(xy)
            area, memb = _core_area(xy, trx, lin, names, labels, q, rng)
            rows["nichearea"].append({**base, **_printed(area, TABLE4C_LINEAGES[cohort])})
            # the two simple niche summaries of Supplementary Tables 3c and 3d, from the same
            # niches (30 um links, an isolated TrxR1-positive cell a niche of one) and niche areas (15 um)
            conf["niche_area_share"] = float((memb >= 0).mean())
            conf["niche_mean_size"] = float(niches.sizes(labels).mean())
            rows["confluence"].append(conf)
            # cancer-cell TrxR1 of the core
            rows["cancer"].append({**base, "n_ck": int(n_ck.get(core, 0)), "pos_ck": int(pos_ck.get(core, 0))})
            done += 1
            if done % 100 == 0:
                print(f"   s1 {cohort}: {done} cores, {time.time() - t0:.0f} s", flush=True)
        print(f"s1 {cohort}: {done} qualifying cores of {stroma.core.nunique()}, {time.time() - t0:.0f} s",
              flush=True)
    _write_percore(rows)


def _lineage_note():
    return ("lineage suffix: T, Mac, CK, Fib, Other in BOMI2; T, CK, Other in BOMI1, where Other holds the "
            "CD3-negative, cytokeratin-negative cells; a lineage a cohort lacks is empty")


def _area_note():
    return ("lineage suffix: the lineages Supplementary Table 3b prints, T, Fib, Mac, Other in BOMI2 and T, Other "
            "in BOMI1 (Other: the CD3-negative, cytokeratin-negative cells), every readout; a lineage or "
            "readout not written for a cohort is empty")


def _write_percore(rows):
    base = [("cohort", "BOMI2 or BOMI1"), ("core", "core identifier (inForm Sample Name)"),
            ("Case", "patient identifier"), ("n_stroma", "stromal cells of the core with a TrxR1 call"),
            ("n_pos", "of those, TrxR1-positive")]
    rule = ("One row per core passing RULE_CH4 (>= 60 stromal cells with a TrxR1 call, >= 25 positive, "
            ">= 25 negative), BOMI2 cores then BOMI1 cores in the order of cells.tsv.")
    _write(pd.DataFrame(rows["confluence"]), "ch4_confluence_percore.csv",
           rule + " Confluence at each link distance, and the two simple niche summaries at the manuscript's "
           "link distance (Supplementary Tables 3c and 3d).",
           base + [(f"biggest_{r}", f"confluence with links at <= {r} um: cells in the largest group / "
                    "all TrxR1-positive stromal cells") for r in RADII]
           + [("niche_area_share", "stromal cells within 15 um (inclusive) of a niche cell, the niche cells "
               "included, / all stromal cells of the core with a TrxR1 call (the niche areas of s1)"),
              ("niche_mean_size", f"TrxR1-positive stromal cells / niches, niches linked at <= {LINK_UM} um, an "
               "isolated TrxR1-positive cell counting as a niche of one")])
    _write(pd.DataFrame(rows["makeup"]), "ch4_makeup_percore.csv",
           rule + " Lineage of the niche cells (every TrxR1-positive stromal cell) and of all stromal cells.",
           base + [("mk_pos_<lineage>", "proportion of the TrxR1-positive stromal cells in this lineage"),
                   ("mk_all_<lineage>", "proportion of all stromal cells in this lineage"),
                   ("(lineages)", _lineage_note())])
    _write(pd.DataFrame(rows["nichearea"]), "ch4_nichearea_percore.csv",
           rule + " Niche area (every stromal cell within 15 um of a niche cell, niches linked at 30 um, "
           "any size) against the rest of the stroma of the core. The difference columns exist only "
           "for cores with >= 20 stromal cells inside and outside; the TrxR1-negative ones need >= 20 "
           "TrxR1-negative cells on each side. Values are fractions, not percent.",
           base + [("inall_<lineage>", "proportion of the stromal cells inside the niche area"),
                   ("outall_<lineage>", "proportion of the stromal cells outside"),
                   ("all_<lineage>", "inside minus outside, all stromal cells"),
                   ("dens_<lineage>", "the same within fifths of local stromal density (cells within "
                    "30 um), weighted, fifths with >= 5 cells on each side"),
                   ("shuf_net_<lineage>", "all minus its mean over 20 shuffles of the TrxR1 label within "
                    "lineage, the niche area matched in size (seed SEED)"),
                   ("neg_<lineage>", "inside minus outside, TrxR1-negative stromal cells"),
                   ("negdens_<lineage>", "the same within density fifths"),
                   ("negshuf_net_<lineage>", "neg minus its shuffled mean, the area matched in size "
                    "among the TrxR1-negative cells"),
                   ("(lineages)", _area_note())])
    _write(pd.DataFrame(rows["cancer"]), "ch4_cancer_percore.csv",
           rule + " Cancer-cell TrxR1 of the core.",
           base + [("n_ck", "cytokeratin-positive cells of the Tumor category with a TrxR1 call"),
                   ("pos_ck", "of those, TrxR1-positive")])


# ---------------------------------------------------------------- s2 per-patient tables

def _wavg(g, c):
    ok = g[c].notna()
    return float(np.average(g.loc[ok, c], weights=g.loc[ok, "n_stroma"])) if ok.any() else np.nan


def _percase(pc, cols, extra=None):
    """Case means weighted by stromal cells, each column over the cores that carry it."""
    def one(g):
        out = {c: _wavg(g, c) for c in cols}
        if extra:
            out.update(extra(g))
        return pd.Series(out)
    return pc.groupby(["cohort", "Case"]).apply(one, include_groups=False).reset_index()


def _with_histology(frame, name, keep=None):
    """Attach the histology of clinical.csv per cohort; keep only the listed histologies if given."""
    cl = pd.concat([_clinical(c)[["Case", "histology"]].assign(cohort=c) for c in paths.COHORTS])
    m = _merge(frame, cl, name, "clinical.csv", on=["cohort", "Case"])
    if keep:
        m = m[m.histology.isin(keep)]
    return m.reset_index(drop=True)


def _value_cols(pc, keys):
    return [c for c in pc.columns if c.rsplit("_", 1)[0] in keys and c.rsplit("_", 1)[1] in
            set(cellsmod.LINEAGES["BOMI2"])]


def s2_percase():
    """Per-patient means of every per-core table (case values: see the header)."""
    base = [("cohort", "BOMI2 or BOMI1"), ("Case", "patient identifier"),
            ("histology", "LUSC, LUAD or other, from clinical.csv")]
    how = "Mean over the patient's qualifying cores weighted by stromal cells, each column over the cores that carry it."

    conf = _read("ch4_confluence_percore.csv")
    cols = [f"biggest_{r}" for r in RADII] + NICHE_SUMMARIES
    t = _percase(conf, cols, lambda g: {"n_stroma": float(g.n_stroma.sum())})
    _write(_with_histology(t, "ch4_confluence_percase.csv"), "ch4_confluence_percase.csv",
           "One row per patient with a qualifying core. " + how,
           base + [(f"biggest_{r}", f"confluence at <= {r} um links") for r in RADII]
           + [("niche_area_share", "proportion of the stromal cells in niche areas (Supplementary Tables 3c and 3d)"),
              ("niche_mean_size", "TrxR1-positive stromal cells per niche, an isolated one a niche of one "
               "(Supplementary Tables 3c and 3d)"),
              ("n_stroma", "stromal cells over the qualifying cores")])

    mk = _read("ch4_makeup_percore.csv")
    cols = [c for c in mk.columns if c.startswith("mk_")]
    t = _percase(mk, cols)
    _write(_with_histology(t, "ch4_makeup_percase.csv"), "ch4_makeup_percase.csv",
           "One row per patient with a qualifying core. " + how,
           base + [("mk_pos_<lineage>", "proportion of the niche cells in this lineage"),
                   ("mk_all_<lineage>", "proportion of all stromal cells in this lineage"),
                   ("(lineages)", _lineage_note())])

    area = _read("ch4_nichearea_percore.csv")
    cols = _value_cols(area, AREA_KEYS)
    t = _percase(area, cols)
    _write(_with_histology(t, "ch4_nichearea_percase.csv", keep=HISTOLOGIES), "ch4_nichearea_percase.csv",
           "One row per LUSC or LUAD patient with a qualifying core. " + how + " Fractions, not percent.",
           base + [("<readout>_<lineage>", "the per-core readouts of ch4_nichearea_percore.csv"),
                   ("(lineages)", _area_note())])

    ca = _read("ch4_cancer_percore.csv")

    def pooled(g):
        use = g.n_ck >= MIN_CK
        if not use.any():
            return {"tum_ck": np.nan}
        return {"tum_ck": float(g.loc[use, "pos_ck"].sum() / g.loc[use, "n_ck"].sum())}
    t = _percase(ca, [], pooled)
    _write(_with_histology(t, "ch4_cancer_percase.csv"), "ch4_cancer_percase.csv",
           "One row per patient with a qualifying core.",
           base + [("tum_ck", "TrxR1-positive fraction of the cytokeratin-positive Tumor cells, pooled over the "
                    "qualifying cores with >= 20 such cells; empty when no core has them")])


# ---------------------------------------------------------------- s3 statistics

def _survival_set(cohort):
    """Figure 5e to 5h: connectivity, the stromal TrxR1 fraction and follow-up, LUSC and LUAD.

    The TrxR1 fraction comes from percase.csv (mif/03_percase.py), the patients with >= 20 stromal
    cells, and the clinical columns from clinical.csv; the rows keep the order of the confluence
    table, which is the order the bootstrap draws from.
    """
    conf = _read("ch4_confluence_percase.csv")
    conf = conf[conf.cohort == cohort].drop(columns=["histology", "n_stroma"])
    pc = pd.read_csv(paths.intermediate(cohort, "percase.csv"), dtype={"Case": str},
                     usecols=["Case", "n_stroma", "stromal_trxr1"], float_precision="round_trip")
    pc = pc[pc.n_stroma >= cellsmod.MIN_STROMAL_CELLS].rename(columns={"stromal_trxr1": "frac"})
    d = _merge(conf, pc[["Case", "frac"]], "ch4_confluence_percase.csv", f"{cohort} percase.csv", on="Case")
    d = _merge(d, _clinical(cohort), "ch4_confluence_percase.csv", f"{cohort} clinical.csv", on="Case")
    d = d[d.histology.isin(HISTOLOGIES) & d.OS_weeks.notna() & d.event.notna() & (d.OS_weeks > 0)].copy()
    d["conf"] = d[f"biggest_{LINK_UM}"]
    return d.reset_index(drop=True)


def _table4a(mk):
    rows = []
    for cohort in paths.COHORTS:
        g = mk[(mk.cohort == cohort) & mk.histology.isin(HISTOLOGIES)]
        for nm in _names(cohort):
            rho, p = stats.spearman(g[f"mk_pos_{nm}"], g[f"mk_all_{nm}"])
            rows.append(dict(cohort=cohort, lineage=nm, niche_pct=100 * g[f"mk_pos_{nm}"].mean(),
                             all_pct=100 * g[f"mk_all_{nm}"].mean(), rho=rho, p=p, n=len(g)))
    return pd.DataFrame(rows)


def _table4a_nichearea(mk):
    """The niche-area row of Supplementary Table 3a and the range the Results paragraph on niche
    areas states: the percentage of stromal cells in niche areas per patient (niche_area_share of
    ch4_confluence_percase.csv), median and interquartile range per cohort with LUSC and LUAD together (the printed
    cell). The run stops when the patients are not those of Table 3a (mk, ch4_makeup_percase.csv)."""
    conf = _read("ch4_confluence_percase.csv")
    rows = []
    for cohort in paths.COHORTS:
        c = conf[(conf.cohort == cohort) & conf.histology.isin(HISTOLOGIES)]
        t3a = mk[(mk.cohort == cohort) & mk.histology.isin(HISTOLOGIES)]
        if sorted(c.Case) != sorted(t3a.Case) or c.niche_area_share.isna().any():
            raise ValueError(f"{cohort}: the patients of ch4_confluence_percase.csv with a niche-area percentage "
                             "are not those of Supplementary Table 3a (ch4_makeup_percase.csv)")
        v = c.niche_area_share
        q1, med, q3 = v.quantile(.25), v.median(), v.quantile(.75)
        rows.append(dict(cohort=cohort, histology="both", n=len(v), median=med, q1=q1, q3=q3,
                         text=f"{100 * med:.0f} % ({100 * q1:.0f} to {100 * q3:.0f} %)"))
    return pd.DataFrame(rows)


def _table4c(area):
    """The rows of Supplementary Table 3b.

    Every readout has its patient mean; the differences have their Wilcoxon p; the patient count
    is kept for the two readouts whose row groups the table heads with it (all, neg).
    """
    rows = []
    for cohort, names in TABLE4C_LINEAGES.items():
        z = area[(area.cohort == cohort) & area.histology.isin(HISTOLOGIES)]
        for nm in names:
            r = dict(cohort=cohort, lineage=nm)
            for k in AREA_KEYS:
                v = z[f"{k}_{nm}"].dropna().to_numpy()
                r[f"{k}_mean"] = 100 * v.mean()
                if k in ("all", "neg"):
                    r[f"{k}_n"] = len(v)
                if k not in ("inall", "outall"):
                    r[f"{k}_p"] = stats.wilcoxon_p(v)
            rows.append(r)
    return pd.DataFrame(rows).astype({"all_n": "Int64", "neg_n": "Int64"})   # counts stay integers


def _table4e():
    """Supplementary Table 3c and the Results paragraph on niche connectivity: the
    correlations.

    Spearman rho per cohort and histology among the TrxR1-positive fraction of stromal cells
    (all primary cores, as Figure 1), confluence (biggest_30 of ch4_confluence_percase.csv, the
    case mean over qualifying cores, as Figure 5f and 5g), the T-cell and fibroblast proportions
    of all stromal cells on the qualifying cores (ch4_makeup_percase.csv, as Supplementary
    Table 3a), and stromal cell density (stromal cells per mm2 of stroma, all primary cores).
    Patients: LUSC or LUAD with follow-up (OS_weeks > 0), at least MIN_STROMAL_CELLS stromal
    cells, at least one qualifying core. Three of the pairs are the fraction against the
    proportion of stromal cells in niche areas (share) and against the mean
    niche size (size), both case means over the qualifying cores as confluence, and confluence
    against the TrxR1-positive fraction of cancer cells (tum, tum_ck of ch4_cancer_percase.csv, the
    measure of ch4_cancer_trxr1.csv), on the same patients.
    """
    pairs = [("frac", "conf"), ("frac", "T"), ("frac", "Fib"), ("conf", "T"), ("conf", "Fib"),
             ("T", "Fib"), ("conf", "dens"), ("frac", "dens"), ("frac", "share"), ("frac", "size"),
             ("conf", "tum")]
    conf = _read("ch4_confluence_percase.csv")
    mk = _read("ch4_makeup_percase.csv")
    ca = _read("ch4_cancer_percase.csv")
    rows = []
    for cohort in paths.COHORTS:
        names = _names(cohort)
        pc = pd.read_csv(paths.intermediate(cohort, "percase.csv"), dtype={"Case": str},
                         float_precision="round_trip")
        pc = pc.rename(columns={"stromal_trxr1": "frac"})
        pc["dens"] = pc[[f"dens_{nm}_S" for nm in names]].sum(axis=1)
        c = conf[conf.cohort == cohort][["Case", f"biggest_{LINK_UM}"] + NICHE_SUMMARIES].rename(
            columns={f"biggest_{LINK_UM}": "conf", "niche_area_share": "share", "niche_mean_size": "size"})
        mcols = ["Case", "mk_all_T"] + (["mk_all_Fib"] if cohort != "BOMI1" else [])
        m = mk[mk.cohort == cohort][mcols].rename(columns={"mk_all_T": "T", "mk_all_Fib": "Fib"})
        t = ca[ca.cohort == cohort][["Case", "tum_ck"]].rename(columns={"tum_ck": "tum"})
        d = (pc[["Case", "n_stroma", "frac", "dens"]].merge(c, on="Case").merge(m, on="Case")
             .merge(t, on="Case", how="left"))
        d = _merge(d, _clinical(cohort)[["Case", "histology", "OS_weeks"]], "ch4_table4e patients",
                  f"{cohort} clinical.csv", on="Case")
        d = d[(d.n_stroma >= cellsmod.MIN_STROMAL_CELLS) & d.histology.isin(HISTOLOGIES)
             & d.OS_weeks.notna() & (d.OS_weeks > 0)]
        for h in HISTOLOGIES:
            z = d[d.histology == h]
            for a, b in pairs:
                if cohort == "BOMI1" and "Fib" in (a, b):
                    rows.append(dict(cohort=cohort, histology=h, pair=f"{a}~{b}", n=len(z),
                                     rho=np.nan, p=np.nan))
                    continue
                u = z[[a, b]].dropna()
                rho, p = stats.spearman(u[a], u[b])
                rows.append(dict(cohort=cohort, histology=h, pair=f"{a}~{b}", n=len(u), rho=rho, p=p))
    return pd.DataFrame(rows)


def _core_geometry(stroma, tumor, cohort, case, key, conf_pc):
    """One example core: centered positions, TrxR1 flags, niche labels and the largest niche.

    The core is the one of this patient whose identifier contains key (a short core-name fragment
    from examples.json); it must be unique within the patient, and its confluence must equal the
    value of ch4_confluence_percore.csv, so that the drawing shows the core the tables use.
    """
    row = conf_pc[(conf_pc.cohort == cohort) & (conf_pc.Case == str(case))
                  & conf_pc.core.str.contains(key, regex=False)]
    if len(row) != 1:
        raise ValueError(f"{cohort} patient {case}: {len(row)} qualifying cores match {key}, expected one")
    row = row.iloc[0]
    g = stroma[stroma.core == row.core]
    xy0 = g[["x_um", "y_um"]].to_numpy(dtype=float)
    off = (xy0.max(0) + xy0.min(0)) / 2
    xy = xy0 - off
    trx = g.TrxR1.to_numpy(dtype=int).astype(bool)
    labels = niches.link(xy, trx, float(LINK_UM))
    big_label = int(np.argmax(niches.sizes(labels)))
    big = trx & (labels == big_label)
    share = niches.confluence(labels, trx)
    if len(g) != row.n_stroma or abs(share - row[f"biggest_{LINK_UM}"]) > 1e-9:
        raise ValueError(f"{row.core}: recomputed confluence {share} does not match the per-core table")
    tu = tumor[tumor.core == row.core].dropna(subset=["x_um", "y_um"])
    return dict(cohort=cohort, Case=str(case), core=row.core, xy=xy, trx=trx, labels=labels, big=big,
                share=share, rate=float(trx.mean()), n=len(g),
                tumor=tu[["x_um", "y_um"]].to_numpy(dtype=float) - off)


def _window(geo, width):
    """Figure 5a and 5c: a square window centered on the largest niche, which it must hold whole.

    Returns the window center and the niche-area membership of every stromal cell of the core.
    """
    xy, trx, big, labels = geo["xy"], geo["trx"], geo["big"], geo["labels"]
    c0 = xy[big].mean(0)
    inwin = (np.abs(xy - c0) <= width / 2).all(1)
    in_area = niches.area_membership(xy, trx, labels, ZONE_UM) >= 0
    if int((inwin & big).sum()) != int(big.sum()):
        raise ValueError(f"{geo['core']}: the {width} um window does not hold the whole largest niche")
    return c0, in_area


def s3_statistics():
    """Supplementary Tables 3a to 3c, the Figure 5 statistics, the example cores and the text values."""
    mk = _read("ch4_makeup_percase.csv")
    area = _read("ch4_nichearea_percase.csv")

    t4a = _table4a(mk)
    _write(t4a, "supptable4a.csv",
           "Supplementary Table 3a and the Results paragraph on niche composition: lineage "
           "of the niche "
           "cells and of all stromal cells, "
           "mean over LUSC and LUAD patients with a qualifying core, and their Spearman correlation.",
           [("cohort", "BOMI2 or BOMI1"), ("lineage", "T, Mac, Fib, Other, CK (BOMI1: T, CK, Other)"),
            ("niche_pct", "mean over patients of the proportion of niche cells in this lineage, percent"),
            ("all_pct", "the same for all stromal cells, percent"),
            ("rho", "Spearman rho between the two proportions, one patient per point"),
            ("p", "its p value"), ("n", "patients")])

    _write(_table4a_nichearea(mk), "supptable4a_nichearea.csv",
           "Supplementary Table 3a, niche-area row, and the Results paragraph on niche areas "
           "(the medians it states): the percentage of stromal cells that lay in niche areas, per "
           "patient (niche cells and every stromal cell within 15 um of one, over the stromal cells with a TrxR1 "
           "call, mean over the patient's qualifying cores weighted by stromal cells; niche_area_share of "
           "ch4_confluence_percase.csv), on the patients of Table 3a: the printed cell, median (IQR) per cohort "
           "with LUSC and LUAD together.",
           [("cohort", "BOMI2 or BOMI1"), ("histology", "both (LUSC and LUAD together, the printed cell)"),
            ("n", "patients"), ("median", "median over patients, a fraction"),
            ("q1, q3", "25th and 75th percentiles, a fraction"),
            ("text", "the printed cell, whole percent")])

    t4c = _table4c(area)
    _write(t4c, "supptable4c.csv",
           "Supplementary Table 3b, Figure 5b and the Results paragraph on niche "
           "composition: niche areas against the rest of the stroma, mean over LUSC and LUAD patients, percent or "
           "percentage points; BOMI2 T, Fib, Mac, Other and BOMI1 T, Other (Other: the CD3-negative, "
           "cytokeratin-negative cells).",
           [("cohort", "BOMI2 or BOMI1"), ("lineage", "T, Fib, Mac, Other (BOMI1: T, Other)"),
            ("inall_mean, outall_mean", "percent of the stromal cells inside and outside the niche areas"),
            ("all_mean", "inside minus outside, all stromal cells, percentage points"),
            ("dens_mean", "the same at the same local density"),
            ("shuf_net_mean", "the same minus its label-shuffled expectation"),
            ("neg_mean, negdens_mean, negshuf_net_mean", "the three readouts for TrxR1-negative cells"),
            ("all_n, neg_n", "patients with the readout, all stromal cells and TrxR1-negative cells"),
            ("<readout>_p", "Wilcoxon signed-rank p of the patient values against zero, for the six "
             "differences")])

    t4e = _table4e()
    _write(t4e, "supptable4e.csv",
           "Supplementary Table 3c and the Results paragraph on niche connectivity: "
           "Spearman rho per "
           "cohort and histology among "
           "the TrxR1-positive fraction of stromal cells, confluence, the T-cell and fibroblast "
           "proportions of stromal cells and stromal cell density; also the fraction against "
           "the proportion of stromal cells in niche areas and the mean niche size, and confluence against the "
           "TrxR1-positive fraction of cancer cells.",
           [("cohort", "BOMI2 or BOMI1"), ("histology", "LUSC or LUAD"),
            ("pair", "frac (TrxR1-positive fraction, all cores), conf (confluence, biggest_30, "
             "qualifying cores), T, Fib (proportion of all stromal cells, qualifying cores; "
             "BOMI1 has no Fib), dens (stromal cells per mm2 of stroma, all cores), share (proportion of the "
             "stromal cells in niche areas, qualifying cores), size (TrxR1-positive cells per niche, an "
             "isolated one a niche of one, qualifying cores), tum (TrxR1-positive fraction of the "
             "cytokeratin-positive Tumor cells, qualifying cores with >= 20 such cells)"),
            ("n", "patients with both values (BOMI1 Fib pairs: all patients of the histology, "
             "rho and p empty)"),
            ("rho", "Spearman rho"), ("p", "its p value")])

    _fig4_stats(area)
    _cancer()
    _numbers()


def _fig4_stats(area):
    """Every statistic Figure 5 draws or its legend prints, in the row layout of fig4_stats.csv."""
    D = {c: _survival_set(c) for c in paths.COHORTS}
    stat = []

    def fit(h, var):
        out = {c: survival.hr_per_sd(D[c][D[c].histology == h], var, covars=("stage",)) for c in paths.COHORTS}
        out["both"] = survival.meta_fixed([out["BOMI2"], out["BOMI1"]])
        return out

    F = {h: fit(h, "conf") for h in HISTOLOGIES}
    prof = {(h, r): fit(h, f"biggest_{r}") for h in HISTOLOGIES for r in RADII}
    inter = {c: survival.histology_interaction(D[c], "conf", covars=("stage",)) for c in paths.COHORTS}
    inter["both"] = survival.meta_fixed([inter["BOMI2"], inter["BOMI1"]])
    for h in HISTOLOGIES:
        for k in ("BOMI2", "BOMI1", "both"):
            r = F[h][k]
            stat.append(dict(panel="g", group=f"{h} {k}", stat="HR per SD", value=r["HR"], lo=r["lo"],
                             hi=r["hi"], p=r["p"], n=r["n"], events=r["events"]))
            if h == "LUAD" and k != "both":         # panel h draws LUAD for both cohorts only
                continue
            for rad in RADII:
                q = prof[(h, rad)][k]
                stat.append(dict(panel="h", group=f"{h} {k} {rad} um", stat="HR per SD", value=q["HR"],
                                 lo=q["lo"], hi=q["hi"], p=q["p"], n=q["n"], events=q["events"]))
    for c in ("BOMI2", "BOMI1", "both"):
        stat.append(dict(panel="g", group=f"interaction {c}", stat="Wald p",
                         p=inter[c]["p_wald"] if c != "both" else inter[c]["p"]))

    ex = _examples()
    d_rows, s_rows = [], []
    if ex is None:
        for name in ("fig4ac_window.csv", "fig4d_examples.csv", "suppfig6_stats.csv"):
            for p in (_out(name), _out(name.rsplit(".", 1)[0] + ".md")):
                if p.exists():
                    p.unlink()
                    print(f"   examples: removed stale {p.name}")
    if ex is not None:
        conf_pc = _read("ch4_confluence_percore.csv")
        cells = {c: _stroma_and_tumor(c) for c in {e["cohort"] for e in ex["fig4d"] + ex["suppfig6"]}
                 | {ex["fig4a_window"]["cohort"]}}
        p2 = D["BOMI2"][D["BOMI2"].histology == "LUSC"].set_index("Case")
        pcase = {c: pd.read_csv(paths.intermediate(c, "percase.csv"), dtype={"Case": str},
                                usecols=["Case", "stromal_trxr1"], float_precision="round_trip")
                 .set_index("Case").stromal_trxr1 for c in paths.COHORTS}
        cpc = _read("ch4_confluence_percase.csv").set_index(["cohort", "Case"])
        for e in ex["fig4d"]:
            geo = _core_geometry(*cells[e["cohort"]], e["cohort"], e["Case"], e["core"], conf_pc)
            case = str(e["Case"])
            d_rows.append(dict(role=e["role"], cohort=e["cohort"], patient=case, core=geo["core"],
                               stromal_cells=geo["n"], core_TrxR1_fraction=geo["rate"], core_confluence=geo["share"],
                               patient_TrxR1_fraction=p2.loc[case, "frac"], patient_confluence=p2.loc[case, "conf"]))
        w = ex["fig4a_window"]
        geo = _core_geometry(*cells[w["cohort"]], w["cohort"], w["Case"], w["core"], conf_pc)
        _window(geo, float(w["width_um"]))                   # stops the run when the window cuts the niche
        stat.append(dict(panel="a", group="window", stat="window side, um", value=float(w["width_um"]), n=1))
        _write(pd.DataFrame([dict(cohort=w["cohort"], patient=str(w["Case"]), core=geo["core"],
                                  width_um=float(w["width_um"]), core_confluence=geo["share"])]),
               "fig4ac_window.csv", "The window of Figure 5a and 5c (examples.json).",
               [("cohort", "BOMI2 or BOMI1"), ("patient", "patient identifier"), ("core", "core shown"),
                ("width_um", "side of the square window, um, centered on the largest niche, which it holds "
                 "whole"),
                ("core_confluence", "confluence of the core, printed on panel c")])
        for e in ex["suppfig6"]:
            geo = _core_geometry(*cells[e["cohort"]], e["cohort"], e["Case"], e["core"], conf_pc)
            rf, cc = SUPP_CELLS[e["cell"]]
            case = str(e["Case"])
            s_rows.append(dict(cell=f"{rf} fraction, {cc} confluence", cell_code=e["cell"], cohort=e["cohort"], patient=case,
                               core=geo["core"], stromal_cells=geo["n"], core_TrxR1_fraction=geo["rate"],
                               core_confluence=geo["share"], patient_TrxR1_fraction=float(pcase[e["cohort"]][case]),
                               patient_confluence=float(cpc.loc[(e["cohort"], case), f"biggest_{LINK_UM}"])))
        _check_suppfig6(s_rows, D)

    # b: niche areas against the rest, BOMI2, from the per-patient table
    ac = area[(area.cohort == "BOMI2") & area.histology.isin(HISTOLOGIES)]
    for nm in ("T", "Mac", "Fib", "Other"):
        x = ac[[f"inall_{nm}", f"outall_{nm}", f"all_{nm}"]].dropna()
        stat.append(dict(panel="b", group=f"BOMI2 {nm}", stat="% of all stromal cells in niche areas minus rest "
                         "of the stroma (mean over patients, LUSC and LUAD)", value=100 * x[f"all_{nm}"].mean(),
                         lo=100 * x[f"outall_{nm}"].mean(), hi=100 * x[f"inall_{nm}"].mean(),
                         p=stats.wilcoxon_p(x[f"all_{nm}"]), n=len(x)))
    # e: one generator for the four groups, in table order
    boot = np.random.default_rng(SEED_BOOTSTRAP)
    for c in paths.COHORTS:
        for h in HISTOLOGIES:
            v = 100 * D[c][D[c].histology == h].conf.to_numpy()
            m, lo, hi = stats.bootstrap_median_ci(v, boot, reps=BOOT_REPS)
            stat.append(dict(panel="e", group=f"{c} {h}", stat="median confluence %", value=m, lo=lo, hi=hi, n=len(v)))
        p = stats.mann_whitney_p(D[c][D[c].histology == "LUSC"].conf, D[c][D[c].histology == "LUAD"].conf)
        stat.append(dict(panel="e", group=f"{c} LUSC vs LUAD", stat="Mann-Whitney p", p=p))
    # f: confluence against the TrxR1 fraction, LUSC, and the pair with similar confluence
    for c in ("BOMI1", "BOMI2"):
        s = D[c][D[c].histology == "LUSC"]
        rho, p = stats.spearman(s.frac, s.conf)
        stat.append(dict(panel="f", group=f"{c} LUSC", stat="Spearman rho", value=rho, p=p, n=len(s)))
    if ex is not None:
        pair = [str(e["Case"] if isinstance(e, dict) else e) for e in ex["fig4f_squares"]]
        if abs(p2.loc[pair[0], "conf"] - p2.loc[pair[1], "conf"]) > 0.02:
            raise ValueError("the Figure 5f pair no longer has a similar confluence (within 2 points)")
        for c_ in pair:
            stat.append(dict(panel="f", group=f"horizontal pair, patient {c_}", stat="TrxR1+ fraction / confluence",
                             value=p2.loc[c_, "frac"], hi=p2.loc[c_, "conf"], n=1))

    cols = [("panel", "Figure 5 panel letter"), ("group", "what the row describes"), ("stat", "the statistic"),
            ("value", "its value (panel a: window side in um; b: difference in points; e: median in percent; "
             "f pair: TrxR1+ fraction)"),
            ("lo", "lower 95 % limit (b: percent outside)"),
            ("hi", "upper 95 % limit (b: percent inside; f pair: confluence)"),
            ("p", "p value (Wald, Wilcoxon, Mann-Whitney or Spearman as the stat says)"),
            ("n", "patients (1 for an example)"), ("events", "deaths")]
    _write(pd.DataFrame(stat), "fig4_stats.csv",
           "Every statistic of Figure 5. Panels e to h use the LUSC and LUAD patients with a qualifying core, "
           "follow-up and >= 20 stromal cells; hazard ratios per 1 SD, adjusted for stage; bootstrap seed "
           "SEED_BOOTSTRAP, 2,000 resamples. Panel h holds the four drawn series: LUSC in each cohort, and "
           "LUSC and LUAD with both cohorts pooled.", cols)
    pts = pd.concat([D[c].assign(cohort=c)[["cohort", "Case", "histology", "conf", "frac"]] for c in paths.COHORTS])
    _write(pts, "fig4ef_patients.csv", "The patients of Figure 5e and 5f, one row each, in the order the "
           "bootstrap of panel e draws from.",
           [("cohort", "BOMI2 or BOMI1"), ("Case", "patient identifier"), ("histology", "LUSC or LUAD"),
            ("conf", "confluence at 30 um links, patient mean"), ("frac", "TrxR1-positive fraction of stromal cells")])
    ex_cols = [("cohort", "BOMI2 or BOMI1"), ("patient", "patient identifier"), ("core", "core shown"),
               ("stromal_cells", "stromal cells of the core with a TrxR1 call"),
               ("core_TrxR1_fraction", "TrxR1-positive fraction of those"), ("core_confluence", "confluence of the core"),
               ("patient_TrxR1_fraction", "TrxR1-positive fraction of the patient's stromal cells (percase.csv)"),
               ("patient_confluence", "confluence, patient mean over qualifying cores"),
               ("<column>_pct", "the same four values in percent")]

    def pct(rows):
        t = pd.DataFrame(rows)
        for c in ("core_TrxR1_fraction", "core_confluence", "patient_TrxR1_fraction", "patient_confluence"):
            t[f"{c}_pct"] = 100 * t[c]
        return t
    if d_rows:
        _write(pct(d_rows), "fig4d_examples.csv", "The two example cores of Figure 5d (examples.json).",
               [("role", "fragmented or confluent")] + ex_cols)
    if s_rows:
        _write(pct(s_rows), "suppfig6_stats.csv", "The four cores of Supplementary Figure 6 (examples.json).",
               [("cell", "the matrix cell, low or high fraction and confluence"),
                ("cell_code", "the same as in examples.json: LL, LH, HL, HH (fraction, then confluence)")] + ex_cols)


def _check_suppfig6(s_rows, D):
    """Stop the run when a core of Supplementary Figure 6 no longer fits the rule its legend states.

    Low and high are the bottom and top thirds of the cohort's LUSC patients (the patients of Figure
    5f), for the TrxR1-positive fraction and for confluence, except that the fraction is split at
    the median for LH and HL, the two combinations opposite to the correlation. The patient means
    and the core's own values must both fall in the matrix cell of the entry, and the core needs
    >= SUPP_MIN_CELLS stromal cells and a TrxR1-positive fraction >= SUPP_MIN_CORE_FRACTION.
    """
    def level(v, lo, hi):
        return "low" if v < lo else "high" if v > hi else "middle"

    for r in s_rows:
        lusc = D[r["cohort"]][D[r["cohort"]].histology == "LUSC"]
        tag = f"Supplementary Figure 6 {r['cell_code']} ({r['cohort']} patient {r['patient']})"
        if r["patient"] not in set(lusc.Case):
            raise ValueError(f"{tag}: not a LUSC patient of Figure 5f")
        want = SUPP_CELLS[r["cell_code"]]
        lo_f, hi_f = np.quantile(lusc.frac, [0.5, 0.5] if r["cell_code"] in ("LH", "HL") else [1 / 3, 2 / 3])
        lo_c, hi_c = np.quantile(lusc.conf, [1 / 3, 2 / 3])
        for what, f, c in (("patient", r["patient_TrxR1_fraction"], r["patient_confluence"]),
                           ("core", r["core_TrxR1_fraction"], r["core_confluence"])):
            if (level(f, lo_f, hi_f), level(c, lo_c, hi_c)) != want:
                raise ValueError(f"{tag}: the {what} values are no longer {want[0]} fraction, {want[1]} confluence")
        if r["stromal_cells"] < SUPP_MIN_CELLS or r["core_TrxR1_fraction"] < SUPP_MIN_CORE_FRACTION:
            raise ValueError(f"{tag}: the core has fewer than {SUPP_MIN_CELLS} stromal cells or a TrxR1-positive "
                             f"fraction below {SUPP_MIN_CORE_FRACTION}")


def _cancer():
    """The Results paragraph on niche connectivity: confluence against cancer-cell TrxR1, LUSC."""
    ca = _read("ch4_cancer_percase.csv")
    conf = _read("ch4_confluence_percase.csv")[["cohort", "Case", f"biggest_{LINK_UM}"]]
    cs = _merge(ca, conf, "ch4_cancer_percase.csv", "ch4_confluence_percase.csv", on=["cohort", "Case"])
    rows = []
    for cohort in paths.COHORTS:
        u = cs[(cs.cohort == cohort) & (cs.histology == "LUSC")][["tum_ck", f"biggest_{LINK_UM}"]].dropna()
        if len(u) < MIN_CASES:
            raise ValueError(f"{cohort} LUSC: {len(u)} patients with cancer-cell TrxR1 and confluence")
        rho, p = stats.spearman(u.tum_ck, u[f"biggest_{LINK_UM}"])
        rows.append(dict(cohort=cohort, histology="LUSC", metric=CANCER_METRIC, name="confluence", n=len(u),
                         rho=rho, p=p))
    _write(pd.DataFrame(rows), "ch4_cancer_trxr1.csv",
           "The Results paragraph on niche connectivity: Spearman correlation of confluence with cancer-cell "
           "TrxR1 (TrxR1-positive fraction of the cytokeratin-positive Tumor cells, qualifying cores with >= 20 "
           "such cells), LUSC patients, one patient per point.",
           [("cohort", "BOMI2 or BOMI1"), ("histology", "LUSC"),
            ("metric", "connectivity: niche connectivity (biggest_30, patient mean)"),
            ("name", "the measure in words"), ("n", "patients with both values"),
            ("rho", "Spearman rho"), ("p", "its p value")])


def _numbers():
    """The text values that are a maximum, a count or a percentage of the tables above."""
    fs = _read("fig4_stats.csv")
    rows = []

    def add(key, value, what):
        rows.append(dict(key=key, value=float(value), what=what))
    h = fs[(fs.panel == "h") & fs.group.str.startswith("LUSC both ")]
    h = h[h.group.str.extract(r"(\d+) um")[0].astype(int) >= 15]
    add("sweep_LUSC_p_max", h.p.max(), "largest p of the pooled LUSC HR at link distances of 15 to 100 um (the "
        "Results paragraph on survival, printed as all p < 0.05)")
    b = fs[(fs.panel == "b") & (fs.group == "BOMI2 T")]
    if len(b):
        add("fig4b_patients", b.n.iloc[0], "patients of Figure 5b")
    pair = fs[(fs.panel == "f") & fs.group.str.startswith("horizontal pair")]
    if len(pair) == 2:
        add("fig4f_pair_first_fraction_pct", 100 * pair.value.iloc[0], "Figure 5f squares, first patient, TrxR1-positive fraction, percent")
        add("fig4f_pair_second_fraction_pct", 100 * pair.value.iloc[1], "Figure 5f squares, second patient, TrxR1-positive fraction, percent")
        add("fig4f_pair_first_confluence_pct", 100 * pair.hi.iloc[0], "Figure 5f squares, first patient, confluence, percent")
        add("fig4f_pair_second_confluence_pct", 100 * pair.hi.iloc[1], "Figure 5f squares, second patient, confluence, percent")
    a = fs[(fs.panel == "a") & (fs.stat == "window side, um")]
    if len(a):
        add("fig4a_window_um", a.value.iloc[0], "side of the Figure 5a and 5c window, um")
    if _out("fig4ac_window.csv").exists():
        add("fig4c_core_confluence_pct", 100 * _read("fig4ac_window.csv").core_confluence.iloc[0],
            "confluence of the window core of Figure 5c, percent")
    else:
        print("   examples: no fig4ac_window.csv; fig4c_core_confluence_pct skipped")
    _write(pd.DataFrame(rows), "ch4_numbers.csv", "Text values derived from the Chapter 4 tables.",
           [("key", "what the value is, short"), ("value", "the value"), ("what", "the value in words, with its source")])


# ---------------------------------------------------------------- s4 panels

def _grid_mask(points, grid, shape, radius):
    if not len(points):
        return np.zeros(shape, bool)
    return (cKDTree(points).query(grid)[0] <= radius).reshape(shape)


def _draw_core(ax, geo, half, step):
    """A core as in Figure 5d: niche areas filled and outlined, the largest red, cancer cells beige."""
    from matplotlib.colors import to_rgba
    from scipy.ndimage import binary_closing, binary_fill_holes, binary_opening
    xy, trx, big = geo["xy"], geo["trx"], geo["big"]
    gd = np.arange(-half, half + step, step)
    GX, GY = np.meshgrid(gd, gd)
    GP = np.column_stack([GX.ravel(), GY.ravel()])
    m_tu = _grid_mask(geo["tumor"], GP, GX.shape, 12.0)
    m_tu = binary_opening(binary_fill_holes(binary_closing(m_tu, iterations=3)), iterations=2)
    img = np.zeros(m_tu.shape + (4,))
    img[m_tu] = to_rgba(D_TUM)
    ax.imshow(img, origin="lower", extent=(-half, half, -half, half), interpolation="antialiased", zorder=0.5)
    if m_tu.any():
        ax.contour(GX, GY, m_tu.astype(float), levels=[0.5], colors=[D_TUM_EDGE], linewidths=0.4, zorder=0.6)
    ax.scatter(xy[~trx, 0], xy[~trx, 1], s=0.55, color=D_NEG, lw=0, zorder=1)
    m_all = _grid_mask(xy[trx], GP, GX.shape, ZONE_UM)
    m_big = _grid_mask(xy[big], GP, GX.shape, ZONE_UM)
    m_oth = m_all & ~m_big
    for m, col, z in ((m_oth, OTHER, 2), (m_big, D_BIG, 3)):
        img = np.zeros(m.shape + (4,))
        img[m] = to_rgba(col)
        ax.imshow(img, origin="lower", extent=(-half, half, -half, half), interpolation="antialiased", zorder=z)
    if m_oth.any():
        ax.contour(GX, GY, m_oth.astype(float), levels=[0.5], colors=[D_OTH_EDGE], linewidths=0.45, zorder=2.5)
    ax.contour(GX, GY, binary_fill_holes(m_big).astype(float), levels=[0.5], colors=[D_BIG_EDGE],
               linewidths=0.9, zorder=4)
    ax.set_xlim(-half, half)
    ax.set_ylim(half, -half)
    ax.axis("off")
    ax.set_aspect("equal")


def _core_legend(fig, loc, anchor, ncol=4, fontsize=5.9):
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    fig.legend(handles=[Patch(fc=D_BIG, ec=D_BIG_EDGE, lw=0.8, label="largest niche area"),
                        Patch(fc=OTHER, ec=D_OTH_EDGE, lw=0.5, label="other niche areas"),
                        Line2D([], [], marker="o", ls="", color=D_NEG, ms=3.0, label="TrxR1- stromal cells"),
                        Patch(fc=D_TUM, ec=D_TUM_EDGE, lw=0.5, label="cancer cells")],
               loc=loc, bbox_to_anchor=anchor, ncol=ncol, frameon=False, fontsize=fontsize,
               handletextpad=0.4, columnspacing=1.0, handlelength=1.3)


def _window_axes(ax, width):
    ax.set_aspect("equal")
    ax.set_xlim(-width / 2, width / 2)
    ax.set_ylim(width / 2, -width / 2)
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_color(plotting.GRID)
        sp.set_linewidth(0.6)


def _panels_examples(plt, ex, fs):
    """Figure 5a, 5c, 5d and Supplementary Figure 6, drawn from cells.tsv and the example tables."""
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    conf_pc = _read("ch4_confluence_percore.csv")
    need = {e["cohort"] for e in ex["fig4d"] + ex["suppfig6"]} | {ex["fig4a_window"]["cohort"]}
    cells = {c: _stroma_and_tumor(c) for c in need}

    # ---- a and c: the window
    w = ex["fig4a_window"]
    arow = fs[(fs.panel == "a") & (fs.stat == "window side, um")]
    width = float(arow.value.iloc[0])
    geo = _core_geometry(*cells[w["cohort"]], w["cohort"], w["Case"], w["core"], conf_pc)
    c0, in_area = _window(geo, width)
    xy, trx, big = geo["xy"], geo["trx"], geo["big"]
    pad = LINK_UM + 5
    mw = (np.abs(xy - c0) <= width / 2 + pad).all(1)
    xw, trw, inw = xy[mw] - c0, trx[mw], in_area[mw]
    pw = xy[trx & mw] - c0
    links = cKDTree(pw).query_pairs(float(LINK_UM), output_type="ndarray")
    negm = ~trw
    fig, axes = plt.subplots(1, 2, figsize=(2.75, 1.75))
    ax = axes[0]
    for a_, b_ in links:
        ax.plot(pw[[a_, b_], 0], pw[[a_, b_], 1], color=BIG, lw=0.5, zorder=2, solid_capstyle="round")
    ax.scatter(xw[negm, 0], xw[negm, 1], s=3.0, color=NEGC, lw=0, zorder=1)
    ax.scatter(pw[:, 0], pw[:, 1], s=5.5, color=BIG, lw=0, zorder=3)
    ax.plot([-width / 2 + 14, -width / 2 + 114], [width / 2 - 16] * 2, color=plotting.INK, lw=1.2, zorder=6,
            solid_capstyle="butt")
    ax.text(-width / 2 + 64, width / 2 - 22, "100 µm", ha="center", va="bottom", fontsize=5.8, zorder=6,
            bbox=dict(fc="white", ec="none", pad=0.5))
    ax = axes[1]
    for pt in pw:
        ax.add_patch(plt.Circle(pt, ZONE_UM, color=HALO, lw=0, zorder=0))
    ax.scatter(xw[negm & ~inw, 0], xw[negm & ~inw, 1], s=3.0, color=NEG_OUT, lw=0, zorder=1)
    ax.scatter(xw[negm & inw, 0], xw[negm & inw, 1], s=3.0, color=NEG_IN, lw=0, zorder=1)
    ax.scatter(pw[:, 0], pw[:, 1], s=5.5, color=BIG, lw=0, zorder=3)
    for ax, t in zip(axes, ["1  link TrxR1+ stromal cells\n    within 30 µm: niches",
                            "2  niche cells and the stroma\n    within 15 µm: niche area"]):
        _window_axes(ax, width)
        ax.text(0.0, 1.03, t, transform=ax.transAxes, ha="left", va="bottom", fontsize=6.2, linespacing=1.15)
    fig.legend(handles=[Line2D([], [], marker="o", ls="", color=BIG, ms=3.0, label="niche cells"),
                        Line2D([], [], color=BIG, lw=0.8, label="links (within 30 µm)"),
                        Patch(fc=HALO, ec="none", label="niche area"),
                        Line2D([], [], marker="o", ls="", color=NEG_IN, ms=3.0, label="TrxR1- stromal cells in a niche area"),
                        Line2D([], [], marker="o", ls="", color=NEG_OUT, ms=3.0, label="TrxR1- stromal cells outside")],
               loc="upper center", bbox_to_anchor=(0.5, 0.06), ncol=2, frameon=False, fontsize=5.8,
               handletextpad=0.4, columnspacing=0.9, handlelength=1.2, labelspacing=0.35)
    plotting.save(fig, CHAPTER, "fig4a_niches")

    cval = float(_read("fig4ac_window.csv").core_confluence.iloc[0])
    fig, ax = plt.subplots(figsize=(1.6, 2.2))
    for pt in pw:
        ax.add_patch(plt.Circle(pt, ZONE_UM, color=HALO, lw=0, zorder=0))
    ax.scatter(xw[negm & ~inw, 0], xw[negm & ~inw, 1], s=3.0, color=NEG_OUT, lw=0, zorder=1)
    ax.scatter(xw[negm & inw, 0], xw[negm & inw, 1], s=3.0, color=NEG_IN, lw=0, zorder=1)
    ax.scatter(pw[:, 0], pw[:, 1], s=5.5, color=BIG, lw=0, zorder=3)
    gx = np.linspace(-width / 2, width / 2, 721)
    GX, GY = np.meshgrid(gx, gx)
    dbig = cKDTree(xy[big] - c0).query(np.column_stack([GX.ravel(), GY.ravel()]))[0].reshape(GX.shape)
    from scipy import ndimage
    mask = ndimage.binary_fill_holes(dbig <= ZONE_UM)   # a clean outer boundary, no interior rings
    ax.contour(GX, GY, mask.astype(float), levels=[0.5], colors=[BIG], linewidths=1.1, zorder=2.5)
    _window_axes(ax, width)
    ax.set_xlabel(f"niche connectivity of this core = {100 * cval:.0f} %\n"
                  "cells in the largest niche / cells in all niches", fontsize=6.0, color=plotting.INK)
    fig.legend(handles=[Line2D([], [], marker="o", ls="", color=BIG, ms=3.0, label="niche cells"),
                        Patch(fc=HALO, ec=BIG, lw=1.0, label="largest niche area"),
                        Patch(fc=HALO, ec="none", label="other niche areas"),
                        Line2D([], [], marker="o", ls="", color=NEG_IN, ms=3.0, label="TrxR1- stromal cells in a niche area"),
                        Line2D([], [], marker="o", ls="", color=NEG_OUT, ms=3.0, label="TrxR1- stromal cells outside")],
               loc="upper center", bbox_to_anchor=(0.5, 0.02), ncol=2, frameon=False, fontsize=5.6,
               handletextpad=0.4, columnspacing=0.9, labelspacing=0.35)
    plotting.save(fig, CHAPTER, "fig4c_confluence")

    # ---- d: the two example cores
    dx = _read("fig4d_examples.csv")
    geos = [_core_geometry(*cells[e["cohort"]], e["cohort"], e["Case"], e["core"], conf_pc) for e in ex["fig4d"]]
    half = max(np.abs(g["xy"]).max() for g in geos) + 15
    fig, axes = plt.subplots(1, 2, figsize=(3.2, 2.1))
    for ax, g, e in zip(axes, geos, ex["fig4d"]):
        r = dx[dx.role == e["role"]].iloc[0]
        _draw_core(ax, g, half, 1.5)
        ax.set_title(f"{SHOWN_ROLE.get(e['role'], e['role'])}\nconnectivity {r.core_confluence_pct:.0f} %\n"
                     f"TrxR1+ {r.core_TrxR1_fraction_pct:.0f} % of stromal cells", fontsize=5.8,
                     color=plotting.INK2, linespacing=1.2)
        if g["core"] == geo["core"]:                  # the window of a and c, boxed and tagged as the legend says
            ax.add_patch(plt.Rectangle(c0 - width / 2, width, width, fill=False, ec=plotting.INK, lw=0.8, zorder=5))
            ax.text(c0[0] + width / 2 + 14, c0[1] - width / 2 + 6, "a, c", fontsize=6.5, fontweight="bold", ha="left",
                    va="top", zorder=6, bbox=dict(fc="white", ec="none", pad=0.6, alpha=0.9))
    axes[0].plot([-half + 10, -half + 260], [half - 20] * 2, color=plotting.INK, lw=1.1)
    axes[0].text(-half + 135, half - 34, "250 µm", ha="center", va="bottom", fontsize=6)
    _core_legend(fig, "upper center", (0.5, 0.06), ncol=2)
    plotting.save(fig, CHAPTER, "fig4d_cores")

    # ---- Supplementary Figure 6: the four cores
    sx = _read("suppfig6_stats.csv")
    geos = [_core_geometry(*cells[e["cohort"]], e["cohort"], e["Case"], e["core"], conf_pc) for e in ex["suppfig6"]]
    half = max(max(np.abs(g["xy"]).max(), np.abs(g["tumor"]).max() if len(g["tumor"]) else 0) for g in geos) + 15
    fig, axes = plt.subplots(2, 2, figsize=(5.4, 5.9))
    pos = {("low", "low"): (0, 0), ("low", "high"): (0, 1), ("high", "low"): (1, 0), ("high", "high"): (1, 1)}
    for g, e in zip(geos, ex["suppfig6"]):
        rf, cc = SUPP_CELLS[e["cell"]]
        ax = axes[pos[(rf, cc)]]
        r = sx[(sx.cohort == e["cohort"]) & (sx.patient == str(e["Case"]))].iloc[0]
        _draw_core(ax, g, half, 2.0)
        ax.set_title(f"TrxR1+ {r.core_TrxR1_fraction_pct:.0f} % of stromal cells, connectivity "
                     f"{r.core_confluence_pct:.0f} %", fontsize=6.3, color=plotting.INK)
        ax.text(0.5, -0.02, f"{plotting.COHORT_LABEL[e['cohort']]}", transform=ax.transAxes, ha="center", va="top",
                fontsize=5.9, color=plotting.INK2)
    axes[0, 0].plot([-half + 10, -half + 260], [half - 30] * 2, color=plotting.INK, lw=1.1)
    axes[0, 0].text(-half + 135, half - 45, "250 µm", ha="center", va="bottom", fontsize=6)
    for col, t in ((0, "low niche connectivity"), (1, "high niche connectivity")):
        axes[0, col].annotate(t, (0.5, 1.12), xycoords="axes fraction", ha="center", fontsize=7.5, fontweight="bold")
    for rw, t in ((0, "low TrxR1-positive\nfraction"), (1, "high TrxR1-positive\nfraction")):
        axes[rw, 0].annotate(t, (-0.12, 0.5), xycoords="axes fraction", ha="center", va="center", rotation=90,
                             fontsize=7.5, fontweight="bold")
    _core_legend(fig, "lower center", (0.5, 0.0), fontsize=6.3)
    plotting.save(fig, CHAPTER, "suppfig6_matrix")
    return dx


def s4_panels():
    """Every panel of Figure 5 and Supplementary Figure 6, from the tables re-read from disk."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import to_rgb
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    from matplotlib.transforms import blended_transform_factory
    plotting.style()
    fs = _read("fig4_stats.csv")
    SQ, AD = plotting.HC["LUSC"], plotting.HC["LUAD"]
    ex = _examples()
    dx = None
    if ex is not None and _out("fig4d_examples.csv").exists():
        dx = _panels_examples(plt, ex, fs)
    else:
        print("s4: examples missing, panels a, c, d and Supplementary Figure 6 skipped")

    # ---- b: niche areas against the rest of the stroma, BOMI2
    def tint(c, f):
        return tuple(1 - (1 - np.array(to_rgb(c))) * f)
    fig, ax = plt.subplots(figsize=(2.25, 1.14))
    plotting.clean(ax)
    ax.tick_params(axis="y", length=0)
    rows = [("T", "T cells"), ("Mac", "macrophages"), ("Fib", "fibroblasts"), ("Other", "no lineage marker")]
    bh = 0.36
    for i, (nm, lab) in enumerate(rows):
        r = fs[(fs.panel == "b") & (fs.group == f"BOMI2 {nm}")].iloc[0]
        y0 = -i
        ax.barh(y0 + bh / 2 + 0.02, r.hi, height=bh, color=plotting.LC[nm], lw=0)
        ax.barh(y0 - bh / 2 - 0.02, r.lo, height=bh, color=tint(plotting.LC[nm], 0.35), lw=0)
        ax.text(r.hi + 0.8, y0 + bh / 2 + 0.02, f"{r.hi:.0f} %" + ("   niche areas" if i == 0 else ""),
                va="center", fontsize=5.9)
        ax.text(r.lo + 0.8, y0 - bh / 2 - 0.02, f"{r.lo:.0f} %" + ("   rest of the stroma" if i == 0 else ""),
                va="center", fontsize=5.9, color=plotting.INK2)
        ax.text(1.31, y0, f"{r.value:+.1f}", transform=blended_transform_factory(ax.transAxes, ax.transData),
                ha="right", va="center", fontsize=7.2, fontweight="bold")
    ax.text(1.31, 0.62, "difference\n(points)", transform=blended_transform_factory(ax.transAxes, ax.transData),
            ha="right", va="bottom", fontsize=5.9, color=plotting.INK2, linespacing=1.1)
    ax.set_yticks([-i for i in range(len(rows))])
    ax.set_yticklabels([l for _, l in rows], fontsize=6.4)
    ax.set_ylim(-(len(rows) - 1) - 0.55, 0.62)
    ax.set_xlim(0, 50)
    ax.set_xticks([0, 10, 20, 30, 40, 50])
    ax.set_xticklabels(["0", "10", "20", "30", "40", "50 %"])
    ax.set_xlabel(f"% of all stromal cells, {plotting.COHORT_LABEL['BOMI2']}\n"
                  "(mean over patients, LUSC and LUAD together)", fontsize=6.0)
    plotting.save(fig, CHAPTER, "fig4b_bars")

    pts = _read("fig4ef_patients.csv")
    # ---- e: confluence by histology
    jit = np.random.default_rng(SEED_JITTER)
    fig, ax = plt.subplots(figsize=(1.45, 1.45))
    plotting.clean(ax, left=True)
    xp = {("BOMI2", "LUSC"): 0, ("BOMI2", "LUAD"): 1.15, ("BOMI1", "LUSC"): 2.55, ("BOMI1", "LUAD"): 3.7}
    for c in paths.COHORTS:
        for h, col in (("LUSC", SQ), ("LUAD", AD)):
            v = 100 * pts[(pts.cohort == c) & (pts.histology == h)].conf.to_numpy()
            x0 = xp[(c, h)]
            ax.scatter(x0 + jit.uniform(-0.28, 0.28, len(v)), v, s=1.6, marker=plotting.HM[h],
                       color=col, lw=0, zorder=2)
            r = fs[(fs.panel == "e") & (fs.group == f"{c} {h}")].iloc[0]
            ax.plot([x0 - 0.3, x0 + 0.3], [r.value] * 2, color=plotting.INK, lw=1.4, solid_capstyle="butt", zorder=4)
            ax.plot([x0, x0], [r.lo, r.hi], color=plotting.INK, lw=1.0, zorder=4)
        p = fs[(fs.panel == "e") & (fs.group == f"{c} LUSC vs LUAD")].p.iloc[0]
        x0, x1 = xp[(c, "LUSC")], xp[(c, "LUAD")]
        ax.plot([x0, x0, x1, x1], [101, 104, 104, 101], color=plotting.INK2, lw=0.6)
        ax.text((x0 + x1) / 2, 105.5, f"p = {p:.3f}", ha="center", va="bottom", fontsize=5.6, color=plotting.INK2)
    ax.set_xticks(list(xp.values()))
    ax.set_xticklabels(["LUSC", "LUAD"] * 2, fontsize=5.8)
    ax.tick_params(axis="x", length=0)
    tr = blended_transform_factory(ax.transData, ax.transAxes)
    for c, xm in (("BOMI2", 0.575), ("BOMI1", 3.125)):
        ax.text(xm, -0.15, plotting.COHORT_LABEL[c], transform=tr, ha="center", va="top", fontsize=6.3)
    ax.set_xlim(-0.6, 4.3)
    ax.set_ylim(0, 112)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_ylabel("niche connectivity (%)", fontsize=6.3)
    plotting.save(fig, CHAPTER, "fig4e_histology")

    # ---- f: confluence against the TrxR1 fraction, LUSC
    fig, ax = plt.subplots(figsize=(1.35, 1.45))
    plotting.clean(ax, left=True)
    rho = {}
    for c, col in (("BOMI1", LUSC_LIGHT), ("BOMI2", SQ)):
        s = pts[(pts.cohort == c) & (pts.histology == "LUSC")]
        rho[c] = fs[(fs.panel == "f") & (fs.group == f"{c} LUSC")].value.iloc[0]
        ax.scatter(100 * s.frac, 100 * s.conf, s=4.0, color=col, lw=0, zorder=2 if c == "BOMI2" else 1)
    if dx is not None:
        vp = [(r.patient_TrxR1_fraction_pct, r.patient_confluence_pct, r.role) for r in dx.itertuples()]
        ax.plot([vp[0][0], vp[1][0]], [vp[0][1], vp[1][1]], color=plotting.INK2, lw=0.6, ls=(0, (2.5, 2)), zorder=3)
        for x, y_, tag in vp:
            ax.scatter(x, y_, s=24, facecolor="none", edgecolor=plotting.INK, lw=0.8, zorder=4)
            ax.annotate(SHOWN_ROLE.get(tag, tag), (x, y_), xytext=(x + 10, y_ + (12 if tag == "confluent" else 6)), fontsize=5.6,
                        ha="left", va="center", bbox=dict(fc="white", ec="none", pad=0.6), zorder=5,
                        arrowprops=dict(arrowstyle="-", color=plotting.INK, lw=0.5, shrinkB=3))
    pair = fs[(fs.panel == "f") & fs.group.str.startswith("horizontal pair")]
    if len(pair) == 2:
        hx, hy = list(100 * pair.value), list(100 * pair.hi)
        ax.plot(hx, hy, color=plotting.INK2, lw=0.6, ls=(0, (2.5, 2)), zorder=3)
        ax.scatter(hx, hy, s=22, marker="s", facecolor="none", edgecolor=plotting.INK, lw=0.8, zorder=4)
        ax.text(sum(hx) / 2, max(hy) + 4, "similar connectivity", ha="center", va="bottom", fontsize=5.6, zorder=5,
                bbox=dict(fc="white", ec="none", pad=0.6))
    ax.set_xlim(0, 72)
    ax.set_ylim(0, 100)
    ax.set_xticks([0, 20, 40, 60])
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_xlabel("TrxR1-positive fraction of\nstromal cells (%)", fontsize=6.3)
    ax.set_ylabel("niche connectivity (%)", fontsize=6.3)
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color=SQ, ms=3.0,
                              label=f"{plotting.COHORT_LABEL['BOMI2']}, rho = {rho['BOMI2']:.2f}"),
                       Line2D([], [], marker="o", ls="", color=LUSC_LIGHT, ms=3.0,
                              label=f"{plotting.COHORT_LABEL['BOMI1']}, rho = {rho['BOMI1']:.2f}")],
              loc="upper left", bbox_to_anchor=(-0.02, 1.05), frameon=True, facecolor="white", edgecolor="none",
              framealpha=0.9, fontsize=5.6, handletextpad=0.05, borderpad=0.2, labelspacing=0.2)
    plotting.save(fig, CHAPTER, "fig4f_fraction")

    # ---- g: survival forest
    fig, ax = plt.subplots(figsize=(1.25, 1.45))
    plotting.clean(ax)
    ax.tick_params(axis="y", length=0)
    y, ticks, labs = 0.0, [], []
    for gname, col, h in (("squamous", SQ, "LUSC"), ("adenocarcinoma", AD, "LUAD")):
        ax.text(-0.04, y, gname, transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=6.6,
                fontweight="bold")
        y -= 1
        cm = plotting.HM[h]                 # cohort rows take the histology marker; the pooled row stays a diamond
        for lab_, k, mk in ((plotting.COHORT_LABEL["BOMI2"], "BOMI2", cm),
                            (plotting.COHORT_LABEL["BOMI1"], "BOMI1", cm),
                            ("combined", "both", "D")):
            r = fs[(fs.panel == "g") & (fs.group == f"{h} {k}")].iloc[0]
            sig = r.p < 0.05
            ax.plot([r.lo, r.hi], [y, y], color=col, lw=1.2, solid_capstyle="round", zorder=3)
            ax.plot(r.value, y, mk, color=col, mfc=col if sig else "white", mec=col, mew=1.0,
                    ms=4.6 if mk == "D" else 4.2, zorder=4)
            ticks.append(y)
            labs.append(lab_)
            y -= 1
        y -= 0.2
    ax.plot([1, 1], [y + 0.7, 0.6], color=plotting.MUTED, lw=0.7, ls=(0, (2.5, 2.5)), zorder=1)
    ax.set_yticks(ticks)
    ax.set_yticklabels(labs, fontsize=6.2)
    ax.set_ylim(y - 0.4, 0.6)
    ax.set_xscale("log")
    ax.set_xlim(0.45, 1.6)
    ax.minorticks_off()
    ax.set_xticks([0.5, 0.7, 1, 1.4])
    ax.set_xticklabels(["0.5", "0.7", "1", "1.4"])
    ax.set_xlabel("HR per 1 SD of niche connectivity", fontsize=6.3)
    ip = {c: fs[(fs.panel == "g") & (fs.group == f"interaction {c}")].p.iloc[0] for c in ("BOMI2", "BOMI1", "both")}
    ax.text(1.0, -0.42, f"histology interaction p = {ip['BOMI2']:.3f} ({plotting.COHORT_LABEL['BOMI2']}), "
            f"{ip['BOMI1']:.3f} ({plotting.COHORT_LABEL['BOMI1']}),\n"
            f"{ip['both']:.3f} (combined)", transform=ax.transAxes, ha="right", va="top", fontsize=5.6,
            color=plotting.INK2, linespacing=1.2)
    # the key of the assembled figure: the diamond for the combined estimate, an open mark for p >= 0.05
    ax.legend(handles=[Line2D([], [], marker="D", ls="", color=plotting.INK2, ms=3.8, label="combined"),
                       Line2D([], [], marker="o", ls="", color=plotting.INK2, mfc="white", ms=3.6,
                              label="p ≥ 0.05")],
              loc="upper right", bbox_to_anchor=(1.0, -0.62), ncol=2, frameon=False, fontsize=5.8,
              handletextpad=0.2, columnspacing=0.9)
    plotting.save(fig, CHAPTER, "fig4g_forest")

    # ---- h: hazard ratio across link distances
    fig, ax = plt.subplots(figsize=(2.15, 1.45))
    plotting.clean(ax, left=True)
    ax.set_xscale("log")
    ax.set_yscale("log")

    def prof(h, k):
        return [fs[(fs.panel == "h") & (fs.group == f"{h} {k} {r} um")].iloc[0] for r in RADII]
    band = tint(SQ, 0.2)
    rs = prof("LUSC", "both")
    ax.fill_between(RADII, [q.lo for q in rs], [q.hi for q in rs], color=band, lw=0, zorder=1)
    for h, col in (("LUAD", AD), ("LUSC", SQ)):
        rs = prof(h, "both")
        ax.plot(RADII, [q.value for q in rs], color=col, lw=1.2, zorder=3)
        for r, q in zip(RADII, rs):
            ax.plot(r, q.value, "D", color=col, mfc=col if q.p < 0.05 else "white", mec=col, mew=0.9, ms=4.0, zorder=4)
    for c, mk, f in (("BOMI2", "o", 0.93), ("BOMI1", "s", 1.075)):
        for r, q in zip(RADII, prof("LUSC", c)):
            ax.plot(r * f, q.value, mk, color=SQ, mfc=SQ if q.p < 0.05 else "white", mec=SQ, mew=0.7, ms=2.6, zorder=2)
    ax.axhline(1, color=plotting.MUTED, lw=0.7, ls=(0, (2.5, 2.5)), zorder=0)
    ax.axvline(LINK_UM, color=plotting.GRID, lw=0.8, zorder=0)
    ax.set_xlim(8.3, 122)
    ax.set_ylim(0.42, 1.9)
    ax.minorticks_off()
    ax.set_xticks(RADII)
    ax.set_xticklabels([str(r) for r in RADII])
    for lab in ax.get_xticklabels():
        if lab.get_text() == str(LINK_UM):
            lab.set_fontweight("bold")
    ax.set_yticks([0.5, 0.7, 1])
    ax.set_yticklabels(["0.5", "0.7", "1"])
    ax.set_xlabel("link distance between TrxR1-positive stromal cells (µm)", fontsize=6.1)
    ax.set_ylabel("HR per 1 SD of\nniche connectivity", fontsize=6.3)
    ax.legend(handles=[Line2D([], [], color=SQ, marker="D", ms=3.4, lw=1.2, label="LUSC, combined"),
                       Line2D([], [], color=AD, marker="D", ms=3.4, lw=1.2, label="LUAD, combined"),
                       Line2D([], [], color=SQ, marker="o", ls="", ms=2.8,
                              label=f"LUSC, {plotting.COHORT_LABEL['BOMI2']}"),
                       Line2D([], [], color=SQ, marker="s", ls="", ms=2.8,
                              label=f"LUSC, {plotting.COHORT_LABEL['BOMI1']}"),
                       Patch(fc=band, ec="none", label="95 % CI, LUSC")],
              loc="upper left", bbox_to_anchor=(0.03, 1.0), ncol=2, frameon=False, fontsize=5.6, handletextpad=0.3,
              columnspacing=0.8, handlelength=1.5)
    plotting.save(fig, CHAPTER, "fig4h_sweep")
    print("s4: panels written to", paths.output(CHAPTER))


# ---------------------------------------------------------------- s5 niche-area densities

def _ball_mask(tree, pts, r):
    """True for every grid point within r of a point of the tree."""
    d, _ = tree.query(pts, k=1, distance_upper_bound=r)
    return np.isfinite(d)


def _density_core(g):
    """Counts and stromal areas of one core, niche areas against the rest of the stroma.

    g holds every cell of the core (Tumor and Stroma) in file order, with x, y in um, lin, pos
    and stroma. Returns none when the core does not pass RULE_CH4 with a missing call read as
    negative. With niches of any size every TrxR1-positive stromal cell is a niche cell, so the niche areas are the
    stroma within ZONE_UM of any of them and need no linking.
    """
    s = g[g.stroma]
    n, npos = len(s), int(s.pos.sum())
    nneg = n - npos
    min_n, min_pos, min_neg = cellsmod.RULE_CH4
    if n < min_n or npos < min_pos or nneg < min_neg:
        return None
    XY = s[["x", "y"]].to_numpy()
    P = s.pos.to_numpy()
    L = s.lin.to_numpy()
    # niche-area cells: the niche cells and every stromal cell within ZONE_UM of a niche cell
    ip = np.where(P)[0]
    tp = cKDTree(XY[ip])
    d, _ = tp.query(XY, k=1, distance_upper_bound=ZONE_UM)
    inside = np.isfinite(d)
    inside[ip] = True
    # a grid of DENSITY_GRID_UM; the stroma is mask A below
    ts = cKDTree(XY)
    x0, y0 = XY.min(0) - 20
    x1, y1 = XY.max(0) + 20
    gx, gy = np.meshgrid(np.arange(x0, x1, DENSITY_GRID_UM), np.arange(y0, y1, DENSITY_GRID_UM))
    G = np.c_[gx.ravel(), gy.ravel()]
    mask_a = _ball_mask(ts, G, 10.0)                      # within 10 um of a stromal cell
    gin = _ball_mask(tp, G, ZONE_UM)                     # grid points in a niche area
    rec = dict(n=n, npos=npos, n_in=int(inside.sum()), n_out=int((~inside).sum()))
    area = DENSITY_GRID_UM ** 2
    for mname, M in (("A", mask_a),):
        rec[f"area_in_{mname}"] = (M & gin).sum() * area
        rec[f"area_out_{mname}"] = (M & ~gin).sum() * area
    for reg, R in (("in", inside), ("out", ~inside)):
        for nm in DENSITY_LINEAGES:
            for st, S in (("pos", P), ("neg", ~P)):
                rec[f"c_{reg}_{nm}_{st}"] = int((R & (L == nm) & S).sum())
    return rec


def _median_or_empty(v):
    """The median of the values that are not missing, or empty (NaN) without a warning when none
    is left, as for a positive-cell density whose outside value is zero in every patient."""
    v = v.dropna()
    return float(v.median()) if len(v) else np.nan


def s5_niche_density():
    """Niche areas as densities (the density statements of the Results paragraph on niche
    composition and of the Figure 5b legend): cells per 1,000 um2 of stroma inside and outside the
    niche areas, DENSITY_COHORT only.

    Every cell of every core is read; a missing TrxR1 call reads as negative (so a core's count
    of stromal cells includes the few without a call), and a missing marker as absent. Stroma is
    measured on a 5 um grid as the grid points within 10 um of a stromal cell (mask A). Per patient,
    the mean over cores of each core's value; over the LUSC and LUAD patients, the median of the inside and outside values and of their ratio,
    and the Wilcoxon p of their difference, for the DENSITY_FEATURES of Supplementary Table 3b and
    the text. Then the printed cells of the density block of Supplementary Table 3b, mask
    DENSITY_PRINTED_MASK, from the summary re-read from disk (_table4c_density).
    """
    cohort = DENSITY_COHORT
    c = cellsmod.load_cells(cohort)
    c["x"], c["y"] = c.x_um.to_numpy(), c.y_um.to_numpy()
    c["lin"] = c.lineage
    c["pos"] = c.TrxR1.fillna(0).to_numpy() == 1
    c["stroma"] = c.tissue == "Stroma"
    rows = []
    for (case, core), g in c.groupby(["Case", "core"], sort=False):
        rec = _density_core(g)
        if rec is not None:
            rows.append(dict(cohort=cohort, Case=str(case), core=core, **rec))
    print(f"s5 {cohort}: {len(rows)} cores", flush=True)
    _write(pd.DataFrame(rows), "ch4_density_percore.csv",
           f"Per core of {cohort}, niche areas (niche cells and every stromal cell within 15 um of them) "
           "against the rest of the stroma, as cell counts and stromal areas. Cores passing RULE_CH4 with a "
           "missing TrxR1 call read as negative. Written by chapter4.py s5.",
           [("cohort", "BOMI2"), ("Case", "patient identifier, as a string"),
            ("core", "core identifier"), ("n", "stromal cells"), ("npos", "of those TrxR1-positive"),
            ("n_in, n_out", "stromal cells inside and outside the niche areas"),
            ("area_<region>_A", "stromal area in um2 of region in (niche areas) or out (the rest), "
                                "mask A (grid points within 10 um of a stromal cell)"),
            ("c_<region>_<lineage>_<state>", "stromal cells of that lineage (T, Fib, Mac) and TrxR1 state "
                                             "(pos, neg) in that region")])
    pc = _read("ch4_density_percore.csv")
    hist = dict(zip(*[_clinical(cohort)[k] for k in ("Case", "histology")]))
    pc["histo"] = pc.Case.map(hist)
    base = pc[pc.histo.isin(HISTOLOGIES) & (pc.n_in >= MIN_SIDE) & (pc.n_out >= MIN_SIDE)]
    out = []
    for mk in (DENSITY_PRINTED_MASK,):
        for feature, nm, states in DENSITY_FEATURES:
            if nm is None:                     # every stromal cell, whatever its lineage and TrxR1 state
                ca, cb = base.n_in, base.n_out
            else:
                ca = sum(base[f"c_in_{nm}_{s}"] for s in states)
                cb = sum(base[f"c_out_{nm}_{s}"] for s in states)
            va = ca / (base[f"area_in_{mk}"] / 1000.0)
            vb = cb / (base[f"area_out_{mk}"] / 1000.0)
            cd = pd.DataFrame({"Case": base.Case, "histo": base.histo, "va": va, "vb": vb})
            cd = cd.replace([np.inf, -np.inf], np.nan).dropna()
            q = cd.groupby(["Case", "histo"])[["va", "vb"]].mean().reset_index()
            d = q.va - q.vb
            out.append(dict(cohort=cohort, comparison="niche_vs_out", measure=f"dens{mk}", feature=feature,
                            hist="both", n_pat=len(q), median_a=q.va.median(), median_b=q.vb.median(),
                            median_ratio=_median_or_empty((q.va / q.vb).replace([np.inf], np.nan)),
                            p=stats.wilcoxon_p(d) if (d != 0).any() else 1.0))
    _write(pd.DataFrame(out), "ch4_density_summary.csv",
           "Niche areas against the rest of the stroma as densities, the features of the density block of "
           "Supplementary Table 3b and the text, LUSC and LUAD patients of BOMI2 together, stroma as mask A "
           "(grid points within 10 um of a stromal cell). Written by chapter4.py s5 from "
           "ch4_density_percore.csv.",
           [("cohort", "BOMI2"), ("comparison", "niche_vs_out: the niche areas against the rest of the stroma"),
            ("measure", "densA: cells per 1,000 um2 of stroma, mask A"),
            ("feature", "Fib_any (all fibroblasts), T_neg (TrxR1-negative T cells), Fib_neg (TrxR1-negative "
                        "fibroblasts), T_any (all T cells), Mac_any (all macrophages), Mac_neg (TrxR1-negative "
                        "macrophages), all_cells (all stromal cells)"),
            ("hist", "both: LUSC and LUAD together"), ("n_pat", "patients"),
            ("median_a", "median over patients inside the niche areas"),
            ("median_b", "median over patients outside"),
            ("median_ratio", "median of the patient ratios inside / outside"),
            ("p", "Wilcoxon signed-rank p of the patient differences inside minus outside")])
    _table4c_density()


def _table4c_density():
    """The density block of Supplementary Table 3b: for each row of TABLE4C_DENSITY_ROWS, the
    median cells per mm2 of stroma inside and outside the niche areas, the median of the patient ratios and the
    Wilcoxon p, mask DENSITY_PRINTED_MASK, from ch4_density_summary.csv re-read from disk. The block header
    prints the patients once, so the run stops when the rows differ in patients."""
    sm = _read("ch4_density_summary.csv")
    sm = sm[(sm.comparison == "niche_vs_out") & (sm["hist"] == "both")   # sm.hist is the DataFrame method
            & (sm.measure == f"dens{DENSITY_PRINTED_MASK}")]
    rows = []
    for feature, label in TABLE4C_DENSITY_ROWS:
        r = sm[sm.feature == feature]
        if len(r) != 1:
            raise ValueError(f"ch4_density_summary.csv holds {len(r)} rows for {feature}, mask "
                             f"{DENSITY_PRINTED_MASK}; expected one")
        r = r.iloc[0]
        inside, outside = PER_MM2 * r.median_a, PER_MM2 * r.median_b
        rows.append(dict(cohort=r.cohort, feature=feature, row=label, mask=DENSITY_PRINTED_MASK,
                         inside=inside, outside=outside, ratio=r.median_ratio, p=r.p, patients=int(r.n_pat),
                         inside_text=f"{inside:,.0f}", outside_text=f"{outside:,.0f}",
                         ratio_text=f"{r.median_ratio:.2f}", p_text=_fmt_p(r.p)))
    if len({r["patients"] for r in rows}) != 1:
        raise ValueError("the rows of the Supplementary Table 3b density block differ in patients")
    _write(pd.DataFrame(rows), "supptable4c_density.csv",
           "Supplementary Table 3b, density block: the printed cells, in the order of the table. "
           "BOMI2 (Uppsala II), LUSC and LUAD patients together; the block header prints the patients once, as "
           "the other blocks of the table do (for example \"Uppsala II, cells per mm2 of stroma "
           f"({rows[0]['patients']} patients)\"); columns inside the niche areas, outside, ratio inside / outside "
           "and p. Cores passing RULE_CH4 with a missing TrxR1 call read "
           "as negative, >= 20 stromal cells inside and outside; per patient the mean over cores; densities are "
           "medians over patients, the ratio is the median of the patient ratios, p is the Wilcoxon signed-rank "
           "p of the patient differences. Stroma is mask A (grid points within 10 um of a stromal cell). "
           "Written by chapter4.py s5 from ch4_density_summary.csv.",
           [("cohort", "BOMI2"), ("feature", "the feature of ch4_density_summary.csv"),
            ("row", "the row label as the table prints it"), ("mask", "A, the stromal mask printed"),
            ("inside, outside", "median cells per mm2 of stroma inside and outside the niche areas"),
            ("ratio", "median of the patient ratios inside / outside"),
            ("p", "Wilcoxon signed-rank p of the patient differences inside minus outside"),
            ("patients", "LUSC and LUAD patients (the same in every row)"),
            ("inside_text, outside_text", "the printed densities, whole cells per mm2"),
            ("ratio_text", "the printed ratio, 2 decimals"),
            ("p_text", "the printed p, 3 decimals, p < 0.001 below")])


# ---------------------------------------------------------------- s6 niche measures and survival

def _niche_survival_set(cohort):
    """The Figure 5g patients with every exposure and covariate of Supplementary Table 3d.

    _survival_set (LUSC and LUAD patients with follow-up, >= 20 stromal cells and a qualifying core) with the
    two niche summaries of ch4_confluence_percase.csv (share, size); from percase.csv the stromal cell density
    (dens, stromal cells per mm2 of stroma over all primary cores, as Supplementary Table 3c) and the T-cell
    proportion of all stromal cells (prop_T, all primary cores, as Figure 1e); age, sex and WHO performance
    status (redox.clinical.load_covariates); and the TMA slide holding most of the patient's cells
    (redox.staining.patient_slide on cells.tsv, as tumor_vs_stroma.py s1 assigns it).
    """
    d = _survival_set(cohort).rename(columns={"niche_area_share": "share", "niche_mean_size": "size"})
    pc = pd.read_csv(paths.intermediate(cohort, "percase.csv"), dtype={"Case": str}, float_precision="round_trip")
    pc["dens"] = pc[[f"dens_{nm}_S" for nm in _names(cohort)]].sum(axis=1)
    d = _merge(d, pc[["Case", "dens", "prop_T"]], "the Figure 5g patients", f"{cohort} percase.csv", on="Case")
    d = d.merge(clinicalmod.load_covariates(cohort), on="Case", how="left")
    slide = staining.patient_slide(cellsmod.load_cells(cohort, columns=["core", "Case"]), cohort)
    slide = slide.rename("slide").reset_index()
    slide["Case"] = slide.Case.astype(str)
    return d.merge(slide, on="Case", how="left")


def _fit_model(x, exposures, extra):
    """One Cox model of overall survival: every exposure per 1 SD, stage, and the extra covariates.

    Rows missing an exposure, follow-up, stage or an extra covariate are dropped, and follow-up of zero or
    less. "slide" stands for one indicator per slide among the analyzed patients, the first in sorted order
    the reference. The hazard ratio of each exposure is
    survival.hr_per_sd of that exposure with the other exposures (standardized over the same rows) and the
    covariates beside it, so that every term is per 1 SD of the analyzed set and comes from
    survival.hr_per_sd. Returns {exposure: the dictionary of hr_per_sd}.
    """
    plain = [c for c in extra if c != "slide"]
    x = x.dropna(subset=list(exposures) + ["OS_weeks", "event", "stage"] + list(extra))
    x = x[x.OS_weeks > 0].copy()
    covars = ["stage"] + plain
    if "slide" in extra:
        S = pd.get_dummies(x.slide.astype(str), prefix="slide", drop_first=True).astype(float)
        x = pd.concat([x, S], axis=1)
        covars += list(S.columns)
    for e in exposures:
        x[f"_z_{e}"] = (x[e] - x[e].mean()) / x[e].std()
    return {e: survival.hr_per_sd(x, e, covars=[f"_z_{o}" for o in exposures if o != e] + covars)
            for e in exposures}


def _p_cell(p):
    """p as Supplementary Tables 1a and 3d print it: a bound below 0.001, 3 decimals below 0.1, else 2."""
    p = float(p)
    return "< 0.001" if p < 0.001 else (f"{p:.3f}" if p < 0.1 else f"{p:.2f}")


def s6_niche_survival():
    """Supplementary Table 3d: the niche measures and overall survival.

    Each niche measure alone and beside the TrxR1-positive fraction of stromal cells, the joint models of
    connectivity with the fraction and with stromal density or the T-cell proportion, and connectivity fully
    adjusted, per cohort and histology on the Figure 5g patients, the two
    cohorts combined per term by fixed effect. Writes ch4_niche_survival.csv (every model and term) and
    supptable3d.csv (the printed cells).
    """
    D = {c: _niche_survival_set(c) for c in paths.COHORTS}
    rows = []
    for exposures, extra in SURVIVAL_MODELS:
        model, adjusted = "+".join(exposures), "all" if extra else "stage"
        for h in HISTOLOGIES:
            per = {c: _fit_model(D[c][D[c].histology == h], exposures, extra) for c in paths.COHORTS}
            for term in exposures:
                fits = {c: per[c][term] for c in paths.COHORTS}
                fits["both"] = survival.meta_fixed([fits[c] for c in paths.COHORTS])
                for cohort, r in fits.items():
                    rows.append(dict(model=model, adjusted=adjusted, term=term, histology=h, cohort=cohort,
                                     **{k: r[k] for k in ("HR", "lo", "hi", "p", "logHR", "se", "n", "events")}))
    T = pd.DataFrame(rows)
    # Checks: every model of a histology and cohort uses the same patients (no exposure or covariate is
    # missing), and the stage model of connectivity is Figure 5g.
    n = T.groupby(["histology", "cohort"]).n.nunique()
    if (n > 1).any():
        raise ValueError(f"s6: the models of {n[n > 1].index.tolist()} use different patients")
    f4 = _read("fig4_stats.csv")
    for h in HISTOLOGIES:
        for c in ("BOMI2", "BOMI1", "both"):
            ref = f4[(f4.panel == "g") & (f4.group == f"{h} {c}")].iloc[0]
            got = T[(T.model == "conf") & (T.adjusted == "stage") & (T.histology == h) & (T.cohort == c)].iloc[0]
            if abs(got.HR - ref.value) > 1e-9 or int(got.n) != int(ref.n):
                raise ValueError(f"s6: connectivity {h} {c} gives {got.HR:.6f} (n {got.n}), Figure 5g "
                                 f"{ref.value:.6f} (n {ref.n})")
    _write(T, "ch4_niche_survival.csv",
           "Supplementary Table 3d: Cox models of overall survival on the patients of Figure 5g (LUSC and LUAD, "
           "follow-up, >= 20 stromal cells, a qualifying core), per cohort and histology, every exposure per 1 SD "
           "of the analyzed set, adjusted for stage (adjusted = stage) or for stage, age, sex, WHO performance "
           "status and TMA slide (adjusted = all); the two cohorts combined per term by inverse-variance fixed "
           "effect (cohort = both). Written by chapter4.py s6.",
           [("model", "the exposures of the model joined by +: frac (TrxR1-positive fraction of stromal cells, "
             "percase.csv), conf (niche connectivity, biggest_30), share (proportion of stromal cells in niche areas), "
             "size (mean niche size, an isolated TrxR1-positive cell a niche of one), dens (stromal cells per mm2 of "
             "stroma, all primary cores), prop_T (T-cell proportion of all stromal cells, all primary cores)"),
            ("adjusted", "stage, or all (stage, age, sex, WHO performance status and one indicator per TMA slide)"),
            ("term", "the exposure this row's hazard ratio belongs to"),
            ("histology", "LUSC or LUAD"), ("cohort", "BOMI2, BOMI1 or both (fixed-effect combination)"),
            ("HR", "hazard ratio per 1 SD of the term"), ("lo", "lower 95 % Wald confidence limit"),
            ("hi", "upper 95 % Wald confidence limit"), ("p", "Wald p value"),
            ("logHR", "the coefficient on the log hazard scale, which the combination averages"),
            ("se", "its standard error"), ("n", "patients in the model (both: summed over the two cohorts)"),
            ("events", "deaths among them (both: summed)")])
    T = _read("ch4_niche_survival.csv")
    pick = T.set_index(["model", "adjusted", "term", "histology", "cohort"]).sort_index()
    groups = [(h, c) for h in HISTOLOGIES for c in ("BOMI2", "BOMI1", "both")]

    def cell(r):
        return f"{r.HR:.2f} ({r.lo:.2f} to {r.hi:.2f}), p {_p_cell(r.p)}"
    shown = []
    for block, label, model, adjusted, term in TABLE3D_ROWS:
        for h, c in groups:
            r = pick.loc[(model, adjusted, term, h, c)]
            shown.append(dict(block=block, row=label, model=model, adjusted=adjusted, term=term, histology=h,
                              cohort=c, text=cell(r)))
    for h, c in groups:
        r = pick.loc[("frac", "stage", "frac", h, c)]
        shown.append(dict(block="n", row="patients (deaths)", model="", adjusted="", term="", histology=h, cohort=c,
                          text=f"{int(r.n)} ({int(r.events)})", patients=int(r.n), deaths=int(r.events)))
    for label, exposures in TABLE3D_JOINT:
        model = "+".join(exposures)
        for term in exposures:
            for h, c in TABLE3D_JOINT_GROUPS:
                r = pick.loc[(model, "stage", term, h, c)]
                shown.append(dict(block="joint", row=f"{label}: {TERM_LABEL[term]}", model=model, adjusted="stage",
                                  term=term, histology=h, cohort=c, text=cell(r)))
    _write(pd.DataFrame(shown).astype({"patients": "Int64", "deaths": "Int64"}), "supptable3d.csv",
           "Supplementary Table 3d, Niche measures and overall survival: the printed cells, "
           "in the order of the table. Hazard ratio per 1 SD (95 % CI), p; adjusted for stage unless the row says "
           "otherwise. The joint models are printed for LUSC with both cohorts combined only, under the block "
           "header \"Joint models, squamous cell carcinoma, both cohorts combined\". Written by chapter4.py s6 from "
           "ch4_niche_survival.csv.",
           [("block", "single (one measure), full (connectivity fully adjusted), beside (the niche measure beside the TrxR1-positive fraction), n (patients and "
             "deaths) or joint (the joint models)"),
            ("row", "the row label as the supplement prints it; for a joint model, the model and the term "
             "(\"model: term\")"),
            ("model, adjusted, term", "the row of ch4_niche_survival.csv the cell is taken from"),
            ("histology", "LUSC or LUAD"), ("cohort", "BOMI2, BOMI1 or both (the two cohorts combined)"),
            ("text", "the printed cell: HR (95 % CI), p; patients (deaths) on the n row"),
            ("patients", "n row: patients in every model of the column"), ("deaths", "n row: deaths among them")])


SECTIONS = {"s1": s1_percore, "s2": s2_percase, "s3": s3_statistics, "s4": s4_panels,
            "s5": s5_niche_density, "s6": s6_niche_survival}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("RUN:")[0].strip(),
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", nargs="*", choices=sorted(SECTIONS),
                    help="run these sections only; the default runs all of them in order")
    a = ap.parse_args()
    for name in (a.only or sorted(SECTIONS)):
        t0 = time.time()
        print(f"== {name}", flush=True)
        SECTIONS[name]()
        print(f"== {name} done, {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
