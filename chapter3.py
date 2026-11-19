#!/usr/bin/env python3
"""
Chapter 3 analysis: the spatial organization of TrxR1-positive stromal cells.

WHAT: six analysis steps that make Figure 4a-c, Supplementary Figure 5a-c and Supplementary Table 2.
      Each step writes its tables first and draws its panels from those tables, re-read from disk.
      The chapter and file names (chapter3, fig3*, suppfig5*, supptable3) are historical; the items they
      hold are the ones named here.
FROM: work_root/intermediate/<cohort>/cells.tsv (positions in micrometers, lineage, TrxR1 call) and
      the histology column of work_root/intermediate/<cohort>/percase.csv. Uppsala II (BOMI2) positions
      are converted from pixels (0.4964 um per pixel) once, in redox.cells.load_cells; Uppsala I (BOMI1)
      positions are micrometers.
TO:   work_root/output/chapter3/, every table with a sibling .md:
        s1 clustering (Figure 4a)
          fig3a_moran_percore.csv        Moran's I of the TrxR1 label per core, both cohorts, with the
                                         within-lineage shuffled value and its permutation p
          fig3a_clustering.png/pdf       Figure 4a
        s2 neighborhoods (Figure 4b, Supplementary Figure 5c)
          fig3b_neighborhood_percore.csv proportion of each lineage among stromal cells within 30 um,
                                         around TrxR1-positive minus around TrxR1-negative cells of
                                         one center lineage, per core
          fig3b_neighborhood_percase.csv the same, mean over the cores of a patient
          fig3b_fibroblast_center.png/pdf       Figure 4b (BOMI2, fibroblast center)
          suppfig5c_bomi1_unstained_center.png/pdf  Supplementary Figure 5c (BOMI1, CD3- CK- center)
        s3 same-type distances (Supplementary Figure 5b)
          suppfig5b_distance_percore.csv distance to the nearest TrxR1-positive cell of the same type,
                                         TrxR1-positive minus TrxR1-negative cells, per core
          suppfig5b_distance_percase.csv the same, mean over the cores of a patient
          suppfig5b_same_type_distance.png/pdf  Supplementary Figure 5b
        s4 example core (Figure 4c)
          fig3c_candidates.csv           the cores that pass the rule, with the tie-break value
          fig3c_core_cells.csv           every cell of the selected core, shifted to its corner
          fig3c_zoom.csv                 the zoom window and the two circled fibroblasts
          fig3c_map.png/pdf, fig3c_zoom.png/pdf  Figure 4c
        s5 same-label neighbors (Supplementary Figure 5a)
          suppfig5a_same_label_percore.csv  per core, the proportion of same-label neighbors of
                                         TrxR1-positive and of TrxR1-negative cells, observed and after
                                         20 within-lineage shuffles
          suppfig5a_same_label.png/pdf   Supplementary Figure 5a
        s6 cell counts around fibroblasts (Supplementary Table 2)
          supptable3_percase.csv         per patient, the cells within 30 um of TrxR1-positive and of
                                         TrxR1-negative fibroblasts, by type and by the neighbor's state
          supptable3.csv                 Supplementary Table 2
        fig3_stats.csv                   medians, bootstrap 95 % CI and Wilcoxon p for every dot panel,
                                         the Moran summary of Figure 4a and the same-label medians of
                                         Supplementary Figure 5a; each section replaces its own rows
WHY:  one script, one cell table: the statistics are computed and written to tables before any panel is
      drawn, and each random draw has its own generator, so no p value depends on the order of the cores.

Definitions:
  clustering    stromal cells (tissue Stroma) with a TrxR1 call and a position, cores with at least
                50 such cells, at least 5 positive and 5 negative (redox.cells.RULE_CH3). Moran's I
                of the 0/1 label on a row-standardized six-nearest-neighbor graph; 199 shuffles of
                the label among cells of the same lineage in the core; p = (1 + shuffles at or
                above the observed I) / 200. Lineage priority CD3 > CD68 > CK > alpha-SMA in BOMI2,
                CD3 > CK in BOMI1; a cell positive for more than one lineage marker is assigned to
                the first lineage in that order. Each core draws its shuffles from its own generator,
                seeded by SEED and a checksum (CRC-32) of the core identifier, so p does not depend on
                core order.
  neighborhoods stromal cells of cores with at least 100 of them; the neighbors of a cell are the
                other stromal cells within 30 um (inclusive). Center lineages with at least 30 cells
                that have a neighbor, at least 5 of them TrxR1-positive and 5 negative. For each
                neighbor lineage with at least 5 cells in the core, the lineage's proportion of the
                neighborhood, averaged over positive centers minus over negative centers. It is formed
                as the sum of the same difference for the TrxR1-positive and for the TrxR1-negative
                cells of the lineage, in that order, with float32 arithmetic for the adjacency and the
                proportions (kept fixed so that the values are reproducible to the last bit).
                A stromal cell without a TrxR1 call (its intensities are missing) is left out.
  distances     every cell of the core (Tumor and Stroma), a missing TrxR1 call counted negative.
                For a lineage with at least 10 positive and 10 negative cells: mean distance from
                positive cells to the nearest other TrxR1-positive cell of the same lineage minus the
                same from negative cells. BOMI2 lineages T, Fib, Mac by priority; BOMI1 T cells are
                CD3-positive and CK-negative (the BOMI2 T lineage does not exclude CK).
  statistics    one value per patient (mean over cores), LUSC and LUAD only; the median with a 95 %
                percentile bootstrap interval (2000 resamples, redox.stats.bootstrap_median_ci, a fresh
                generator seeded with SEED for every row) and the two-sided Wilcoxon signed-rank p
                against zero.
  example core  see s4_example_core; the zoom proportions count every cell of the core (Tumor and
                Stroma) within 30 um.
  same label    the graph of s1 (six nearest stromal neighbors, cores of RULE_CH3). For the cells of
                each label, the proportion of their neighbors that carry the same label, observed and
                averaged over 20 shuffles of the label among cells of the same lineage in the core;
                their ratio is drawn. Positions are the raw export positions (a nearest-neighbor
                graph does not depend on the unit, and the raw grid keeps ties as they are in the
                export). Each core draws its shuffles from numpy.random.default_rng(CRC-32 of the core
                identifier), with no constant added.
  cell counts   BOMI2 stromal cells, cores with at least 100 of them; a fibroblast counts when it has
                a stromal neighbor within 30 um (inclusive); at least 30 such fibroblasts, 5
                TrxR1-positive and 5 TrxR1-negative. A missing TrxR1 call reads as negative. Per core,
                the mean over TrxR1-positive and over TrxR1-negative fibroblasts; per patient, the mean
                over cores; per histology, medians over patients, the median of the patient
                differences with a bootstrap 95 % CI (seed SEED_TABLE3 + CRC-32 of
                "<histology>_<measure>", one generator per row), the Wilcoxon p and the percent of
                patients with a higher value around TrxR1-positive fibroblasts. The measure keys keep
                the marker names ASMA, CD3 and CD68 for fibroblasts, T cells and macrophages because the
                seeds are computed from them. The blocks are the fibroblast, T-cell and macrophage
                blocks.

Seeds:
SEED = 20260914 (the Moran shuffles of s1, the bootstrap intervals and the dot jitter of the panels).
s5 does not use SEED: its shuffles are seeded per core by the CRC-32 checksum of the core identifier.
SEED_TABLE3 = 20260922 (the bootstrap of Supplementary Table 2, one generator per row; the constant's
name is historical).

RUN: python3 chapter3.py [--only s1 s2 s3 s4 s5 s6]
     s4 reads the tables of s1 and s2 from disk, so run those first.
"""
import argparse
import zlib

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial import cKDTree

from redox import cells as cellsmod
from redox import paths, plotting, stats

CHAPTER = 3
SEED = 20260914

# s1, clustering
K_NEIGHBORS = 6
N_SHUFFLES = 199

# s2, neighborhoods
RADIUS_UM = 30.0
MIN_CORE_STROMA = 100          # stromal cells a core needs before any center is considered
MIN_CENTER = 30                # center cells with at least one neighbor
MIN_CENTER_STATE = 5           # TrxR1-positive and TrxR1-negative centers, each
MIN_SOURCE = 5                 # cells a neighbor lineage needs in the core
# Center lineages: BOMI2 uses the fibroblast center only; BOMI1 has no fibroblast marker, so its
# center is the CD3- CK- group ("Other").
CENTERS = {"BOMI2": ["Fib"], "BOMI1": ["Other"]}

# s3, same-type distances
MIN_SAME_TYPE = 10
SAME_TYPE = {"BOMI2": ["T", "Fib", "Mac"], "BOMI1": ["T"]}

# s5, same-label neighbors
N_SAME_LABEL_SHUFFLES = 20
BAND = (0.4, 0.6)              # TrxR1-positive fraction of the cores of the legend's band, (lo, hi]

# s6, cell counts around fibroblasts (Supplementary Table 2)
SEED_TABLE3 = 20260922
# the printed rows in the order of the table: measure key and row label. Every row has its own
# bootstrap generator.
TABLE3_ROWS = [("A_ASMA", "fibroblasts"), ("S_ASMA_pos", "TrxR1-positive"),
               ("S_ASMA_neg", "TrxR1-negative"), ("A_CD3", "T cells"),
               ("S_CD3_pos", "TrxR1-positive"), ("S_CD3_neg", "TrxR1-negative"),
               ("A_CD68", "macrophages"), ("S_CD68_pos", "TrxR1-positive"),
               ("S_CD68_neg", "TrxR1-negative")]
# the lineages of the printed blocks, by the marker names the measure keys use (the seeds use them)
OLD_NAME = {lin: old for lin, old in {"Fib": "ASMA", "T": "CD3", "Mac": "CD68"}.items()
            if f"A_{old}" in dict(TABLE3_ROWS)}

# s4, example core: the selection rule of Figure 4c
EX_HISTOLOGY = "LUSC"
EX_MORAN_WINDOW = 0.02         # clustering index within this of the BOMI2 median
EX_STROMA_RANGE = (3000, 9000) # stromal cells with a TrxR1 call
EX_RATE_RANGE = (0.2, 0.4)     # TrxR1-positive rate of those
EX_MIN_POS_FIB = 60            # TrxR1-positive fibroblasts in the core
ZOOM_HALF_UM = 90.0            # half width of the zoom window
ZOOM_SIDE_UM = (65.0, None)    # horizontal offset of the negative fibroblast; upper bound below
ZOOM_VERTICAL_UM = 25.0
ZOOM_TOP_POSITIVE = 40         # positive fibroblasts tried, closest to the core mean first
# The zoom proportions count every cell of the core within 30 um (Tumor and Stroma), and any
# fibroblast of the core can be a center.

HISTOLOGIES = ("LUSC", "LUAD")

# How each lineage is labeled in the dot panels; BOMI1 "Other" holds fibroblasts and macrophages.
ROW_LABEL = {"Fib": "Fibroblasts", "T": "T cells", "Mac": "Macrophages", "CK": "Cytokeratin+ cells",
             "Other": "No lineage marker"}
BOMI1_ROW_LABEL = {"T": "T cells", "CK": "Cytokeratin+ cells", "Other": "CD3- CK- cells"}

STATS_FILE = "fig3_stats.csv"
STATS_COLUMNS = ["panel", "cohort", "center", "row", "histology", "patients", "cores", "median",
                 "ci_lo", "ci_hi", "wilcoxon_p", "shuffled_median", "pct_cores_p05",
                 "same_obs_median", "same_shuf_median"]


# ---------------------------------------------------------------- shared helpers

def _core_rng(core, stream):
    """A generator for one core: SEED, a checksum of the identifier and a stream number.

    The checksum (CRC-32) is stable across runs and platforms, unlike Python's own string hash, so
    the shuffles of a core are the same whatever else is run and in whatever order.
    """
    return np.random.default_rng([SEED, zlib.crc32(str(core).encode()), stream])


def _histology(cohort):
    """Case to histology, from percase.csv (the one source of patient-level quantities)."""
    path = paths.intermediate(cohort, "percase.csv")
    if not path.exists():
        raise FileNotFoundError(f"{path} is missing. Run mif/03_percase.py first.")
    t = pd.read_csv(path, usecols=["Case", "histology"], dtype={"Case": str, "histology": str})
    return dict(zip(t.Case, t.histology))


def _attach_histology(frame, cohort, cells_name="cells.tsv"):
    """Add the histology column, and stop when no row finds its patient."""
    h = _histology(cohort)
    frame["histology"] = frame.Case.map(h)
    if len(frame) and frame.histology.notna().sum() == 0:
        raise ValueError(f"no patient of {paths.intermediate(cohort, cells_name)} is in "
                         f"{paths.intermediate(cohort, 'percase.csv')}")
    return frame


def _sidecar(name, what, columns):
    """Write the sibling .md that describes every column of an output table."""
    lines = [f"# {name}.csv", "", what, "", "| column | meaning |", "|---|---|"]
    lines += [f"| {c} | {m} |" for c, m in columns] + [""]
    paths.output(CHAPTER, f"{name}.md").write_text("\n".join(lines))


def _read(name):
    """A chapter table, re-read from disk with identifiers as strings."""
    path = paths.output(CHAPTER, name)
    if not path.exists():
        raise FileNotFoundError(f"{path} is missing; run the section that writes it first.")
    return pd.read_csv(path, dtype={"Case": str, "core": str}, float_precision="round_trip")


def _dot_stats(v):
    """Median, bootstrap 95 % CI and Wilcoxon p."""
    v = np.asarray(v, dtype=float)
    med, lo, hi = stats.bootstrap_median_ci(v, SEED)
    return dict(patients=int(len(v)), median=med, ci_lo=lo, ci_hi=hi,
                wilcoxon_p=stats.wilcoxon_p(v))


def _update_stats(panels, rows):
    """Replace the rows of the named panels in fig3_stats.csv, keeping every other panel."""
    path = paths.output(CHAPTER, STATS_FILE)
    new = pd.DataFrame(rows, columns=STATS_COLUMNS)
    if path.exists():
        old = pd.read_csv(path, float_precision="round_trip")
        old = old[~old.panel.isin(panels)]
        new = pd.concat([old, new], ignore_index=True) if len(old) else new
    order = {"Fig3a": 0, "Fig3b": 1, "SuppFig5a": 2, "SuppFig5a_band": 3, "SuppFig5b": 4,
             "SuppFig5c": 5}
    new = new.sort_values("panel", key=lambda s: s.map(order), kind="stable")
    new.to_csv(path, index=False)
    _sidecar("fig3_stats",
             "Summary statistics of Figure 4 and Supplementary Figure 5, one row per drawn group. "
             "Written by chapter3.py; s1 writes panel Fig3a, s2 Fig3b and SuppFig5c, s3 SuppFig5b, "
             "s5 SuppFig5a and SuppFig5a_band, each replacing only its own rows.",
             [("panel", "Fig3a, Fig3b (Figure 4a, 4b), SuppFig5a, SuppFig5a_band, SuppFig5b or "
                        "SuppFig5c (Supplementary Figure 5a, the band of the 5a legend, 5b, 5c)"),
              ("cohort", "BOMI2 or BOMI1; both for SuppFig5a_band"),
              ("center", "the lineage of the TrxR1-positive and TrxR1-negative centers (dot panels)"),
              ("row", "the lineage the value is about: neighbor lineage (Fig3b, SuppFig5c) or "
                      "same-type lineage (SuppFig5b); Moran for Fig3a; the label, TrxR1+ or TrxR1-, "
                      "for SuppFig5a and SuppFig5a_band"),
              ("histology", "LUSC or LUAD; all for every core of the cohort (Fig3a, SuppFig5a)"),
              ("patients", "patients with a value"),
              ("cores", "cores (Fig3a, SuppFig5a, SuppFig5a_band)"),
              ("median", "median over patients (dot panels), over cores of the observed Moran's I "
                         "(Fig3a) or over cores of the same-label ratio, observed / shuffled "
                         "(SuppFig5a; SuppFig5a_band: cores whose TrxR1-positive fraction is above "
                         "0.4 and at most 0.6, both cohorts, every histology)"),
              ("ci_lo", "lower limit of the bootstrap 95 % CI of the median"),
              ("ci_hi", "upper limit"),
              ("wilcoxon_p", "two-sided Wilcoxon signed-rank p against zero"),
              ("shuffled_median", "Fig3a: median over cores of the mean shuffled Moran's I"),
              ("pct_cores_p05", "Fig3a: percent of cores with permutation p below 0.05"),
              ("same_obs_median", "SuppFig5a: median over cores of the observed proportion of "
                                  "same-label neighbors (a fraction)"),
              ("same_shuf_median", "SuppFig5a: the same after the shuffles (a fraction)")])
    return path


def _stats_table():
    return pd.read_csv(paths.output(CHAPTER, STATS_FILE), float_precision="round_trip")


# ---------------------------------------------------------------- s1

def _moran_core(xy, x, lab, rng):
    """Moran's I of a 0/1 label on a six-neighbor graph, and its within-lineage shuffles.

    Vectorized over the shuffles: each column of Xp is one shuffle, in which the label is permuted
    among the cells of each lineage separately, keeping every lineage's positive count. The shuffles
    of a lineage are drawn as the rank order of uniform numbers.
    Returns the observed I, the mean I over the shuffles and the permutation p.
    """
    n = len(x)
    kk = min(K_NEIGHBORS, n - 1)
    _, nn = cKDTree(xy).query(xy, k=kk + 1)
    nn = nn[:, 1:]
    ri = np.repeat(np.arange(n), kk)
    W = sparse.csr_matrix((np.full(ri.size, 1.0 / kk), (ri, nn.ravel())), shape=(n, n))
    mean0 = x.mean()
    denom0 = float(((x - mean0) ** 2).sum())

    def moran(Xp):
        Xc = Xp - mean0
        return (Xc * (W @ Xp - mean0)).sum(0) / denom0   # row-standardized: W @ const = const

    i_obs = float(moran(x[:, None])[0])
    Xp = np.repeat(x[:, None], N_SHUFFLES, axis=1)
    for L in pd.unique(lab):
        idx = np.where(lab == L)[0]
        if len(idx) > 1:
            Xp[idx, :] = x[idx][np.argsort(rng.random((len(idx), N_SHUFFLES)), axis=0)]
    ib = moran(Xp)
    p = (1 + int((ib >= i_obs).sum())) / (N_SHUFFLES + 1)
    return i_obs, float(ib.mean()), float(p)


def s1_clustering():
    """Figure 4a: TrxR1 labels cluster beyond lineage composition, per core, both cohorts."""
    rows = []
    for cohort in paths.COHORTS:
        d = cellsmod.load_cells(cohort, tissue="Stroma",
                                columns=["core", "Case", "x_um", "y_um", "lineage", "TrxR1"])
        d = d.dropna(subset=["TrxR1", "x_um", "y_um"])
        ok = cellsmod.qualifying_cores(d, *cellsmod.RULE_CH3)
        done = 0
        for core, idx in d.groupby("core", sort=True).indices.items():
            if core not in ok:
                continue
            g = d.iloc[idx]
            x = g.TrxR1.to_numpy(dtype=float)
            i_obs, i_shuf, p = _moran_core(g[["x_um", "y_um"]].to_numpy(dtype=float), x,
                                           g.lineage.to_numpy(), _core_rng(core, 1))
            rows.append(dict(cohort=cohort, core=core, Case=g.Case.iloc[0], n=len(x),
                             npos=int(x.sum()), rate=float(x.mean()), I_obs=i_obs,
                             I_shuffled=i_shuf, p=p))
            done += 1
        print(f"s1 {cohort}: {done} cores")
    table = pd.DataFrame(rows)
    table = pd.concat([_attach_histology(table[table.cohort == c].copy(), c)
                       for c in paths.COHORTS], ignore_index=True)
    table = table[["cohort", "core", "Case", "histology", "n", "npos", "rate", "I_obs",
                   "I_shuffled", "p"]]
    path = paths.output(CHAPTER, "fig3a_moran_percore.csv")
    table.to_csv(path, index=False)
    _sidecar("fig3a_moran_percore",
             "Moran's I of the TrxR1 label among stromal cells, one row per core that passes the "
             "Chapter 3 rule (50 stromal cells with a call, 5 positive, 5 negative). Six-nearest-"
             f"neighbor graph, {N_SHUFFLES} shuffles of the label among cells of the same lineage, "
             "each core seeded by SEED and a checksum of its identifier. Written by chapter3.py s1.",
             [("cohort", "BOMI2 or BOMI1"), ("core", "core identifier"),
              ("Case", "patient identifier, as a string"),
              ("histology", "LUSC, LUAD or other, from percase.csv"),
              ("n", "stromal cells with a TrxR1 call and a position"),
              ("npos", "of those, TrxR1-positive"), ("rate", "npos / n"),
              ("I_obs", "Moran's I of the observed label"),
              ("I_shuffled", "mean Moran's I over the shuffles"),
              ("p", "(1 + shuffles with I at or above I_obs) / (shuffles + 1)")])
    print(f"s1 -> {path}")

    t = pd.read_csv(path, dtype={"Case": str, "core": str}, float_precision="round_trip")
    srows = []
    for cohort in paths.COHORTS:
        c = t[t.cohort == cohort]
        for histology in ("all",) + HISTOLOGIES:
            s = c if histology == "all" else c[c.histology == histology]
            srows.append(dict(panel="Fig3a", cohort=cohort, center="", row="Moran",
                              histology=histology, patients=int(s.Case.nunique()),
                              cores=int(len(s)), median=float(s.I_obs.median()),
                              shuffled_median=float(s.I_shuffled.median()),
                              pct_cores_p05=float(100 * (s.p < 0.05).mean())))
    spath = _update_stats(["Fig3a"], srows)
    print(f"s1 -> {spath}")
    print(pd.DataFrame(srows)[["cohort", "histology", "cores", "median", "shuffled_median",
                                "pct_cores_p05"]].to_string(index=False))
    _panel_3a(pd.read_csv(path, dtype={"Case": str, "core": str}, float_precision="round_trip"),
              _stats_table())


def _panel_3a(table, st):
    """Figure 4a, drawn from fig3a_moran_percore.csv and the Fig3a rows of fig3_stats.csv."""
    plotting.style()
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    rng = np.random.default_rng(SEED)
    fig = plt.figure(figsize=(1.6, 2.6))
    ax = fig.add_axes([0.30, 0.22, 0.66, 0.64])
    plotting.clean(ax, left=True)
    xpos = {"BOMI2": 0.0, "BOMI1": 1.25}
    for cohort, x0 in xpos.items():
        t = table[table.cohort == cohort]
        ax.scatter(x0 - 0.22 + rng.uniform(-0.12, 0.12, len(t)), t.I_shuffled, s=1.6,
                   color=plotting.SHUF, lw=0, zorder=2)
        for h in HISTOLOGIES:
            u = t[t.histology == h]
            ax.scatter(x0 + 0.22 + rng.uniform(-0.12, 0.12, len(u)), u.I_obs, s=1.6,
                       marker=plotting.HM[h], color=plotting.HC[h], alpha=0.55, lw=0, zorder=3)
        r = st[(st.panel == "Fig3a") & (st.cohort == cohort) & (st.histology == "all")].iloc[0]
        for dx, v in ((-0.22, r.shuffled_median), (0.22, r["median"])):
            ax.plot([x0 + dx - 0.17, x0 + dx + 0.17], [v, v], color=plotting.INK, lw=1.3,
                    zorder=4, solid_capstyle="butt")
        ax.text(x0 + 0.22, 0.62, f"{r.pct_cores_p05:.0f} %\nof cores", ha="center", va="bottom",
                fontsize=5.4, color=plotting.INK2)
    ax.set_xticks(list(xpos.values()))
    ax.set_xticklabels([plotting.COHORT_LABEL[c] for c in xpos], fontsize=6.4)
    ax.tick_params(axis="x", length=0)
    ax.set_ylim(-0.1, 0.75)
    ax.set_yticks([0, 0.2, 0.4, 0.6])
    ax.axhline(0, color=plotting.GRID, lw=0.6, zorder=1)
    ax.set_ylabel("clustering index (Moran's I)", fontsize=6.2, labelpad=2)
    fig.text(0.02, 0.94, "a", fontsize=9.5, fontweight="bold")
    fig.text(0.14, 0.94, "TrxR1+ cells cluster", fontsize=7.3, fontweight="bold")
    fig.legend(handles=[Line2D([], [], marker=plotting.HM["LUSC"], ls="", color=plotting.HC["LUSC"], ms=3.4,
                               label="LUSC"),
                        Line2D([], [], marker=plotting.HM["LUAD"], ls="", color=plotting.HC["LUAD"], ms=3.4,
                               label="LUAD"),
                        Line2D([], [], marker="o", ls="", color=plotting.SHUF, ms=3.4,
                               label="labels shuffled")],
               loc="lower left", bbox_to_anchor=(0.0, 0.0), ncol=2, frameon=False, fontsize=5.6,
               handletextpad=0.2, columnspacing=0.6)
    for p in plotting.save(fig, CHAPTER, "fig3a_clustering"):
        print(f"s1 -> {p}")


# ---------------------------------------------------------------- s2

def _neighborhood_core(xy, lin, trx, lineages, centers):
    """Observed neighborhood differences for one core.

    Returns a list of (center, source, n_center, n_center_pos, obs). The arithmetic follows the
    reference arithmetic to the bit: a float32 adjacency, float32 neighbor counts, the per-center
    proportions divided in float32 and averaged with float64 weights +1/k over positive centers and
    -1/(m - k) over negative ones, separately for the TrxR1-positive and the TrxR1-negative cells of
    the source lineage; obs is the sum of the two, added in that order.
    """
    n = len(xy)
    nbrs = cKDTree(xy).query_ball_point(xy, RADIUS_UM)
    indptr = np.zeros(n + 1, dtype=np.int64)
    indices = []
    for i, lst in enumerate(nbrs):
        a = np.asarray(lst)
        a = a[a != i]
        indices.append(a)
        indptr[i + 1] = indptr[i] + len(a)
    indices = np.concatenate(indices) if len(indices) else np.array([], dtype=np.int64)
    A = sparse.csr_matrix((np.ones(len(indices), dtype=np.float32), indices, indptr), shape=(n, n))
    tot = np.asarray(A.sum(1)).ravel()
    keep = tot > 0
    masks = {S: (lin == S) for S in lineages}
    sources = [S for S in lineages if masks[S].sum() >= MIN_SOURCE]

    def counts(t):
        return {S: (A @ (masks[S] & t).astype(np.float32), A @ (masks[S] & ~t).astype(np.float32))
                for S in sources}

    def diff(num, tm, t_center):
        """Mean proportion over positive centers minus over negative centers."""
        k = int(t_center.sum())
        m = len(t_center)
        if k < 1 or k >= m:
            return np.nan
        w = np.full(m, -1.0 / (m - k))
        w[t_center] = 1.0 / k
        dd = tot[tm]
        f = np.zeros(m)
        ok = dd > 0
        f[ok] = num[tm][ok] / dd[ok]
        return float(w @ f)

    v0 = counts(trx)
    out = []
    for T in centers:
        tm = masks[T] & keep
        m = int(tm.sum())
        if m < MIN_CENTER:
            continue
        tt = trx[tm]
        k = int(tt.sum())
        if k < MIN_CENTER_STATE or m - k < MIN_CENTER_STATE:
            continue
        if (tot[tm] > 0).sum() < MIN_CENTER:
            continue
        for S in sources:
            pos_num, neg_num = v0[S]
            # the two halves are kept apart and then added, in this order
            o_pos, o_neg = diff(pos_num, tm, tt), diff(neg_num, tm, tt)
            out.append((T, S, m, k, o_pos + o_neg))
    return out


def s2_neighborhoods():
    """Figure 4b and Supplementary Figure 5c: who sits within 30 um of TrxR1-positive cells.

    A stromal cell without a TrxR1 call is left out of the cores before any neighborhood is formed.
    """
    rows = []
    for cohort in paths.COHORTS:
        d = cellsmod.load_cells(cohort, tissue="Stroma",
                                columns=["core", "Case", "x_um", "y_um", "lineage", "TrxR1"])
        lin = d.lineage.to_numpy(dtype=object).copy()
        missing = d.TrxR1.isna().to_numpy()
        trx = d.TrxR1.fillna(0).to_numpy().astype(bool)
        keep = ~missing
        d, lin, trx = d[keep].reset_index(drop=True), lin[keep], trx[keep]
        lineages = cellsmod.LINEAGES[cohort]
        xy_all = d[["x_um", "y_um"]].to_numpy(dtype=float)
        done = 0
        for core, idx in d.groupby("core", sort=True).indices.items():
            if len(idx) < MIN_CORE_STROMA:
                continue
            xy = xy_all[idx]
            res = _neighborhood_core(xy, lin[idx], trx[idx], lineages, CENTERS[cohort])
            case = d.Case.iloc[idx[0]]
            for T, S, m, k, obs in res:
                rows.append(dict(cohort=cohort, core=core, Case=case, center=T, source=S,
                                 n_center=m, n_center_pos=k, obs=obs))
            done += 1
        print(f"s2 {cohort}: {done} cores with at least {MIN_CORE_STROMA} stromal cells")
    table = pd.DataFrame(rows)
    table = pd.concat([_attach_histology(table[table.cohort == c].copy(), c)
                       for c in paths.COHORTS], ignore_index=True)
    table = table[["cohort", "core", "Case", "histology", "center", "source", "n_center",
                   "n_center_pos", "obs"]]
    path = paths.output(CHAPTER, "fig3b_neighborhood_percore.csv")
    table.to_csv(path, index=False)
    _sidecar("fig3b_neighborhood_percore",
             "Composition of the stromal cells within 30 um of TrxR1-positive against TrxR1-negative "
             "center cells of one lineage on the same core, one row per core, center and neighbor "
             "lineage. Written by chapter3.py s2.",
             [("cohort", "BOMI2 or BOMI1"), ("core", "core identifier"),
              ("Case", "patient identifier, as a string"), ("histology", "LUSC, LUAD or other"),
              ("center", "lineage of the center cells"),
              ("source", "lineage of the neighbors counted"),
              ("n_center", "center cells with at least one stromal neighbor within 30 um"),
              ("n_center_pos", "of those, TrxR1-positive"),
              ("obs", "the source lineage's proportion of the neighborhood, mean around positive "
                      "minus around negative centers (a fraction; computed for the TrxR1-positive "
                      "and the TrxR1-negative cells of the lineage and added, as published)")])
    print(f"s2 -> {path}")

    t = pd.read_csv(path, dtype={"Case": str, "core": str}, float_precision="round_trip")
    t = t[t.histology.isin(HISTOLOGIES)]
    pc = (t.groupby(["cohort", "center", "source", "histology", "Case"])
          .agg(cores=("obs", "size"), diff=("obs", "mean")).reset_index())
    ppath = paths.output(CHAPTER, "fig3b_neighborhood_percase.csv")
    pc.to_csv(ppath, index=False)
    _sidecar("fig3b_neighborhood_percase",
             "One dot of Figure 4b and Supplementary Figure 5c: per patient (LUSC and LUAD), the "
             "mean over cores of the difference in a lineage's proportion of stromal cells within "
             "30 um, around TrxR1-positive minus around TrxR1-negative centers. Written by "
             "chapter3.py s2 from fig3b_neighborhood_percore.csv.",
             [("cohort", "BOMI2 or BOMI1"), ("center", "lineage of the center cells"),
              ("source", "lineage whose proportion is compared"),
              ("histology", "LUSC or LUAD"), ("Case", "patient identifier, as a string"),
              ("cores", "cores averaged"), ("diff", "mean over those cores of obs")])
    print(f"s2 -> {ppath}")

    pc = pd.read_csv(ppath, dtype={"Case": str}, float_precision="round_trip")
    panels = [("Fig3b", "BOMI2", "Fib", ["Fib", "T", "Mac", "CK", "Other"]),
              ("SuppFig5c", "BOMI1", "Other", ["T", "CK", "Other"])]
    srows = []
    for panel, cohort, center, sources in panels:
        for S in sources:
            for h in HISTOLOGIES:
                v = pc[(pc.cohort == cohort) & (pc.center == center) & (pc.source == S)
                       & (pc.histology == h)].sort_values("Case")["diff"].to_numpy()
                if len(v):
                    srows.append(dict(panel=panel, cohort=cohort, center=center, row=S,
                                      histology=h, **_dot_stats(v)))
    spath = _update_stats([p[0] for p in panels], srows)
    print(f"s2 -> {spath}")
    print(pd.DataFrame(srows)[["panel", "row", "histology", "patients", "median", "ci_lo",
                                "ci_hi", "wilcoxon_p"]].to_string(index=False))

    pc = pd.read_csv(ppath, dtype={"Case": str}, float_precision="round_trip")
    st = _stats_table()
    _panel_dots(pc.rename(columns={"diff": "value"}), st, "Fig3b", "BOMI2", "Fib",
                ["Fib", "T", "Mac", "CK", "Other"], ROW_LABEL, (-0.25, 0.20),
                [-0.2, -0.1, 0, 0.1, 0.2], ["-0.2", "-0.1", "0", "+0.1", "+0.2"],
                "stromal cells within 30 µm: proportion around\nTrxR1+ minus around TrxR1- fibroblasts",
                "b", "Around TrxR1-positive fibroblasts", "fig3b_fibroblast_center", (2.3, 2.6))
    _panel_dots(pc.rename(columns={"diff": "value"}), st, "SuppFig5c", "BOMI1", "Other",
                ["T", "CK", "Other"], BOMI1_ROW_LABEL, (-0.15, 0.15),
                [-0.1, 0, 0.1], ["-0.1", "0", "+0.1"],
                "proportion within 30 µm,\naround TrxR1+ minus TrxR1-\nCD3- CK- cells",
                "c", f"Around TrxR1+ CD3- CK- cells, {plotting.COHORT_LABEL['BOMI1']}",
                "suppfig5c_bomi1_unstained_center", (2.1, 2.4))


def _panel_dots(values, st, panel, cohort, center, rows, labels, xlim, xticks, xticklabels,
                xlabel, letter, title, name, size, key="source"):
    """One dot panel: a dot per patient, the median tick and its bootstrap CI from fig3_stats.csv.

    values holds the per-patient values in a column named value; st is fig3_stats.csv re-read.
    Values beyond the axis are drawn at its edge.
    """
    plotting.style()
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    rng = np.random.default_rng(SEED)
    xl, xr = xlim
    fig = plt.figure(figsize=size)
    ax = fig.add_axes([0.36, 0.28, 0.60, 0.58])
    plotting.clean(ax)
    clipped = 0
    for i, r in enumerate(rows):
        for h, off in (("LUSC", 0.18), ("LUAD", -0.18)):
            sel = (values.cohort == cohort) & (values[key] == r) & (values.histology == h)
            if "center" in values.columns:
                sel &= values.center == center
            v = values[sel].sort_values("Case").value.to_numpy(dtype=float)
            s = st[(st.panel == panel) & (st.cohort == cohort) & (st.row == r)
                   & (st.histology == h)]
            if not len(v) or s.empty:
                continue
            s = s.iloc[0]
            y = -i + off
            clipped += int(((v < xl) | (v > xr)).sum())
            ax.scatter(np.clip(v, xl, xr), y + rng.uniform(-0.08, 0.08, len(v)), s=1.8,
                       marker=plotting.HM[h], color=plotting.HC[h], alpha=0.4, lw=0, zorder=2)
            ax.plot([s.ci_lo, s.ci_hi], [y, y], color=plotting.INK, lw=2.0,
                    solid_capstyle="butt", zorder=3)
            ax.plot([s["median"]] * 2, [y - 0.12, y + 0.12], color=plotting.INK, lw=1.2, zorder=4)
    ax.axvline(0, color=plotting.MUTED, lw=0.7, zorder=1)
    ax.set_yticks([-i for i in range(len(rows))])
    ax.set_yticklabels([labels[r] for r in rows], color=plotting.INK, fontsize=6.6)
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(-len(rows) + 0.45, 0.45)
    ax.set_xlim(xl - 0.02 * (xr - xl), xr + 0.02 * (xr - xl))
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels)
    ax.set_xlabel(xlabel, fontsize=6.0, labelpad=2)
    fig.text(0.02, 0.94, letter, fontsize=9.5, fontweight="bold")
    fig.text(0.10, 0.94, title, fontsize=7.3, fontweight="bold")
    fig.legend(handles=[Line2D([], [], marker=plotting.HM["LUSC"], ls="", color=plotting.HC["LUSC"], ms=3.4,
                               label="LUSC"),
                        Line2D([], [], marker=plotting.HM["LUAD"], ls="", color=plotting.HC["LUAD"], ms=3.4,
                               label="LUAD"),
                        Line2D([], [], color=plotting.INK, lw=2.0, label="median, 95 % CI")],
               loc="lower left", bbox_to_anchor=(0.0, 0.0), ncol=3, frameon=False, fontsize=5.6,
               handletextpad=0.3, columnspacing=0.8)
    for p in plotting.save(fig, CHAPTER, name):
        print(f"{panel}: {clipped} patient values beyond the axis, drawn at its edge -> {p}")


# ---------------------------------------------------------------- s3

def _nearest_partner_contrast(xy, trx):
    """Distance to the nearest TrxR1-positive partner: d_pos, d_neg and their difference.

    The nearest partner of a cell is the nearest OTHER TrxR1-positive cell of the lineage: two
    neighbors are asked for, and the second is used when the first is the cell itself.
    """
    dist, idx = cKDTree(xy[trx]).query(xy, k=2)
    pid = np.where(trx)[0]
    self_first = pid[idx[:, 0]] == np.arange(len(xy))
    nn = np.where(self_first, dist[:, 1], dist[:, 0])
    d_pos, d_neg = float(nn[trx].mean()), float(nn[~trx].mean())
    return d_pos, d_neg, float(nn[trx].mean() - nn[~trx].mean())


def s3_same_type_distances():
    """Supplementary Figure 5b: TrxR1-positive cells sit nearer TrxR1-positive cells of their type."""
    rows = []
    for cohort in paths.COHORTS:
        cols = ["core", "Case", "x_um", "y_um", "lineage", "TrxR1", "CD3", "CK"]
        d = cellsmod.load_cells(cohort, columns=cols)
        trx_all = d.TrxR1.fillna(0).to_numpy().astype(bool)   # missing call counted negative
        xy_all = d[["x_um", "y_um"]].to_numpy(dtype=float)
        if cohort == "BOMI1":
            # BOMI1 T cells: CD3-positive and CK-negative, in any compartment.
            member = {"T": ((d.CD3.fillna(0) == 1) & (d.CK.fillna(0) == 0)).to_numpy()}
        else:
            member = {L: (d.lineage == L).to_numpy() for L in SAME_TYPE[cohort]}
        for core, idx in d.groupby("core", sort=True).indices.items():
            case = d.Case.iloc[idx[0]]
            for L in SAME_TYPE[cohort]:
                sel = idx[member[L][idx]]
                if len(sel) < 2 * MIN_SAME_TYPE:
                    continue
                trx = trx_all[sel]
                if trx.sum() < MIN_SAME_TYPE or (~trx).sum() < MIN_SAME_TYPE:
                    continue
                d_pos, d_neg, contrast = _nearest_partner_contrast(xy_all[sel], trx)
                rows.append(dict(cohort=cohort, core=core, Case=case, lineage=L,
                                 n_pos=int(trx.sum()), n_neg=int((~trx).sum()), d_pos=d_pos,
                                 d_neg=d_neg, contrast_um=contrast))
        print(f"s3 {cohort}: done")
    table = pd.DataFrame(rows)
    table = pd.concat([_attach_histology(table[table.cohort == c].copy(), c)
                       for c in paths.COHORTS], ignore_index=True)
    table = table[["cohort", "core", "Case", "histology", "lineage", "n_pos", "n_neg", "d_pos",
                   "d_neg", "contrast_um"]]
    path = paths.output(CHAPTER, "suppfig5b_distance_percore.csv")
    table.to_csv(path, index=False)
    _sidecar("suppfig5b_distance_percore",
             "Distance to the nearest TrxR1-positive cell of the same type, one row per core and "
             "lineage; all cells of the core (Tumor and Stroma), a missing TrxR1 call counted "
             f"negative, at least {MIN_SAME_TYPE} positive and {MIN_SAME_TYPE} negative cells of the "
             "lineage. Positions in micrometers (BOMI2 pixels x 0.4964). Written by chapter3.py s3.",
             [("cohort", "BOMI2 or BOMI1"), ("core", "core identifier"),
              ("Case", "patient identifier, as a string"), ("histology", "LUSC, LUAD or other"),
              ("lineage", "T, Fib or Mac (BOMI1: T, meaning CD3-positive and CK-negative)"),
              ("n_pos", "TrxR1-positive cells of the lineage in the core"),
              ("n_neg", "TrxR1-negative cells of the lineage in the core"),
              ("d_pos", "mean distance (um) from the positive cells to the nearest other "
                        "TrxR1-positive cell of the lineage"),
              ("d_neg", "the same from the negative cells"),
              ("contrast_um", "d_pos - d_neg (um); below zero, the positive cells sit nearer")])
    print(f"s3 -> {path}")

    t = pd.read_csv(path, dtype={"Case": str, "core": str}, float_precision="round_trip")
    pc = (t.groupby(["cohort", "lineage", "histology", "Case"]).contrast_um.mean()
          .rename("to_pos_partner").reset_index())
    cores = t.groupby(["cohort", "lineage", "Case"]).core.nunique().rename("cores").reset_index()
    pc = pc.merge(cores, on=["cohort", "lineage", "Case"], how="left")
    pc = pc[["cohort", "lineage", "histology", "Case", "cores", "to_pos_partner"]]
    ppath = paths.output(CHAPTER, "suppfig5b_distance_percase.csv")
    pc.to_csv(ppath, index=False)
    _sidecar("suppfig5b_distance_percase",
             "Per patient, the mean over cores of contrast_um of suppfig5b_distance_percore.csv; the "
             "LUSC and LUAD rows are the dots of Supplementary Figure 5b. Written by chapter3.py s3.",
             [("cohort", "BOMI2 or BOMI1"), ("lineage", "T, Fib or Mac"),
              ("histology", "LUSC, LUAD or other"), ("Case", "patient identifier, as a string"),
              ("cores", "cores averaged"),
              ("to_pos_partner", "distance to the nearest TrxR1-positive cell of the same type, "
                                 "TrxR1-positive minus TrxR1-negative cells (um)")])
    print(f"s3 -> {ppath}")

    pc = pd.read_csv(ppath, dtype={"Case": str}, float_precision="round_trip")
    order = [("BOMI2", "T"), ("BOMI1", "T"), ("BOMI2", "Fib"), ("BOMI2", "Mac")]
    srows = []
    for cohort, L in order:
        for h in HISTOLOGIES:
            v = pc[(pc.cohort == cohort) & (pc.lineage == L) & (pc.histology == h)]
            v = v.sort_values("Case").to_pos_partner.dropna().to_numpy()
            if len(v):
                srows.append(dict(panel="SuppFig5b", cohort=cohort, center=L, row=L,
                                  histology=h, **_dot_stats(v)))
    spath = _update_stats(["SuppFig5b"], srows)
    print(f"s3 -> {spath}")
    print(pd.DataFrame(srows)[["cohort", "row", "histology", "patients", "median", "ci_lo",
                                "ci_hi", "wilcoxon_p"]].to_string(index=False))
    _panel_suppfig5b(pd.read_csv(ppath, dtype={"Case": str}, float_precision="round_trip"),
                     _stats_table(), order)


def _panel_suppfig5b(pc, st, order):
    """Supplementary Figure 5b, drawn from suppfig5b_distance_percase.csv and fig3_stats.csv."""
    plotting.style()
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    rng = np.random.default_rng(SEED)
    xl, xr = -100, 40
    names = {"T": "T cells", "Fib": "Fibroblasts", "Mac": "Macrophages"}
    fig = plt.figure(figsize=(2.75, 2.4))
    ax = fig.add_axes([0.37, 0.28, 0.55, 0.58])
    plotting.clean(ax)
    clipped = 0
    for i, (cohort, L) in enumerate(order):
        for h, off in (("LUSC", 0.18), ("LUAD", -0.18)):
            v = pc[(pc.cohort == cohort) & (pc.lineage == L) & (pc.histology == h)]
            v = v.sort_values("Case").to_pos_partner.dropna().to_numpy()
            s = st[(st.panel == "SuppFig5b") & (st.cohort == cohort) & (st.row == L)
                   & (st.histology == h)]
            if not len(v) or s.empty:
                continue
            s = s.iloc[0]
            y = -i + off
            clipped += int(((v < xl) | (v > xr)).sum())
            ax.scatter(np.clip(v, xl, xr), y + rng.uniform(-0.08, 0.08, len(v)), s=1.8,
                       marker=plotting.HM[h], color=plotting.HC[h], alpha=0.4, lw=0, zorder=2)
            ax.plot([s.ci_lo, s.ci_hi], [y, y], color=plotting.INK, lw=2.0,
                    solid_capstyle="butt", zorder=3)
            ax.plot([s["median"]] * 2, [y - 0.12, y + 0.12], color=plotting.INK, lw=1.2, zorder=4)
    ax.axvline(0, color=plotting.MUTED, lw=0.7, zorder=1)
    ax.set_yticks([-i for i in range(len(order))])
    ax.set_yticklabels([f"{names[L]}, {plotting.COHORT_LABEL[c]}" for c, L in order],
                       color=plotting.INK, fontsize=6.6)
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(-len(order) + 0.45, 0.45)
    ax.set_xlim(xl - 0.02 * (xr - xl), xr + 0.02 * (xr - xl))
    ax.set_xticks([-90, -60, -30, 0, 30])
    ax.set_xlabel("distance to the nearest TrxR1+ cell of the same\ntype, TrxR1+ minus TrxR1- "
                  "cells (µm)", fontsize=6.0, labelpad=2)
    fig.text(0.02, 0.94, "b", fontsize=9.5, fontweight="bold")
    fig.text(0.10, 0.94, "TrxR1+ cells lie nearer to each other", fontsize=7.3,
             fontweight="bold")
    fig.legend(handles=[Line2D([], [], marker=plotting.HM["LUSC"], ls="", color=plotting.HC["LUSC"], ms=3.4,
                               label="LUSC"),
                        Line2D([], [], marker=plotting.HM["LUAD"], ls="", color=plotting.HC["LUAD"], ms=3.4,
                               label="LUAD"),
                        Line2D([], [], color=plotting.INK, lw=2.0, label="median, 95 % CI")],
               loc="lower left", bbox_to_anchor=(0.0, 0.0), ncol=3, frameon=False, fontsize=5.6,
               handletextpad=0.3, columnspacing=0.8)
    for p in plotting.save(fig, CHAPTER, "suppfig5b_same_type_distance"):
        print(f"SuppFig5b: {clipped} patient values beyond the axis, drawn at its edge -> {p}")


# ---------------------------------------------------------------- s4

def _fib_props(xy, lin):
    """Fibroblast proportion of all other cells within 30 um, per fibroblast (empty when none)."""
    nb = cKDTree(xy).query_ball_point(xy, RADIUS_UM)
    fi = np.where(lin == "Fib")[0]
    pf = np.full(len(fi), np.nan)
    for j, i in enumerate(fi):
        a = np.array([k for k in nb[i] if k != i], dtype=int)
        if len(a):
            pf[j] = np.mean(lin[a] == "Fib")
    return fi, pf


def s4_example_core():
    """Figure 4c: one BOMI2 squamous core and a zoom, both chosen by rule.

    Core: histology LUSC, clustering index within 0.02 of the BOMI2 median over all cores, 3,000 to
    9,000 stromal cells with a TrxR1 call, a TrxR1-positive rate of 0.2 to 0.4 (all from
    fig3a_moran_percore.csv), and at least 60 TrxR1-positive fibroblasts. Among those, the core whose
    mean fibroblast proportion around TrxR1-positive minus around TrxR1-negative fibroblasts is
    closest to the LUSC patient median of Figure 4b (from fig3b_neighborhood_percase.csv); a tie
    goes to the first core in identifier order.
    Zoom: among TrxR1-positive fibroblasts at least 90 um from the core's bounding box edge, with at
    least the median number of cells in a 90 um radius, and a defined proportion, the 40 whose
    fibroblast proportion is closest to the core mean for positive fibroblasts are tried in that
    order; the first that has a TrxR1-negative fibroblast 65 to 116 um to the side and within 25 um
    vertically is taken, with the negative fibroblast whose proportion is closest to the core mean
    for negative fibroblasts. The window is 180 um wide, centered between the two.
    When examples.json has a fig3c entry with positive_raw and negative_raw (export pixel positions),
    those two fibroblasts are the pair instead. The pin exists because the order of equally close
    candidates depends on the numpy build, so the rule alone may pick another pair on another
    machine; the pinned cells are checked to be fibroblasts of the selected core with the stated
    TrxR1 call.
    Every cell of the core (Tumor and Stroma) counts as a neighbor, and every fibroblast as a
    center candidate.
    """
    mor = _read("fig3a_moran_percore.csv")
    mor = mor[mor.cohort == "BOMI2"]
    med_i = float(mor.I_obs.median())
    cand = mor[(mor.histology == EX_HISTOLOGY) & ((mor.I_obs - med_i).abs() < EX_MORAN_WINDOW)
               & mor.n.between(*EX_STROMA_RANGE) & mor.rate.between(*EX_RATE_RANGE)]
    pc = _read("fig3b_neighborhood_percase.csv")
    lusc_med = float(pc[(pc.cohort == "BOMI2") & (pc.center == "Fib") & (pc.source == "Fib")
                        & (pc.histology == EX_HISTOLOGY)]["diff"].median())
    d = cellsmod.load_cells("BOMI2", columns=["core", "Case", "tissue", "x_um", "y_um", "x_raw",
                                              "y_raw", "lineage", "TrxR1"])
    d = d[d.core.isin(set(cand.core))].copy()
    d["TrxR1"] = d.TrxR1.fillna(0).astype(int)       # missing call counted negative
    crow = []
    best = None
    for core, cc in d.groupby("core", sort=True):
        cc = cc.reset_index(drop=True)
        n_fib_pos = int(((cc.lineage == "Fib") & (cc.TrxR1 == 1)).sum())
        row = dict(core=core, Case=cc.Case.iloc[0], n_fib_pos=n_fib_pos, diff=np.nan, selected=0)
        if n_fib_pos >= EX_MIN_POS_FIB:
            fi, pf = _fib_props(cc[["x_um", "y_um"]].to_numpy(dtype=float), cc.lineage.to_numpy())
            tp = cc.TrxR1.to_numpy().astype(bool)[fi]
            row["diff"] = float(np.nanmean(pf[tp]) - np.nanmean(pf[~tp]))
            if best is None or abs(row["diff"] - lusc_med) < abs(best[1] - lusc_med):
                best = (core, row["diff"])
        crow.append(row)
    if best is None:
        raise ValueError("no core passes the Figure 4c rule")
    cands = pd.DataFrame(crow).merge(cand[["core", "I_obs", "n", "rate"]], on="core", how="left")
    cands.loc[cands.core == best[0], "selected"] = 1
    cands["bomi2_median_I"] = med_i
    cands["lusc_patient_median_diff"] = lusc_med
    cpath = paths.output(CHAPTER, "fig3c_candidates.csv")
    cands.to_csv(cpath, index=False)
    _sidecar("fig3c_candidates",
             "The BOMI2 cores that pass the Figure 4c rule (LUSC, clustering index within 0.02 of "
             "the cohort median, 3,000 to 9,000 stromal cells with a call, positive rate 0.2 to "
             "0.4), and the tie-break value of those with at least 60 TrxR1-positive fibroblasts. "
             "Written by chapter3.py s4.",
             [("core", "core identifier"), ("Case", "patient identifier, as a string"),
              ("n_fib_pos", "TrxR1-positive fibroblasts in the core"),
              ("diff", "mean fibroblast proportion of all cells within 30 um around TrxR1-positive "
                       "minus around TrxR1-negative fibroblasts (a fraction); empty below 60 "
                       "positive fibroblasts"),
              ("selected", "1 for the core drawn: diff closest to lusc_patient_median_diff"),
              ("I_obs", "clustering index of the core"),
              ("n", "stromal cells with a TrxR1 call"), ("rate", "their TrxR1-positive rate"),
              ("bomi2_median_I", "median clustering index over all BOMI2 cores"),
              ("lusc_patient_median_diff", "LUSC patient median of the Figure 4b fibroblast row")])
    print(f"s4 -> {cpath}")

    core = best[0]
    c = d[d.core == core].reset_index(drop=True).copy()
    c["x"] = c.x_um - c.x_um.min()
    c["y"] = c.y_um - c.y_um.min()
    xy = c[["x", "y"]].to_numpy(dtype=float)
    fi, pf = _fib_props(xy, c.lineage.to_numpy())
    tp = c.TrxR1.to_numpy().astype(bool)[fi]
    mpos, mneg = np.nanmean(pf[tp]), np.nanmean(pf[~tp])
    half = ZOOM_HALF_UM
    win = np.array([len(v) for v in cKDTree(xy).query_ball_point(xy[fi], half)])
    inside = ((xy[fi, 0] > half) & (xy[fi, 0] < c.x.max() - half)
              & (xy[fi, 1] > half) & (xy[fi, 1] < c.y.max() - half))
    okp = np.where(inside & tp & np.isfinite(pf) & (win >= np.median(win)))[0]
    okp = okp[np.argsort(np.abs(pf[okp] - mpos))[:ZOOM_TOP_POSITIVE]]
    lim = half - RADIUS_UM - 2
    pair = _pinned_pair(core, c, fi, tp)
    for ip in ([] if pair else okp):
        dd = xy[fi] - xy[fi[ip]]
        okn = np.where(~tp & np.isfinite(pf) & (np.abs(dd[:, 0]) >= ZOOM_SIDE_UM[0])
                       & (np.abs(dd[:, 0]) <= 2 * lim) & (np.abs(dd[:, 1]) <= ZOOM_VERTICAL_UM))[0]
        if len(okn):
            pair = (ip, okn[np.argmin(np.abs(pf[okn] - mneg))])
            break
    if pair is None:
        raise ValueError(f"no zoom pair found in the selected core ({core})")
    ip, i_n = pair
    (cx, cy), (nx, ny) = xy[fi[ip]], xy[fi[i_n]]
    wx, wy = (cx + nx) / 2, (cy + ny) / 2
    zoom = pd.DataFrame([dict(role="positive_center", x=cx, y=cy),
                         dict(role="negative_center", x=nx, y=ny),
                         dict(role="window", x=wx, y=wy)])
    zoom.insert(0, "core", core)
    zoom["radius_um"] = RADIUS_UM
    zoom["half_um"] = half
    zpath = paths.output(CHAPTER, "fig3c_zoom.csv")
    zoom.to_csv(zpath, index=False)
    _sidecar("fig3c_zoom",
             "The zoom of Figure 4c: the two circled fibroblasts and the window center, in "
             "micrometers from the corner of the selected core, as in fig3c_core_cells.csv. Written "
             "by chapter3.py s4.",
             [("core", "the selected core"),
              ("role", "positive_center, negative_center, or window (its center)"),
              ("x", "position in um from the left edge of the core"),
              ("y", "position in um from the top edge"),
              ("radius_um", "circle radius (um)"), ("half_um", "half width of the window (um)")])
    print(f"s4 -> {zpath}")
    cells_out = c[["x", "y", "tissue", "lineage", "TrxR1"]].copy()
    cells_out["in_window"] = (cells_out.x.between(wx - half, wx + half)
                              & cells_out.y.between(wy - half, wy + half)).astype(int)
    kpath = paths.output(CHAPTER, "fig3c_core_cells.csv")
    cells_out.to_csv(kpath, index=False)
    _sidecar("fig3c_core_cells",
             "Every cell of the Figure 4c core, shifted to start at 0. Written by chapter3.py s4.",
             [("x", "position in um from the left edge of the core"),
              ("y", "position in um from the top edge"),
              ("tissue", "Tumor or Stroma"), ("lineage", "T, Mac, CK, Fib or Other"),
              ("TrxR1", "1 positive, 0 negative or missing"),
              ("in_window", "1 inside the zoom window")])
    print(f"s4 -> {kpath}")
    ex = _example_note(core)
    print(f"s4: core {core}, clustering index {cand.set_index('core').loc[core, 'I_obs']:.3f} "
          f"(BOMI2 median {med_i:.3f}); fibroblast proportion around TrxR1+ {mpos:.3f}, TrxR1- "
          f"{mneg:.3f} (difference {mpos - mneg:+.3f}, LUSC patient median {lusc_med:+.3f}); "
          f"fibroblast proportion around the circled TrxR1+ {pf[ip]:.2f}, TrxR1- {pf[i_n]:.2f}. "
          f"{ex}")
    _panel_3c(pd.read_csv(kpath), pd.read_csv(zpath, dtype={"core": str}))


def _pinned_pair(core, c, fi, tp):
    """The zoom pair of the fig3c entry (Figure 4c) of examples.json, as indices into fi, or none without one."""
    import json
    try:
        ex = json.loads(paths.examples().read_text()).get("fig3c") or {}
    except (KeyError, OSError, ValueError):
        return None
    if not ex.get("positive_raw") or not ex.get("negative_raw"):
        return None
    if str(ex.get("core", "")) not in core:
        raise ValueError(f"examples.json pins the Figure 4c pair in {ex.get('core')}, the rule picks {core}")
    raw = c[["x_raw", "y_raw"]].to_numpy(dtype=float)[fi]
    out = []
    for key, want in (("positive_raw", True), ("negative_raw", False)):
        hit = np.where((raw[:, 0] == float(ex[key][0])) & (raw[:, 1] == float(ex[key][1])))[0]
        if len(hit) != 1 or bool(tp[hit[0]]) != want:
            raise ValueError(f"examples.json fig3c {key} {ex[key]} is not one fibroblast with the stated call")
        out.append(int(hit[0]))
    return tuple(out)


def _example_note(core):
    """Whether the rule's core agrees with the fig3c entry (Figure 4c) of examples.json, when there is one."""
    import json
    try:
        path = paths.examples()
        ex = json.loads(path.read_text()).get("fig3c")
    except (KeyError, OSError, ValueError):
        return "examples.json not readable; nothing to compare."
    if not ex:
        return "examples.json has no fig3c entry; the rule alone decides."
    want = str(ex.get("core", ""))
    return ("The rule agrees with examples.json." if want and want in core
            else f"WARNING: examples.json names {want}, the rule picks another core.")


def _panel_3c(cells, zoom):
    """Figure 4c, the whole core and the zoom, drawn from fig3c_core_cells.csv and fig3c_zoom.csv."""
    plotting.style()
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Circle, Rectangle
    lc = plotting.LC
    neg_gray, tum_gray = "#c9c8c2", "#b9b8b2"
    w = zoom.set_index("role")
    half, r30 = float(w.half_um.iloc[0]), float(w.radius_um.iloc[0])
    wx, wy = w.loc["window", "x"], w.loc["window", "y"]
    cells["pos"] = cells.TrxR1 == 1

    fig = plt.figure(figsize=(2.3, 2.6))
    ax = fig.add_axes([0.02, 0.10, 0.96, 0.78])
    tum = cells[cells.lineage == "CK"]
    st = cells[cells.lineage != "CK"]
    ax.scatter(tum.x, tum.y, s=0.3, color=lc["CK"], lw=0, zorder=1)
    ax.scatter(st[~st.pos].x, st[~st.pos].y, s=0.3, color=neg_gray, lw=0, zorder=2)
    for L in ("Other", "Mac", "T", "Fib"):
        q = st[st.pos & (st.lineage == L)]
        ax.scatter(q.x, q.y, s=0.8, color=lc[L], lw=0, zorder=3)
    ax.add_patch(Rectangle((wx - half, wy - half), 2 * half, 2 * half, fill=False,
                           ec=plotting.INK, lw=0.9, zorder=4))
    ymax = cells.y.max()
    ax.set_aspect("equal")
    ax.set_xlim(-20, cells.x.max() + 20)
    ax.set_ylim(ymax + 20, -20)
    ax.axis("off")
    ax.plot([0, 250], [ymax + 12] * 2, color=plotting.INK, lw=1.2)
    ax.text(125, ymax + 4, "250 µm", ha="center", va="bottom", fontsize=5.8)
    fig.text(0.02, 0.94, "c", fontsize=9.5, fontweight="bold")
    fig.text(0.10, 0.94, "TrxR1-positive cells in one core", fontsize=7.3, fontweight="bold")
    fig.legend(handles=[Line2D([], [], marker="o", ls="", color=lc[L], ms=3.6,
                               label=plotting.LABEL[L]) for L in ("Fib", "T", "Mac", "Other")]
               + [Line2D([], [], marker="o", ls="", color=neg_gray, ms=3.6, label="TrxR1- stromal"),
                  Line2D([], [], marker="o", ls="", color=lc["CK"], ms=3.6, label="tumor")],
               loc="lower left", bbox_to_anchor=(0.0, 0.0), ncol=3, frameon=False, fontsize=5.4,
               handletextpad=0.1, columnspacing=0.5)
    for p in plotting.save(fig, CHAPTER, "fig3c_map"):
        print(f"s4 -> {p}")

    fig = plt.figure(figsize=(1.6, 1.8))
    az = fig.add_axes([0.04, 0.16, 0.92, 0.80])
    win = cells[cells.in_window == 1]
    for L, g in win.groupby("lineage"):
        col = lc[L] if L != "CK" else tum_gray
        az.scatter(g[~g.pos].x, g[~g.pos].y, s=8, facecolor="white", edgecolor=col, lw=0.6,
                   zorder=2)
        az.scatter(g[g.pos].x, g[g.pos].y, s=8, color=col, lw=0, zorder=3)
    for role, filled in (("positive_center", True), ("negative_center", False)):
        x_, y_ = w.loc[role, "x"], w.loc[role, "y"]
        az.add_patch(Circle((x_, y_), r30, fill=False, ec=plotting.INK, lw=0.8, ls=(0, (3, 2)),
                            zorder=4))
        az.plot(x_, y_, marker="o", ms=4.4, mfc=lc["Fib"] if filled else "white",
                mec=plotting.INK, mew=0.9, zorder=5)
    az.set_xlim(wx - half, wx + half)
    az.set_ylim(wy + half, wy - half)
    az.set_aspect("equal")
    az.set_xticks([])
    az.set_yticks([])
    for sp in az.spines.values():
        sp.set_color(plotting.INK2)
        sp.set_linewidth(0.8)
    az.plot([wx - half + 8, wx - half + 58], [wy - half + 8] * 2, color=plotting.INK, lw=1.1)
    az.text(wx - half + 33, wy - half + 11, "50 µm", ha="center", va="top", fontsize=5.6,
            bbox=dict(boxstyle="square,pad=0.05", fc="white", ec="none", alpha=0.85))
    fig.text(0.04, 0.06, "filled TrxR1+, open TrxR1-; circles 30 µm", fontsize=5.4,
             color=plotting.INK2)
    for p in plotting.save(fig, CHAPTER, "fig3c_zoom"):
        print(f"s4 -> {p}")


# ---------------------------------------------------------------- s5

def _same_label(lab, nb):
    """Proportion of the neighbors of TrxR1-positive cells that are positive, and of TrxR1-negative
    cells that are negative. nb holds the six neighbor indices of every cell."""
    nl = lab[nb]
    return nl[lab == 1].mean(), 1 - nl[lab == 0].mean()


def s5_same_label_neighbors():
    """Supplementary Figure 5a: same-label neighbors of TrxR1-positive and TrxR1-negative cells."""
    rows = []
    for cohort in paths.COHORTS:
        d = cellsmod.load_cells(cohort, tissue="Stroma",
                                columns=["core", "Case", "x_raw", "y_raw", "lineage", "TrxR1"])
        d = d[d.TrxR1.notna()]
        hist = _histology(cohort)
        done = 0
        for core, g in d.groupby("core", sort=True):
            n, npos = len(g), int(g.TrxR1.sum())
            min_n, min_pos, min_neg = cellsmod.RULE_CH3
            if n < min_n or npos < min_pos or n - npos < min_neg:
                continue
            xy = g[["x_raw", "y_raw"]].to_numpy(dtype=float)
            lab = g.TrxR1.to_numpy().astype(int)
            lin = g.lineage.to_numpy()
            _, nb = cKDTree(xy).query(xy, k=K_NEIGHBORS + 1)
            nb = nb[:, 1:]
            pos_same, neg_same = _same_label(lab, nb)
            # Seed: the CRC-32 of the core identifier, no constant added (see the header).
            rng = np.random.default_rng(zlib.crc32(str(core).encode()))
            shuffled = []
            for _ in range(N_SAME_LABEL_SHUFFLES):
                L = lab.copy()
                for li in np.unique(lin):
                    m = lin == li
                    L[m] = rng.permutation(L[m])
                shuffled.append(_same_label(L, nb))
            shuffled = np.array(shuffled)
            case = g.Case.iloc[0]
            rows.append(dict(cohort=cohort, core=core, Case=case, histology=hist.get(case),
                             n=n, frac_pos=npos / n, pos_same=pos_same,
                             pos_same_shuf=shuffled[:, 0].mean(), neg_same=neg_same,
                             neg_same_shuf=shuffled[:, 1].mean()))
            done += 1
        print(f"s5 {cohort}: {done} cores")
    table = pd.DataFrame(rows)
    table["pos_ratio"] = table.pos_same / table.pos_same_shuf
    table["neg_ratio"] = table.neg_same / table.neg_same_shuf
    path = paths.output(CHAPTER, "suppfig5a_same_label_percore.csv")
    table.to_csv(path, index=False)
    _sidecar("suppfig5a_same_label_percore",
             "Same-label neighbors on the six-nearest-neighbor graph of the stromal cells with a "
             "TrxR1 call, one row per core that passes the Chapter 3 rule (as Figure 4a). Shuffles: "
             f"{N_SAME_LABEL_SHUFFLES} per core, the label permuted among cells of the same lineage, "
             "each core seeded by the CRC-32 of its identifier. Written by chapter3.py s5.",
             [("cohort", "BOMI2 or BOMI1"), ("core", "core identifier"),
              ("Case", "patient identifier, as a string"),
              ("histology", "LUSC, LUAD or other, from percase.csv"),
              ("n", "stromal cells with a TrxR1 call"),
              ("frac_pos", "their TrxR1-positive fraction"),
              ("pos_same", "proportion of the neighbors of TrxR1-positive cells that are positive"),
              ("pos_same_shuf", "the same, mean over the shuffles"),
              ("neg_same", "proportion of the neighbors of TrxR1-negative cells that are negative"),
              ("neg_same_shuf", "the same, mean over the shuffles"),
              ("pos_ratio", "pos_same / pos_same_shuf, the value drawn for TrxR1+"),
              ("neg_ratio", "neg_same / neg_same_shuf, the value drawn for TrxR1-")])
    print(f"s5 -> {path}")

    t = pd.read_csv(path, dtype={"Case": str, "core": str}, float_precision="round_trip")
    srows = []
    for cohort in paths.COHORTS:
        z = t[t.cohort == cohort]
        for state, col in (("TrxR1+", "pos"), ("TrxR1-", "neg")):
            r = (z[f"{col}_same"] / z[f"{col}_same_shuf"]).to_numpy()
            srows.append(dict(panel="SuppFig5a", cohort=cohort, center="", row=state,
                              histology="all", cores=int(len(z)), median=float(np.median(r)),
                              same_obs_median=float(np.median(z[f"{col}_same"])),
                              same_shuf_median=float(np.median(z[f"{col}_same_shuf"]))))
    lo, hi = BAND
    band = t[(t.frac_pos > lo) & (t.frac_pos <= hi)]
    for state, col in (("TrxR1+", "pos"), ("TrxR1-", "neg")):
        srows.append(dict(panel="SuppFig5a_band", cohort="both", center="", row=state,
                          histology="all", cores=int(len(band)),
                          median=float((band[f"{col}_same"] / band[f"{col}_same_shuf"]).median())))
    stpath = _update_stats(["SuppFig5a", "SuppFig5a_band"], srows)
    print(f"s5 -> {stpath}")
    print(pd.DataFrame(srows)[["panel", "cohort", "row", "cores", "median", "same_obs_median",
                                "same_shuf_median"]].to_string(index=False))
    _panel_suppfig5a(pd.read_csv(path, dtype={"Case": str, "core": str},
                                 float_precision="round_trip"), _stats_table())


def _panel_suppfig5a(t, st):
    """Supplementary Figure 5a, drawn from suppfig5a_same_label_percore.csv and fig3_stats.csv.

    One dot per LUSC or LUAD core; the bar and the number above it are the median over every core
    of the cohort, read from fig3_stats.csv. Values beyond the axis are drawn at its edge.
    """
    plotting.style()
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    rng = np.random.default_rng(SEED)
    lo, hi = 0.87, 3.6
    fig = plt.figure(figsize=(2.1, 2.5))
    ax = fig.add_axes([0.30, 0.27, 0.66, 0.56])
    plotting.clean(ax, left=True)
    xs = {("BOMI2", "pos"): 0, ("BOMI2", "neg"): 0.9, ("BOMI1", "pos"): 2.1, ("BOMI1", "neg"): 3.0}
    clipped = 0
    for (cohort, col), x0 in xs.items():
        z = t[t.cohort == cohort]
        r = z[f"{col}_ratio"].to_numpy()
        for h in HISTOLOGIES:
            m = z.histology.to_numpy() == h
            clipped += int(((r[m] < lo) | (r[m] > hi)).sum())
            ax.scatter(x0 + rng.uniform(-0.22, 0.22, m.sum()), np.clip(r[m], lo, hi), s=1.4,
                       marker=plotting.HM[h], color=plotting.HC[h], alpha=0.5, lw=0, zorder=3)
        state = "TrxR1+" if col == "pos" else "TrxR1-"
        med = st[(st.panel == "SuppFig5a") & (st.cohort == cohort) & (st.row == state)]
        med = float(med["median"].iloc[0])
        ax.plot([x0 - 0.3, x0 + 0.3], [med, med], color=plotting.INK, lw=1.3, zorder=4,
                solid_capstyle="butt")
        ax.text(x0, 3.85, f"{med:.2f}", ha="center", va="bottom", fontsize=5.3,
                color=plotting.INK2, clip_on=False)
    ax.set_yscale("log")
    ax.set_ylim(0.85, 3.7)
    ax.set_yticks([1, 1.5, 2, 3])
    ax.set_yticklabels(["1", "1.5", "2", "3"])
    ax.minorticks_off()
    ax.axhline(1, color=plotting.GRID, lw=0.8, zorder=1)
    ax.set_xticks(list(xs.values()))
    ax.set_xticklabels(["+", "-", "+", "-"], fontsize=7)
    ax.tick_params(axis="x", length=0, pad=1)
    for cohort, x in (("BOMI2", 0.45), ("BOMI1", 2.55)):
        ax.text(x, 0.745, plotting.COHORT_LABEL[cohort], ha="center", va="top", fontsize=6.2,
                color=plotting.INK, clip_on=False)
    ax.text(-0.75, 0.80, "TrxR1", ha="right", va="center", fontsize=5.6, color=plotting.INK2,
            clip_on=False)
    ax.set_ylabel("same-label neighbors,\nobserved / shuffled", fontsize=6.0, labelpad=2)
    fig.text(0.02, 0.95, "a", fontsize=9.5, fontweight="bold")
    fig.text(0.10, 0.95, "Neighbors with the same TrxR1 label", fontsize=7.3, fontweight="bold")
    fig.legend(handles=[Line2D([], [], marker=plotting.HM["LUSC"], ls="", color=plotting.HC["LUSC"], ms=3.4,
                               label="LUSC"),
                        Line2D([], [], marker=plotting.HM["LUAD"], ls="", color=plotting.HC["LUAD"], ms=3.4,
                               label="LUAD"),
                        Line2D([], [], color=plotting.INK, lw=1.3,
                               label="median, one dot per core")],
               loc="lower left", bbox_to_anchor=(0.0, 0.0), ncol=2, frameon=False, fontsize=5.6,
               handletextpad=0.3, columnspacing=0.8)
    for p in plotting.save(fig, CHAPTER, "suppfig5a_same_label"):
        print(f"SuppFig5a: {clipped} core values beyond the axis, drawn at its edge -> {p}")


# ---------------------------------------------------------------- s6

def _table3_core(xy, lin, trx):
    """Supplementary Table 2 for one core: per measure, the mean around positive and negative
    fibroblasts, or none when the core does not qualify.

    The counts come from a sparse neighbor matrix instead of a loop over fibroblasts; every count is
    an exact integer, so the means do not depend on the summation order.
    """
    n = len(xy)
    if n < MIN_CORE_STROMA:
        return None
    nbrs = cKDTree(xy).query_ball_point(xy, RADIUS_UM)
    rows = np.repeat(np.arange(n), [len(v) for v in nbrs])
    cols = np.concatenate([np.asarray(v, dtype=np.int64) for v in nbrs])
    keep = rows != cols
    A = sparse.csr_matrix((np.ones(int(keep.sum())), (rows[keep], cols[keep])), shape=(n, n))
    deg = np.asarray(A.sum(1)).ravel()
    f = np.where(lin == "Fib")[0]
    f = f[deg[f] > 0]
    tf = trx[f]
    if len(f) < MIN_CENTER or tf.sum() < MIN_CENTER_STATE or (~tf).sum() < MIN_CENTER_STATE:
        return None

    def around(mask):
        return (A @ mask.astype(float))[f]

    rec = {}
    for L, old in OLD_NAME.items():
        is_l = lin == L
        rec[f"A_{old}"] = around(is_l)
        rec[f"S_{old}_pos"] = around(is_l & trx)
        rec[f"S_{old}_neg"] = around(is_l & ~trx)
    out = {}
    for m, v in rec.items():
        a, b = v[tf], v[~tf]
        out[m] = (np.nanmean(a) if np.isfinite(a).any() else np.nan,
                  np.nanmean(b) if np.isfinite(b).any() else np.nan)
    return out


def s6_fibroblast_neighbor_counts():
    """Supplementary Table 2: cells within 30 um of TrxR1-positive and TrxR1-negative fibroblasts."""
    d = cellsmod.load_cells("BOMI2", tissue="Stroma",
                            columns=["core", "Case", "x_um", "y_um", "lineage", "TrxR1"])
    lin_all = d.lineage.to_numpy(dtype=object)
    trx_all = d.TrxR1.fillna(0).to_numpy().astype(bool)      # a missing call reads as negative
    xy_all = d[["x_um", "y_um"]].to_numpy(dtype=float)
    rows = []
    for core, idx in d.groupby("core", sort=True).indices.items():
        res = _table3_core(xy_all[idx], lin_all[idx], trx_all[idx])
        if res is None:
            continue
        case = d.Case.iloc[idx[0]]
        for m, (pos, neg) in res.items():
            rows.append(dict(core=core, Case=case, measure=m, pos=pos, neg=neg))
    r = pd.DataFrame(rows)
    r["histology"] = r.Case.map(_histology("BOMI2"))
    if r.histology.notna().sum() == 0:
        raise ValueError(f"no patient of {paths.intermediate('BOMI2', 'cells.tsv')} is in "
                         f"{paths.intermediate('BOMI2', 'percase.csv')}")
    print(f"s6: {r.core.nunique()} qualifying BOMI2 cores")
    r = r[r.histology.isin(HISTOLOGIES)].dropna(subset=["pos", "neg"])
    r["diff"] = r.pos - r.neg
    pc = r.groupby(["histology", "measure", "Case"])[["pos", "neg", "diff"]].mean().reset_index()
    ppath = paths.output(CHAPTER, "supptable3_percase.csv")
    pc.to_csv(ppath, index=False)
    _sidecar("supptable3_percase",
             "Per patient (BOMI2, LUSC and LUAD), the mean over qualifying cores of the per-core "
             "means around TrxR1-positive and around TrxR1-negative fibroblasts. Written by "
             "chapter3.py s6.",
             [("histology", "LUSC or LUAD"),
              ("measure", "A_<type>: cells of that type within 30 um; S_<type>_pos, S_<type>_neg: "
                          "the same split by the neighbor's TrxR1 state. Types by the published marker "
                          "names: ASMA fibroblasts, CD3 T cells, CD68 macrophages"),
              ("Case", "patient identifier, as a string"),
              ("pos", "mean around TrxR1-positive fibroblasts (cells)"),
              ("neg", "mean around TrxR1-negative fibroblasts (cells)"),
              ("diff", "mean over cores of pos - neg (cells)")])
    print(f"s6 -> {ppath}")

    pc = pd.read_csv(ppath, dtype={"Case": str}, float_precision="round_trip")
    label = dict(TABLE3_ROWS)
    out = []
    for (h, m), g in pc.groupby(["histology", "measure"]):
        x = g["diff"].to_numpy()
        med, lo, hi = stats.bootstrap_median_ci(x, SEED_TABLE3 + zlib.crc32(f"{h}_{m}".encode()))
        out.append(dict(histology=h, measure=m, row_label=label[m], patients=int(len(x)),
                        around_pos=float(np.median(g.pos)), around_neg=float(np.median(g.neg)),
                        diff_median=med, ci_lo=lo, ci_hi=hi, wilcoxon_p=stats.wilcoxon_p(x),
                        pct_patients_higher=100 * float(np.mean(x > 0))))
    path = paths.output(CHAPTER, "supptable3.csv")
    pd.DataFrame(out).to_csv(path, index=False)
    _sidecar("supptable3",
             "Supplementary Table 2: cells within 30 um of "
             "TrxR1-positive against TrxR1-negative fibroblasts of the same core, BOMI2, medians over "
             "patients per histology (A_* in total, S_* split by the neighbor's TrxR1 state); the "
             "fibroblast, T-cell and macrophage blocks. Written by chapter3.py s6.",
             [("histology", "LUSC or LUAD"), ("measure", "as in supptable3_percase.md"),
              ("row_label", "the row label of the printed table"),
              ("patients", "patients with a value"),
              ("around_pos", "median over patients around TrxR1-positive fibroblasts"),
              ("around_neg", "median over patients around TrxR1-negative fibroblasts"),
              ("diff_median", "median of the patient differences, positive minus negative"),
              ("ci_lo", "lower limit of its bootstrap 95 % CI (2,000 resamples, seed SEED_TABLE3 "
                        "+ CRC-32 of histology_measure)"),
              ("ci_hi", "upper limit"),
              ("wilcoxon_p", "two-sided Wilcoxon signed-rank p of the differences"),
              ("pct_patients_higher", "percent of patients with a higher value around "
                                      "TrxR1-positive fibroblasts")])
    print(f"s6 -> {path}")
    print(pd.DataFrame(out)[["histology", "measure", "around_pos", "around_neg", "diff_median",
                             "ci_lo", "ci_hi", "wilcoxon_p"]].round(3).to_string(index=False))


# ---------------------------------------------------------------- run

SECTIONS = {"s1": s1_clustering,
            "s2": s2_neighborhoods,
            "s3": s3_same_type_distances,
            "s4": s4_example_core,
            "s5": s5_same_label_neighbors,
            "s6": s6_fibroblast_neighbor_counts}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("RUN:")[0].strip(),
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", nargs="*", choices=sorted(SECTIONS),
                    help="run these sections only; the default runs all of them in order")
    a = ap.parse_args()
    for name in sorted(a.only) if a.only else sorted(SECTIONS):
        SECTIONS[name]()


if __name__ == "__main__":
    main()
