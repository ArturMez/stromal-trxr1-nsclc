#!/usr/bin/env python3
"""
Chapter 1: the TrxR1-positive fraction of stromal cells, the stromal composition and overall survival.

WHAT: the seven analysis steps behind Figure 1b (hazard ratio of the stromal TrxR1-positive fraction),
      1c (hazard ratio per cell-type proportion) and 1e (correlations with composition) and
      Supplementary Tables 1b and 1c, each writing its table first and drawing its panel from that
      table. Figure 1d (the TrxR1-positive fraction of each cell type) and Supplementary Table 1a
      (cancer cells against stromal cells) are written by tumor_vs_stroma.py; Figure 1a is a
      schematic, whose patient sums s5 writes.
FROM: work_root/intermediate/<cohort>/percase.csv and clinical.csv, plus the optional tables of
      the published fibroblast data of Pellinen and colleagues, in the folder named by the
      "pellinen" entry of config.json; s6 also reads cells.tsv and the TrxR1 channel of the raw
      cell export, through redox.staining; s7 reads cells.tsv and areas.csv.
TO:   work_root/intermediate/<cohort>/ch1_percase.csv and ch1_tcell_density.csv (s7) (per patient,
      they stay with the other patient-level tables) and work_root/output/chapter1/:
        fig1b_models.csv     hazard ratio per 1 SD of the stromal TrxR1 fraction, per cohort and
                             histology, stage-adjusted, each histology pooled over the two cohorts,
                             and the histology interaction per cohort and pooled (fixed effect of
                             the two cohorts' interaction terms)
        fig1b_forest.png/pdf Figure 1b, the forest plot, drawn from fig1b_models.csv
        supptable1b.csv      the same model for each lineage proportion of stromal cells
        fig1c_tiles.png/pdf  Figure 1c, the hazard ratio tiles, drawn from supptable1b.csv
        supptable1c.csv      Spearman correlation of the TrxR1 fraction with each proportion,
                             measured and reconstructed
        fig1e_pies.png/pdf   Figure 1e, the correlation pies, drawn from supptable1c.csv
        ch1_numbers.csv      the composition means of the paragraph "Stromal composition and
                             survival" and the patient sums of Figure 1a
        supptable1c_staining.csv
                             the correlations of Figure 1e again, within TMA slide and with the
                             background staining of stromal cells removed (the Supplementary
                             Table 1c caption)
        supptable1c_tcell_density.csv
                             the TrxR1 fraction against T cells per mm2 of stroma and of tumor
                             compartment and against their tumor-to-stroma ratio, three rows of
                             Supplementary Table 1c
WHY:  one script and one Cox solver (redox.survival) produce every number of Figure 1b, 1c and 1e,
      with the tables written before the panels so that each panel is drawn from a table.

Definitions:
  exposure     the TrxR1-positive fraction of a patient's stromal cells, no lineage gating, for
               patients with at least 20 stromal cells (redox.cells.stromal_trxr1_fraction)
  proportions  each lineage as a proportion of all stromal cells of the patient (cells of the
               lineage over all stromal cells, so the denominator is all stromal cells)
  models       Cox with Efron ties, hazard ratio per 1 standard deviation of the exposure inside the
               analyzed set, adjusted for stage as one linear term from I to IV
  interaction  exposure standardized once over squamous and adenocarcinoma together, terms z, LUAD,
               z x LUAD and stage; the Wald p on the product term, cohorts pooled by fixed effect
  reconstructed fibroblast proportion
               p_SMA x f_pool, where p_SMA is the alpha-SMA-positive share of the Pellinen
               four-marker pool and f_pool the share of stromal cells carrying neither the T-cell
               nor the cancer-cell lineage marker. It is the only fibroblast measure BOMI1 can
               carry, because that panel has no alpha-SMA channel.
  staining check
               the rank correlation of Figure 1e after TMA slide indicators, and in addition the
               rank of the stromal background index (the mean TrxR1 intensity of the patient's
               TrxR1-negative stromal cells, at least 10), are regressed out of both variables
               (redox.staining.partial_rank_correlation)
  T-cell density
               (s7) T cells (the T-cell lineage, CD3 first in the priority order) of the Stroma or
               Tumor category over the segmented area of that category, both summed over the
               patient's cores; the tumor-to-stroma ratio is the tumor-compartment density over the
               stromal density

Shown labels name the cohorts Uppsala I (BOMI1) and Uppsala II (BOMI2); code, file names and table
columns keep BOMI1 and BOMI2.

No random numbers are used anywhere in this chapter, so it declares no seed.

RUN: python3 chapter1.py [--only s1 s2 ...]
     s5 reads fig1b_models.csv, written by s2; s6 and s7 read supptable1c.csv, written by s4.
"""
import argparse

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from redox import cells as cellsmod
from redox import clinical as clinicalmod
from redox import paths, plotting, staining, stats, survival

CHAPTER = 1

# The three lineages the manuscript reports, in the order of the figure and of both tables.
CELLS = ["T", "Fib", "Mac"]

# The two identifier offsets tried when the table of Pellinen and colleagues is joined to a cohort.
# BOMI2 shares its numbering; the BOMI1 tables were numbered from a different base. The offset that
# matches more patients is the one used, so that no identifier has to be written into the code.
ID_OFFSETS = (0, 1000)

# Panel colors that belong to this figure alone. The histology and lineage colors come from
# redox.plotting, so that they cannot drift from the other chapters.
SURFACE = "#fcfcfb"      # the surface the tints are blended from
RED = "#e34948"          # a hazard ratio above one, in the tiles
NEUTRAL = "#f0efec"      # a hazard ratio of one
DISC = "#efeee9"         # the empty part of a correlation pie


# ---------------------------------------------------------------- shared loading

def _proportions(cohort):
    """Per patient: the stromal TrxR1 fraction and each lineage as a proportion of stromal cells.

    Read once from intermediate/{cohort}/percase.csv, written by mif/03_percase.py, instead of
    recomputing the exposure and the proportions from cells.tsv a second time. n_stroma counts
    every stromal cell, so that the proportions and the exposure share one denominator. Patients
    with fewer than 20 stromal cells are dropped.
    """
    path = paths.intermediate(cohort, "percase.csv")
    if not path.exists():
        raise FileNotFoundError(
            f"{path} is missing. Run mif/01_clinical.py, mif/02_cells.R and mif/03_percase.py first.")
    keep = (["Case", "n_stroma", "n_stroma_trxr1", "stromal_trxr1"]
            + [f"prop_{name}" for name in cellsmod.LINEAGES[cohort]])
    t = pd.read_csv(path, usecols=keep, dtype={"Case": str}, float_precision="round_trip")
    # The pool the reconstruction fills: stromal cells carrying neither lineage marker that both
    # panels share, that is neither the T-cell nor the cancer-cell one.
    t["f_pool"] = 1.0 - t.prop_T - t.prop_CK
    t = t[t.n_stroma >= cellsmod.MIN_STROMAL_CELLS].reset_index(drop=True)
    return t


def _reconstructed(cohort, table):
    """Add the reconstructed fibroblast proportion, or return the table unchanged.

    The reconstruction is p_SMA x f_pool, where p_SMA is the alpha-SMA-positive proportion of the
    four-marker pool of Pellinen and colleagues and f_pool the proportion of stromal cells carrying
    neither the T-cell nor the cancer-cell marker. Returns the table and a note saying what happened, so that
    a caller can record in its output why the reconstructed rows are absent when the optional
    folder is not configured.
    """
    folder = paths.pellinen()
    if folder is None:
        return table, "no pellinen folder in the configuration"
    path = paths.pellinen(f"{cohort}_15comb_ratios_alldata.csv")
    if not path.exists():
        return table, f"{path.name} not found in the pellinen folder"
    caf = pd.read_csv(path)
    # The four-marker combinations of the epithelium-masked pool. The Str100 columns are a second
    # normalization of the same combinations and are excluded.
    combos = [c for c in caf.columns
              if c.endswith("PanEpiMask-") and not c.startswith("Str100")]
    for c in combos:
        caf[c] = pd.to_numeric(caf[c], errors="coerce")
    caf["p_SMA"] = caf[[c for c in combos if "SMA+" in c]].sum(axis=1)
    caf["_id"] = pd.to_numeric(caf[caf.columns[0]], errors="coerce")
    key = pd.to_numeric(table.Case, errors="coerce")
    best, best_n = None, -1
    for offset in ID_OFFSETS:
        joined = table.assign(_j=key + offset).merge(
            caf[["_id", "p_SMA"]], left_on="_j", right_on="_id", how="left")
        n = int(joined.p_SMA.notna().sum())
        if n > best_n:
            best, best_n = joined, n
    out = best.drop(columns=["_j", "_id"])
    out["caf_recon"] = out.p_SMA * out.f_pool
    return out, f"{best_n} of {len(table)} patients matched"


def _percase_path(cohort):
    """Where s1 writes (and the later sections read) the per-patient table of a cohort."""
    return paths.intermediate(cohort, "ch1_percase.csv")


def _load_percase(cohort):
    """The table written by s1, which every later section reads instead of recomputing it."""
    path = _percase_path(cohort)
    if not path.exists():
        raise FileNotFoundError(f"{path} is missing. Run chapter1.py --only s1 first.")
    return pd.read_csv(path, dtype={"Case": str, "histology": str, "os_source": str})


def _analysis_set(frame, histology):
    """The rows of one histology that a survival model can use."""
    return frame[frame.histology == histology]


def _sidecar(name, what, columns):
    """Write the sibling .md that describes an output table: a sentence and one line per column."""
    lines = [f"# {name}.csv", "", what, "", "| column | meaning |", "|---|---|"]
    lines += [f"| {c} | {m} |" for c, m in columns] + [""]
    paths.output(CHAPTER, f"{name}.md").write_text("\n".join(lines))


MODEL_COLUMNS = [
    ("HR", "hazard ratio per 1 standard deviation; an interaction row carries the ratio of the two"),
    ("lo", "lower 95 % Wald confidence limit"),
    ("hi", "upper 95 % Wald confidence limit"),
    ("p", "Wald p value"),
    ("n", "patients in the model"),
    ("events", "deaths among them"),
    ("logHR", "the coefficient on the log hazard scale, which is what the pooling averages"),
    ("se", "its standard error"),
    ("p_lr", "likelihood-ratio p, written for an interaction fitted in one cohort only"),
]


def _fmt_p(p):
    """The house rule for a p value in a printed table: three decimals, below that a threshold."""
    return "< 0.001" if float(p) < 0.001 else f"{float(p):.3f}"


def _fmt_p_sign(p):
    """The same with its sign, for the first p of a list: "= 0.016" or "< 0.001"."""
    return "< 0.001" if float(p) < 0.001 else f"= {float(p):.3f}"


# ---------------------------------------------------------------- s1

def s1_per_patient_table():
    """The per-patient table both cohorts share: exposure, composition, reconstruction, clinical.

    Written next to the other patient-level tables rather than into the chapter output folder,
    because it carries one row per patient and the chapter output folder is the shareable part.
    """
    notes = []
    for cohort in paths.COHORTS:
        t = _proportions(cohort)
        t, note = _reconstructed(cohort, t)
        notes.append(f"{cohort}: {note}")
        percase_path = paths.intermediate(cohort, "percase.csv")
        clinical_path = paths.intermediate(cohort, "clinical.csv")
        t = t.merge(clinicalmod.load_clinical(cohort), on="Case", how="inner")
        if len(t) == 0:
            raise ValueError(
                f"merging {percase_path} with {clinical_path} kept no patient for {cohort}")
        path = _percase_path(cohort)
        t.to_csv(path, index=False)
        _write_percase_sidecar(cohort, t, note)
        print(f"s1 {cohort}: {len(t)} patients -> {path}")
    print("s1 reconstruction: " + "; ".join(notes))


def _write_percase_sidecar(cohort, table, note):
    """Write the sibling .md of ch1_percase.csv, one line per column."""
    lines = [f"# ch1_percase.csv ({cohort})", "",
             "One row per patient with at least 20 stromal cells. Written by chapter1.py section s1",
             "from percase.csv (mif/03_percase.py) and clinical.csv. Patient-level: it stays outside the",
             "repository.",
             f"{len(table)} patients. Reconstructed fibroblast proportion: {note}.", "",
             "| column | meaning |", "|---|---|",
             "| Case | patient identifier, as a string |",
             "| n_stroma | stromal cells of the patient, whatever their TrxR1 intensity |",
             "| n_stroma_trxr1 | of those, the ones carrying a TrxR1 call |",
             "| stromal_trxr1 | the TrxR1-positive fraction of the stromal cells that carry a call |"]
    for name in cellsmod.LINEAGES[cohort]:
        lines.append(f"| prop_{name} | {plotting.LABEL[name]} as a share of all stromal cells |")
    lines += ["| f_pool | stromal cells carrying neither the T-cell nor the cancer-cell marker, as a share |"]
    if "p_SMA" in table.columns:
        lines += ["| p_SMA | alpha-SMA-positive share of the Pellinen four-marker pool, from their table |",
                  "| caf_recon | reconstructed fibroblast proportion, p_SMA x f_pool |"]
    lines += ["| histology | LUSC, LUAD or other |",
              "| OS_weeks | overall survival in weeks |",
              "| event | 1 dead, 0 alive at last follow-up |",
              "| stage | 1 to 4 for stage I to IV |",
              "| os_source | BOMI2 only: the clinical source of the follow-up |", ""]
    paths.intermediate(cohort, "ch1_percase.md").write_text("\n".join(lines))


# ---------------------------------------------------------------- s2

def s2_trxr1_survival():
    """Figure 1b: the exposure against overall survival, per cohort and histology, and pooled."""
    rows = []
    fits = {}
    for histology in ("LUSC", "LUAD"):
        per_cohort = []
        for cohort in paths.COHORTS:
            d = _load_percase(cohort)
            r = survival.hr_per_sd(_analysis_set(d, histology), "stromal_trxr1")
            fits[(histology, cohort)] = r
            per_cohort.append(r)
            rows.append(dict(histology=histology, cohort=cohort, **_model_row(r)))
        pooled = survival.meta_fixed(per_cohort)
        rows.append(dict(histology=histology, cohort="pooled", **_model_row(pooled)))
    interactions = []
    for cohort in paths.COHORTS:
        d = _load_percase(cohort)
        i = survival.histology_interaction(d, "stromal_trxr1")
        interactions.append(i)
        rows.append(dict(histology="interaction", cohort=cohort, **_interaction_row(i)))
    pooled = survival.meta_fixed(interactions)
    rows.append(dict(histology="interaction", cohort="pooled", **_model_row(pooled)))
    table = pd.DataFrame(rows)
    table["exposure"] = "stromal TrxR1 fraction"
    table = table[["histology", "cohort", "exposure", "HR", "lo", "hi", "p", "n", "events",
                   "logHR", "se", "p_lr"]]
    path = paths.output(CHAPTER, "fig1b_models.csv")
    table.to_csv(path, index=False)
    _sidecar("fig1b_models",
             "The stromal TrxR1 fraction against overall survival, one row per model. Written by "
             "chapter1.py section s2; Figure 1b is drawn from it.",
             [("histology", "LUSC, LUAD, or interaction for the difference between the two"),
              ("cohort", "BOMI2, BOMI1, or pooled for the fixed-effect combination"),
              ("exposure", "the variable the hazard ratio belongs to")] + MODEL_COLUMNS)
    print(f"s2 -> {path}")
    _panel_b(pd.read_csv(path, float_precision="round_trip"))


def _model_row(r):
    """The columns a fitted or pooled model contributes to a table."""
    return dict(HR=r["HR"], lo=r["lo"], hi=r["hi"], p=r["p"], n=r["n"], events=r["events"],
                logHR=r["logHR"], se=r["se"], p_lr=np.nan)


def _interaction_row(i):
    """The same columns for a histology interaction: the hazard ratio is the ratio of the two.

    The confidence limits are the Wald ones on the product term, so that an interaction row reads
    like any other row of the table; the likelihood-ratio p is carried beside the Wald p as a check.
    """
    return dict(HR=i["HR_ratio"], lo=np.exp(i["beta"] - 1.96 * i["se"]),
                hi=np.exp(i["beta"] + 1.96 * i["se"]), p=i["p_wald"], n=i["n"],
                events=i["events"], logHR=i["beta"], se=i["se"], p_lr=i["p_lr"])


def _panel_b(table):
    """Figure 1b, the forest plot, drawn from fig1b_models.csv and from nothing else."""
    plotting.style()
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    plt.rcParams.update({"ytick.color": plotting.INK, "ytick.labelsize": 7,
                         "xtick.labelsize": 6.3})

    def pick(histology, cohort):
        """The one row of the model table for a histology and cohort."""
        row = table[(table.histology == histology) & (table.cohort == cohort)]
        return row.iloc[0]

    # Squamous rows are squares and adenocarcinoma rows circles (redox.plotting.HM), so that the
    # histology does not rest on color alone; the pooled "combined" rows keep their diamond.
    name = plotting.COHORT_LABEL
    groups = [(title, histology, plotting.HC[histology],
               [(name["BOMI2"], "BOMI2", plotting.HM[histology]),
                (name["BOMI1"], "BOMI1", plotting.HM[histology]),
                ("combined", "pooled", "D")])
              for title, histology in (("squamous", "LUSC"), ("adenocarcinoma", "LUAD"))]
    fig = plt.figure(figsize=(2.6, 2.55))
    ax = fig.add_axes([0.34, 0.27, 0.62, 0.53])
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.tick_params(axis="y", length=0)
    y, ticks, labels = 0.0, [], []
    for title, histology, color, members in groups:
        ax.text(-0.04, y, title, transform=ax.get_yaxis_transform(), ha="right", va="center",
                fontsize=7, fontweight="bold")
        y -= 1
        for label, cohort, marker in members:
            r = pick(histology, cohort)
            filled = r.p < 0.05
            ax.plot([r.lo, r.hi], [y, y], color=color, lw=1.3, solid_capstyle="round", zorder=3)
            ax.plot(r.HR, y, marker, color=color, mfc=color if filled else SURFACE, mec=color,
                    mew=1.0, ms=5.2 if marker == "D" else 4.8, zorder=4)
            ticks.append(y)
            labels.append(label)
            y -= 1
        y -= 0.2
    ax.set_yticks(ticks)
    ax.set_yticklabels(labels)
    ax.set_ylim(y + 0.5, 0.6)
    ax.set_xscale("log")
    ax.set_xlim(0.47, 2.25)
    ax.minorticks_off()
    ax.set_xticks([0.5, 0.7, 1, 1.4, 2])
    ax.set_xticklabels(["0.5", "0.7", "1", "1.4", "2"])
    ax.axvline(1, color=plotting.MUTED, lw=0.7, ls=(0, (2.5, 2.5)), zorder=1)
    ax.set_xlabel("HR per 1 SD of TrxR1+ fraction", fontsize=6.3, color=plotting.INK2, labelpad=2)
    fig.text(0.01, 0.945, "b", fontsize=9.5, fontweight="bold")
    fig.text(0.075, 0.945, "TrxR1 and survival", fontsize=7.5, fontweight="bold")
    # The histology interaction per cohort and pooled, in the format of Figure 5g: the first p
    # carries its sign ("= 0.016" or "< 0.001"), the others are bare.
    b2 = pick("interaction", "BOMI2").p
    b1 = pick("interaction", "BOMI1").p
    both = pick("interaction", "pooled").p
    fig.text(0.96, 0.825, f"histology interaction p {_fmt_p_sign(b2)} ({name['BOMI2']}), "
             f"{_fmt_p(b1)} ({name['BOMI1']}),\n{_fmt_p(both)} (combined)", ha="right",
             va="bottom", fontsize=5.6, color=plotting.INK2, linespacing=1.2)
    fig.legend(handles=[Line2D([], [], marker=plotting.HM["LUSC"], ls="", color=plotting.HC["LUSC"],
                               ms=4.8, label="squamous"),
                        Line2D([], [], marker=plotting.HM["LUAD"], ls="", color=plotting.HC["LUAD"],
                               ms=4.8, label="adenocarcinoma"),
                        Line2D([], [], marker="D", ls="", color=plotting.INK2, ms=4.5,
                               label="combined"),
                        Line2D([], [], marker="o", ls="", color=plotting.INK2, mfc=SURFACE,
                               ms=4.5, label="p ≥ 0.05")],
               loc="lower left", bbox_to_anchor=(0.02, -0.01), ncol=2, frameon=False, fontsize=6.0,
               handletextpad=0.2, columnspacing=0.9)
    for path in plotting.save(fig, CHAPTER, "fig1b_forest"):
        print(f"s2 -> {path}")


# ---------------------------------------------------------------- s3

def s3_composition_survival():
    """Supplementary Table 1b: each lineage proportion of stromal cells against overall survival."""
    rows = []
    for cell in CELLS:
        column = "prop_" + cell
        interactions = []
        for histology in ("LUSC", "LUAD"):
            per_cohort = []
            for cohort in paths.COHORTS:
                d = _load_percase(cohort)
                if column not in d.columns:
                    continue
                r = survival.hr_per_sd(_analysis_set(d, histology), column)
                per_cohort.append(r)
                rows.append(dict(cell=cell, histology=histology, cohort=cohort, **_model_row(r)))
            if len(per_cohort) > 1:
                pooled = survival.meta_fixed(per_cohort)
                rows.append(dict(cell=cell, histology=histology, cohort="pooled",
                                 **_model_row(pooled)))
        for cohort in paths.COHORTS:
            d = _load_percase(cohort)
            if column not in d.columns:
                continue
            i = survival.histology_interaction(d, column)
            interactions.append(i)
            rows.append(dict(cell=cell, histology="interaction", cohort=cohort,
                             **_interaction_row(i)))
        if len(interactions) > 1:
            pooled = survival.meta_fixed(interactions)
            rows.append(dict(cell=cell, histology="interaction", cohort="pooled",
                             **_model_row(pooled)))
    table = pd.DataFrame(rows)
    table["cell_label"] = table.cell.map(plotting.LABEL)
    # Which row the supplementary table prints: the two cohorts pooled where both carry the
    # lineage, BOMI2 alone where they do not.
    both = {c for c in CELLS if ((table.cell == c) & (table.cohort == "BOMI1")).any()}
    table["printed"] = [1 if r.cohort == ("pooled" if r.cell in both else "BOMI2") else 0
                        for r in table.itertuples()]
    table = table[["cell", "cell_label", "histology", "cohort", "printed", "HR", "lo", "hi", "p",
                   "n", "events", "logHR", "se", "p_lr"]]
    path = paths.output(CHAPTER, "supptable1b.csv")
    table.to_csv(path, index=False)
    _sidecar("supptable1b",
             "Each lineage proportion of stromal cells against overall survival. Written by "
             "chapter1.py section s3; Figure 1c is drawn from the rows marked printed.",
             [("cell", "the lineage, T, Fib or Mac"),
              ("cell_label", "how the lineage is named in the figure and the table"),
              ("histology", "LUSC, LUAD, or interaction for the difference between the two"),
              ("cohort", "BOMI2, BOMI1, or pooled for the fixed-effect combination"),
              ("printed", "1 when the supplementary table prints this row")] + MODEL_COLUMNS)
    print(f"s3 -> {path}")
    _panel_c(pd.read_csv(path, float_precision="round_trip"))


def _panel_c(table):
    """Figure 1c, the hazard ratio tiles, drawn from supptable1b.csv."""
    plotting.style()
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    from matplotlib.colors import to_rgb
    printed = table[table.printed == 1]
    fig = plt.figure(figsize=(1.85, 2.05))
    ax = fig.add_axes([0.42, 0.20, 0.52, 0.64])
    ax.axis("off")
    for i, cell in enumerate(CELLS):
        for j, histology in enumerate(("LUSC", "LUAD")):
            r = printed[(printed.cell == cell) & (printed.histology == histology)].iloc[0]
            v = np.log2(r.HR)
            weight = min(abs(v), 1) * 0.8
            base = np.array(to_rgb(NEUTRAL))
            hue = np.array(to_rgb(plotting.HC["LUSC"] if v < 0 else RED))
            ax.add_patch(Rectangle((j + 0.04, -i - 0.46), 0.92, 0.92,
                                   fc=tuple(base + (hue - base) * weight), lw=0))
            ax.text(j + 0.5, -i, f"{r.HR:.2f}", ha="center", va="center", fontsize=6.8,
                    fontweight="bold" if r.p < 0.05 else "normal",
                    color=plotting.INK if r.p < 0.05 else plotting.INK2)
        ax.text(-0.08, -i, plotting.LABEL[cell], ha="right", va="center", fontsize=7)
    for j, label in enumerate(("LUSC", "LUAD")):
        ax.text(j + 0.5, 0.62, label, ha="center", va="bottom", fontsize=6.3, color=plotting.INK2)
    ax.set_xlim(0, 2)
    ax.set_ylim(-2.5, 0.55)
    ax.set_aspect("equal", adjustable="box")   # square tiles whatever the figure size
    fig.text(0.01, 0.94, "c", fontsize=9.5, fontweight="bold")
    fig.text(0.10, 0.94, "Composition and survival", fontsize=7.5, fontweight="bold")
    fig.text(0.42, 0.15, "HR per 1 SD; blue longer survival,\nred shorter; bold p < 0.05",
             fontsize=5.8, color=plotting.INK2, va="top")
    for path in plotting.save(fig, CHAPTER, "fig1c_tiles"):
        print(f"s3 -> {path}")


# ---------------------------------------------------------------- s4

def s4_trxr1_vs_composition():
    """Supplementary Table 1c: the exposure against each proportion, measured and reconstructed."""
    rows = []
    for cohort in paths.COHORTS:
        d = _load_percase(cohort)
        for histology in ("LUSC", "LUAD"):
            x = _analysis_set(d, histology)
            for cell in CELLS:
                column = "prop_" + cell
                if column in x.columns:
                    rows.append(_rho_row(x, column, cohort, histology, cell, "measured"))
                if cell == "Fib" and "caf_recon" in x.columns:
                    rows.append(_rho_row(x, "caf_recon", cohort, histology, cell,
                                         "reconstructed from Pellinen"))
    table = pd.DataFrame(rows)
    table["cell_label"] = table.cell.map(plotting.LABEL)
    # Which row each cell of the printed table holds: the measured proportion where the panel
    # carries the marker, the reconstruction where it does not.
    measured = {(r.cohort, r.cell) for r in table.itertuples() if r.how == "measured"}
    table["printed"] = [1 if (r.how == "measured" or (r.cohort, r.cell) not in measured) else 0
                        for r in table.itertuples()]
    table = table[["cell", "cell_label", "cohort", "histology", "how", "printed", "rho", "p", "n"]]
    path = paths.output(CHAPTER, "supptable1c.csv")
    table.to_csv(path, index=False)
    _sidecar("supptable1c",
             "Spearman correlation between the stromal TrxR1 fraction and each composition "
             "proportion, per cohort and histology. Written by chapter1.py section s4; Figure 1e is "
             "drawn from the rows marked printed.",
             [("cell", "the lineage, T, Fib or Mac"),
              ("cell_label", "how the lineage is named in the figure and the table"),
              ("cohort", "BOMI2 or BOMI1"),
              ("histology", "LUSC or LUAD"),
              ("how", "measured, or reconstructed from Pellinen where the panel has no alpha-SMA"),
              ("printed", "1 when the supplementary table prints this row"),
              ("rho", "Spearman correlation over patients"),
              ("p", "its two-sided p value"),
              ("n", "patients with both values")])
    print(f"s4 -> {path}")
    _panel_e(pd.read_csv(path, float_precision="round_trip"))


def _rho_row(frame, column, cohort, histology, cell, how):
    """One Spearman correlation between the exposure and one composition column."""
    x = frame[["stromal_trxr1", column]].replace([np.inf, -np.inf], np.nan).dropna()
    rho, p = stats.spearman(x.stromal_trxr1, x[column])
    return dict(cell=cell, cohort=cohort, histology=histology, how=how,
                rho=rho, p=p, n=int(len(x)))


def _panel_e(table):
    """Figure 1e, the correlation pies, drawn from supptable1c.csv."""
    plotting.style()
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle, Patch, Wedge
    from matplotlib.colors import to_rgb
    printed = table[table.printed == 1]
    columns = [(h, c) for h in ("LUSC", "LUAD") for c in paths.COHORTS]

    def tint(color, weight):
        """The color blended from the surface color toward color by weight (0 to 1)."""
        a, b = np.array(to_rgb(SURFACE)), np.array(to_rgb(color))
        return tuple(a + (b - a) * weight)

    fig = plt.figure(figsize=(3.3, 2.05))
    ax = fig.add_axes([0.26, 0.20, 0.70, 0.60])
    ax.axis("off")
    ax.set_xlim(0, 4)
    ax.set_ylim(-2.5, 0.55)
    ax.set_aspect("equal", adjustable="box")   # round pies whatever the figure size
    radius = 0.3
    for i, cell in enumerate(CELLS):
        ax.text(-0.08, -i, plotting.LABEL[cell], ha="right", va="center", fontsize=7)
        for j, (histology, cohort) in enumerate(columns):
            row = printed[(printed.cell == cell) & (printed.cohort == cohort)
                          & (printed.histology == histology)]
            cx, cy = j + 0.5, -i + 0.1
            if row.empty:
                ax.text(cx, -i, "not\nmeasured", ha="center", va="center", fontsize=5.5,
                        color=plotting.MUTED)
                continue
            r = row.iloc[0]
            if pd.isna(r.rho):
                # A Pellinen folder that is configured but matches no patient of this cohort:
                # the row exists but carries no correlation to draw.
                ax.text(cx, -i, "not\navailable", ha="center", va="center", fontsize=5.5,
                        color=plotting.MUTED)
                continue
            strong = r.p < 0.05
            hue = plotting.LC["T"] if r.rho >= 0 else plotting.LC["Fib"]
            ax.add_patch(Circle((cx, cy), radius, fc=DISC, ec=plotting.GRID, lw=0.6))
            extent = 360 * abs(r.rho)
            angles = (90 - extent, 90) if r.rho >= 0 else (90, 90 + extent)
            ax.add_patch(Wedge((cx, cy), radius, *angles,
                               fc=tint(hue, 0.55 if strong else 0.22),
                               ec=tint(hue, 0.9 if strong else 0.45), lw=0.7))
            ax.text(cx, cy - 0.43, f"{r.rho:+.2f}", ha="center", va="center", fontsize=5.8,
                    color=plotting.INK if strong else plotting.MUTED)
    for j, (histology, cohort) in enumerate(columns):
        ax.text(j + 0.5, 0.62, f"{histology}\n{plotting.COHORT_LABEL[cohort]}", ha="center",
                va="bottom", fontsize=6.3, color=plotting.INK2)
    fig.text(0.01, 0.94, "e", fontsize=9.5, fontweight="bold")
    fig.text(0.06, 0.94, "TrxR1 and composition", fontsize=7.5, fontweight="bold")
    fig.text(0.26, 0.15, "filled part of circle = Spearman correlation", fontsize=5.8,
             color=plotting.INK2, va="top")
    fig.legend(handles=[Patch(fc=tint(plotting.LC["T"], 0.55), ec=tint(plotting.LC["T"], 0.9),
                              lw=0.7, label="positive"),
                        Patch(fc=tint(plotting.LC["Fib"], 0.55), ec=tint(plotting.LC["Fib"], 0.9),
                              lw=0.7, label="negative"),
                        Patch(fc=tint(plotting.LC["T"], 0.22), ec=tint(plotting.LC["T"], 0.45),
                              lw=0.7, label="paler p ≥ 0.05")],
               loc="upper left", bbox_to_anchor=(0.24, 0.12), ncol=3, frameon=False, fontsize=5.8,
               handlelength=1.0, handleheight=0.9, handletextpad=0.4, columnspacing=0.9)
    for path in plotting.save(fig, CHAPTER, "fig1e_pies"):
        print(f"s4 -> {path}")


# ---------------------------------------------------------------- s5

def s5_numbers():
    """The Chapter 1 numbers that no other table holds: composition means and the Figure 1a patient sums.

    The first paragraph of "Stromal composition and survival" quotes, for Uppsala II, the mean over all
    patients of the fibroblast, T-cell and macrophage proportions of stromal cells, and the proportion of
    stromal cells assigned to none of them: the stromal cells without a lineage marker together with
    the cytokeratin-positive cells of the stromal compartment (prop_Other + prop_CK; the four proportions
    sum to one per patient). Figure 1a draws the patients of the survival analysis of Figure 1b per cohort,
    LUSC and LUAD summed; they are read back from fig1b_models.csv, written by s2.
    """
    source = paths.output(CHAPTER, "fig1b_models.csv")
    if not source.exists():
        raise FileNotFoundError(f"{source} is missing. Run chapter1.py --only s2 first.")
    rows = []
    d = _load_percase("BOMI2")
    note = "mean over patients of the per-patient proportion"
    for cell in ("Fib", "T", "Mac"):
        rows.append(dict(key=f"mean_prop_{cell}_percent",
                         item=f"mean {plotting.LABEL[cell]} proportion of stromal cells in percent",
                         cohort="BOMI2", subset="all_patients", value=100 * d["prop_" + cell].mean(),
                         n=int(len(d)), note=note))
    rows.append(dict(key="mean_prop_unassigned_percent",
                     item="mean proportion of stromal cells not assigned to T cells, fibroblasts or "
                          "macrophages, in percent",
                     cohort="BOMI2", subset="all_patients",
                     value=100 * (d.prop_Other + d.prop_CK).mean(), n=int(len(d)),
                     note=note + "; stromal cells without a lineage marker (prop_Other) and "
                                 "cytokeratin-positive stromal cells (prop_CK)"))
    models = pd.read_csv(source)
    for cohort in paths.COHORTS:
        m = models[(models.cohort == cohort) & models.histology.isin(["LUSC", "LUAD"])]
        rows.append(dict(key="fig1a_patients", item="patients drawn in Figure 1a", cohort=cohort,
                         subset="LUSC and LUAD, the survival models of Figure 1b",
                         value=int(m.n.sum()), n=int(m.n.sum()),
                         note="the patients of the LUSC and LUAD rows of fig1b_models.csv, summed"))
    table = pd.DataFrame(rows)
    path = paths.output(CHAPTER, "ch1_numbers.csv")
    table.to_csv(path, index=False)
    _sidecar("ch1_numbers",
             "The Chapter 1 numbers that no other table holds: the composition means of the first "
             "paragraph of \"Stromal composition and survival\" and the patient sums drawn in Figure 1a. "
             "Written by chapter1.py section s5.",
             [("key", "a short name for the quantity, used by manifest/check.py"),
              ("item", "the same quantity in words"),
              ("cohort", "BOMI2 or BOMI1"),
              ("subset", "the patients it is computed over"),
              ("value", "the mean in percent, or a count of patients"),
              ("n", "patients"),
              ("note", "how the value is defined")])
    print(f"s5 -> {path}")
    print(table.to_string(index=False))


# ---------------------------------------------------------------- s6

# How the staining check names each composition variable.
VARIABLE = {"T": "T-cell proportion of stromal cells", "Fib": "fibroblast proportion of stromal cells",
            "Mac": "macrophage proportion of stromal cells"}


def s6_staining_check():
    """Supplementary Table 1c caption: the Figure 1e correlations within slide and without background.

    TrxR1 positivity differs between tissue microarray slides, and background staining could move the
    TrxR1-positive fraction and the composition together. Each printed correlation of
    supptable1c.csv is therefore recomputed three ways on the patients with a slide and a stromal
    background index: raw, within slide (the ranks of both variables regressed on slide
    indicators), and within slide with the rank of the stromal background index as a further
    covariate. The raw value on all patients must reproduce supptable1c.csv exactly.
    """
    source = paths.output(CHAPTER, "supptable1c.csv")
    if not source.exists():
        raise FileNotFoundError(f"{source} is missing. Run chapter1.py --only s4 first.")
    printed = pd.read_csv(source, float_precision="round_trip")
    printed = printed[printed.printed == 1]
    rows = []
    for cohort in paths.COHORTS:
        cells = staining.cells_with_intensity(cohort, columns=["core", "Case", "tissue", "TrxR1"])
        _, background = staining.fraction_and_background(cells[cells.tissue == "Stroma"])
        slide = staining.patient_slide(cells, cohort)
        extra = pd.DataFrame({"bg": background, "slide": slide}).rename_axis("Case").reset_index()
        d = _load_percase(cohort).merge(extra, on="Case", how="left")
        print(f"s6 {cohort}: {int(d.bg.notna().sum())} of {len(d)} patients with a background index, "
              f"{int(d.slide.nunique())} slides")
        for histology in ("LUSC", "LUAD"):
            for cell in CELLS:
                sel = printed[(printed.cohort == cohort) & (printed.histology == histology)
                              & (printed.cell == cell)]
                if sel.empty:
                    continue
                r = sel.iloc[0]
                column = "prop_" + cell if r.how == "measured" else "caf_recon"
                x0 = _analysis_set(d, histology)[["stromal_trxr1", column, "bg", "slide"]]
                x0 = x0.replace([np.inf, -np.inf], np.nan).dropna(subset=["stromal_trxr1", column])
                rho0, p0 = stats.spearman(x0.stromal_trxr1, x0[column])
                if abs(rho0 - r.rho) > 1e-12 or len(x0) != r.n:
                    raise ValueError(f"s6 {cohort} {histology} {cell}: rho {rho0} on {len(x0)} "
                                     f"patients, supptable1c.csv {r.rho} on {r.n}")
                x = x0.dropna(subset=["bg", "slide"])
                a, b = x.stromal_trxr1.to_numpy(), x[column].to_numpy()
                S = staining.slide_indicators(x.slide)
                rho, p = stats.spearman(a, b)
                r1, p1 = staining.partial_rank_correlation(a, b, S)
                r2, p2 = staining.partial_rank_correlation(a, b, S + [rankdata(x.bg)])
                rows.append(dict(cohort=cohort, histology=histology, cell=cell,
                                 variable=VARIABLE[cell], how=r.how, n_fig1e=len(x0), rho_fig1e=rho0,
                                 p_fig1e=p0, n=len(x), n_slides=x.slide.nunique(), rho_raw=rho,
                                 p_raw=p, rho_within_slide=r1, p_within_slide=p1,
                                 rho_slide_background=r2, p_slide_background=p2))
    table = pd.DataFrame(rows)
    path = paths.output(CHAPTER, "supptable1c_staining.csv")
    table.to_csv(path, index=False)
    _sidecar("supptable1c_staining",
             "The correlations of Figure 1e and Supplementary Table 1c, recomputed within tissue "
             "microarray slide and with the background staining of stromal cells as a covariate. "
             "Written by chapter1.py section s6; the Supplementary Table 1c caption quotes the last "
             "column pair.",
             [("cohort", "BOMI2 or BOMI1"),
              ("histology", "LUSC or LUAD"),
              ("cell", "the lineage, T, Fib or Mac"),
              ("variable", "the composition variable in words"),
              ("how", "measured, or reconstructed from Pellinen where the panel has no alpha-SMA"),
              ("n_fig1e", "patients of Figure 1e (supptable1c.csv)"),
              ("rho_fig1e", "Spearman rho on them, equal to supptable1c.csv"),
              ("p_fig1e", "its p"),
              ("n", "of those, the patients with a slide and a stromal background index (at least "
                    "10 TrxR1-negative stromal cells)"),
              ("n_slides", "TMA slides among them"),
              ("rho_raw", "Spearman rho on these patients"),
              ("p_raw", "its p"),
              ("rho_within_slide", "correlation of the ranks after slide indicators are regressed out"),
              ("p_within_slide", "t-test p, n - k - 2 degrees of freedom"),
              ("rho_slide_background", "the same with the rank of the stromal background index added"),
              ("p_slide_background", "its p")])
    print(f"s6 -> {path}")


# ---------------------------------------------------------------- s7

# The three rows s7 adds to Supplementary Table 1c, in the order of the table. The labels write mm2
# for the square millimeter, as the other package tables do.
DENSITY_ROWS = [("dens_stroma", "T cells per mm2 of stroma"),
                ("dens_tumor", "T cells per mm2 of tumor compartment (islands)"),
                ("ts_ratio", "tumor-to-stroma ratio of T-cell density")]

# The definition of a T cell in s7, with the printed flag: the T-cell lineage (CD3-positive, first in
# the priority order, whatever the other markers), which the T-cell proportion rows of the table use
# as well.
DENSITY_DEFINITIONS = [("", "T lineage", 1)]


def _p_table1c(p):
    """p as Supplementary Table 1c prints it: below 0.001 a bound, 3 decimals below 0.1, else 2."""
    p = float(p)
    return "< 0.001" if p < 0.001 else (f"{p:.3f}" if p < 0.1 else f"{p:.2f}")


def _tcell_density(cohort):
    """Per patient of ch1_percase.csv: T cells per mm2 of stroma and of tumor compartment.

    Counts and areas are summed over the patient's cores before the division, as percase.csv does
    for its stromal densities (mif/03_percase.py): the cores are those of cells.tsv, and an area row
    is matched to its patient through the cell table. A density with no area (or an area of zero)
    and a tumor-to-stroma ratio with no T cell in the stroma are left empty. The stromal density
    must equal dens_T_S of percase.csv, or the run stops.
    """
    cells = cellsmod.load_cells(cohort, columns=["core", "Case", "tissue", "lineage"])
    areas = pd.read_csv(paths.intermediate(cohort, "areas.csv"), dtype={"core": str, "tissue": str})
    case_of_core = cells.drop_duplicates("core").set_index("core")["Case"]
    areas = areas.assign(Case=areas.core.map(case_of_core))
    areas = areas[areas.Case.notna()]
    d = _load_percase(cohort)[["Case", "histology", "stromal_trxr1"]].set_index("Case")
    area = {"stroma": areas[areas.tissue == "Stroma"].groupby("Case").area_mm2.sum(),
            "tumor": areas[areas.tissue == "Tumor"].groupby("Case").area_mm2.sum()}
    tissue = {"stroma": "Stroma", "tumor": "Tumor"}
    t = cells[cells.lineage == "T"]
    for part in ("stroma", "tumor"):
        d[f"area_{part}_mm2"] = area[part].reindex(d.index)
    for suffix, _, _ in DENSITY_DEFINITIONS:
        group = t
        for part in ("stroma", "tumor"):
            n = group[group.tissue == tissue[part]].groupby("Case").size()
            d[f"n_T_{part}{suffix}"] = n.reindex(d.index).fillna(0).astype(int)
            dens = d[f"n_T_{part}{suffix}"] / area[part].reindex(d.index)
            d[f"dens_{part}{suffix}"] = dens.replace([np.inf, -np.inf], np.nan)
        stroma = d[f"dens_stroma{suffix}"]
        if bool(((stroma == 0) & (d[f"dens_tumor{suffix}"] > 0)).any()):
            raise ValueError(f"s7 {cohort}: a patient has T cells in the tumor compartment but none in the stroma")
        d[f"ts_ratio{suffix}"] = (d[f"dens_tumor{suffix}"] / stroma).where(stroma > 0)
    # The stromal density is the dens_T_S column of percase.csv; the two must agree.
    ref = pd.read_csv(paths.intermediate(cohort, "percase.csv"), usecols=["Case", "dens_T_S"],
                      dtype={"Case": str}, float_precision="round_trip").set_index("Case").dens_T_S
    ref = ref.reindex(d.index)
    gap = (d.dens_stroma - ref).abs() / ref.abs()
    if bool((d.dens_stroma.isna() != ref.isna()).any()) or float(gap.max()) > 1e-12:
        raise ValueError(f"s7 {cohort}: T cells per mm2 of stroma differ from dens_T_S of percase.csv "
                         f"(largest relative difference {float(gap.max())})")
    return d.reset_index()


def s7_tcell_density():
    """Supplementary Table 1c: T cells per mm2 of stroma and of tumor compartment, and their tumor-to-stroma ratio.

    On the patients of the other rows of the table (ch1_percase.csv, LUSC and LUAD, checked against
    the n of the T-cell rows of supptable1c.csv), the Spearman correlation of the TrxR1-positive
    fraction of stromal cells with the T-cell density of each compartment and with their
    tumor-to-stroma ratio (the tumor-compartment density over the stromal density). No random
    numbers are used.
    """
    source = paths.output(CHAPTER, "supptable1c.csv")
    if not source.exists():
        raise FileNotFoundError(f"{source} is missing. Run chapter1.py --only s4 first.")
    t1c = pd.read_csv(source, float_precision="round_trip")
    t1c = t1c[(t1c.cell == "T") & (t1c.how == "measured")]
    rows = []
    for cohort in paths.COHORTS:
        d = _tcell_density(cohort)
        path = paths.intermediate(cohort, "ch1_tcell_density.csv")
        d.to_csv(path, index=False)
        _write_density_sidecar(cohort, d)
        print(f"s7 {cohort}: {len(d)} patients -> {path}")
        d = pd.read_csv(path, dtype={"Case": str, "histology": str}, float_precision="round_trip")
        for histology in ("LUSC", "LUAD"):
            x = _analysis_set(d, histology)
            ref = t1c[(t1c.cohort == cohort) & (t1c.histology == histology)]
            if len(ref) != 1 or int(ref.n.iloc[0]) != len(x):
                raise ValueError(f"s7 {cohort} {histology}: {len(x)} patients, the T-cell row of "
                                 f"{source.name} has {ref.n.tolist()}")
            for suffix, definition, printed in DENSITY_DEFINITIONS:
                for measure, label in DENSITY_ROWS:
                    z = x[["stromal_trxr1", measure + suffix]].replace([np.inf, -np.inf], np.nan).dropna()
                    rho, p = stats.spearman(z.stromal_trxr1, z[measure + suffix])
                    rows.append(dict(measure=measure, measure_label=label, definition=definition,
                                     printed=printed, cohort=cohort, histology=histology, rho=rho, p=p,
                                     n=int(len(z)), n_patients=int(len(x)),
                                     n_undefined=int(len(x) - len(z)),
                                     shown=f"{rho:+.2f}, p {_p_table1c(p)}, n {len(z)}"))
    table = pd.DataFrame(rows)
    order = {m: i for i, (m, _) in enumerate(DENSITY_ROWS)}
    table = table.assign(_d=1 - table.printed, _m=table.measure.map(order),
                         _h=table.histology.map({"LUSC": 0, "LUAD": 1}),
                         _c=table.cohort.map({c: i for i, c in enumerate(paths.COHORTS)}))
    table = table.sort_values(["_d", "_m", "_h", "_c"]).drop(columns=["_d", "_m", "_h", "_c"])
    path = paths.output(CHAPTER, "supptable1c_tcell_density.csv")
    table.to_csv(path, index=False)
    _sidecar("supptable1c_tcell_density",
             "Supplementary Table 1c, the T-cell density rows: Spearman correlation across "
             "patients between the TrxR1-positive fraction of stromal cells and the T-cell density of "
             "the stroma and of the tumor compartment (the cancer cell islands), and the tumor-to-stroma "
             "ratio of T-cell density, on the patients of the other rows (at least 20 stromal cells). Counts and "
             "areas are summed over each patient's cores; a T cell is the package's T-cell lineage "
             "(CD3-positive, first in the priority order). Written by chapter1.py section s7; rows in "
             "the order of the table.",
             [("measure", "dens_stroma (T cells per mm2 of stroma), dens_tumor (per mm2 of tumor "
                          "compartment) or ts_ratio (dens_tumor over dens_stroma)"),
              ("measure_label", "the row label proposed for the table"),
              ("definition", "T lineage"),
              ("printed", "1 when the supplementary table prints this row"),
              ("cohort", "BOMI2 (Uppsala II) or BOMI1 (Uppsala I)"),
              ("histology", "LUSC or LUAD"),
              ("rho", "Spearman correlation over patients"),
              ("p", "its two-sided p value"),
              ("n", "patients with both values"),
              ("n_patients", "patients of the other rows of the table in this cohort and histology"),
              ("n_undefined", "of those, the patients without a value: a tumor-to-stroma ratio is "
                              "undefined when the patient has no T cell in the stroma"),
              ("shown", "the cell as the table prints it: rho, p (3 decimals below 0.1, 2 above, "
                        "< 0.001 below 0.001), n")])
    print(f"s7 -> {path}")
    print(table[table.printed == 1][["measure", "cohort", "histology", "shown", "n_undefined"]]
          .to_string(index=False))


def _write_density_sidecar(cohort, table):
    """The sibling description of the per-patient density table."""
    lines = [f"# ch1_tcell_density.csv ({cohort})", "",
             "One row per patient of ch1_percase.csv. Written by chapter1.py section s7 from cells.tsv,",
             "areas.csv, percase.csv and ch1_percase.csv. Patient-level: it stays outside the repository.",
             f"{len(table)} patients. A T cell is the T-cell lineage of cells.tsv (CD3-positive, first in the",
             "priority order). Counts and areas are summed over the patient's cores.", "",
             "| column | meaning |", "|---|---|",
             "| Case | patient identifier, as a string |",
             "| histology | LUSC, LUAD or other |",
             "| stromal_trxr1 | the TrxR1-positive fraction of the stromal cells that carry a call |",
             "| area_stroma_mm2 | segmented stroma, summed over the patient's cores |",
             "| area_tumor_mm2 | segmented tumor compartment, summed over the patient's cores |",
             "| n_T_stroma | T cells of the Stroma category |",
             "| dens_stroma | n_T_stroma per mm2 of stroma, equal to dens_T_S of percase.csv |",
             "| n_T_tumor | T cells of the Tumor category |",
             "| dens_tumor | n_T_tumor per mm2 of tumor compartment |",
             "| ts_ratio | dens_tumor over dens_stroma; empty without stromal T cells |", ""]
    paths.intermediate(cohort, "ch1_tcell_density.md").write_text("\n".join(lines))


# ---------------------------------------------------------------- run

SECTIONS = {"s1": s1_per_patient_table,
            "s2": s2_trxr1_survival,
            "s3": s3_composition_survival,
            "s4": s4_trxr1_vs_composition,
            "s5": s5_numbers,
            "s6": s6_staining_check,
            "s7": s7_tcell_density}


def main():
    """Command line: run the sections named by --only, or all of them in order."""
    ap = argparse.ArgumentParser(description=__doc__.split("RUN:")[0].strip(),
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", nargs="*", choices=sorted(SECTIONS),
                    help="run these sections only; the default runs all of them in order")
    a = ap.parse_args()
    for name in (a.only or sorted(SECTIONS)):
        SECTIONS[name]()


if __name__ == "__main__":
    main()
