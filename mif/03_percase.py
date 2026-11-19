#!/usr/bin/env python3
"""
The per-patient table every chapter starts from.

WHAT: counts a patient's stromal cells by lineage, turns the counts into densities with the
      patient's segmented stroma area, adds the stromal composition and the TrxR1 fraction that is
      the paper's exposure, and joins the de-identified clinical columns.
      It produces no figure or table of the paper itself; it is the input of chapter1.py (Figure 1,
      Supplementary Table 1b and 1c), chapter3.py, chapter4.py and tumor_vs_stroma.py.
FROM: work_root/intermediate/<cohort>/{cells.tsv, areas.csv, clinical.csv}, written by
      mif/02_cells.R and mif/01_clinical.py.
TO:   work_root/intermediate/<cohort>/percase.csv, with a sibling percase.md that documents the
      columns. The stromal densities are written as dens_<lineage>_S only (stromal cells of a
      lineage per square millimeter of segmented stroma).
WHY:  one table holds the counts, the proportions, the densities and the exposure, computed once,
      so that no chapter script recomputes the TrxR1 fraction from the cell table by itself.

Proportions and the exposure are stated with their denominator: the exposure here averages the
TrxR1 call over every stromal cell that carries one, with no lineage gating (the fraction of all
stromal cells that are TrxR1-positive), and each proportion has all stromal cells of the patient
as its denominator.

RUN: python3 mif/03_percase.py [--cohort BOMI2|BOMI1|both]
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from redox import cells as rc  # noqa: E402
from redox import paths  # noqa: E402

CLINICAL = ["histology", "OS_weeks", "event", "stage", "os_source"]

# How the lineage count of a cohort is spelled out in the sidecar note below.
COUNT_WORDS = {2: "two", 3: "three", 4: "four", 5: "five", 6: "six"}

DOC = """# percase.csv ({cohort})

One row per patient of the analyzed cohort, {n_rows} patients. Written by mif/03_percase.py from
cells.tsv, areas.csv and clinical.csv.

| column | meaning |
|---|---|
| Case | patient identifier, as a string |
| n_stroma | cells of the Stroma tissue category, every lineage, including those whose TrxR1 intensity is missing |
| n_stroma_trxr1 | of those, the cells that carry a TrxR1 call |
| stromal_trxr1 | mean of the TrxR1 call over the stromal cells that carry one; the exposure of Chapters 1 and 4 |
{prop_rows}| dens_LINEAGE_S | cells of that lineage in the Stroma category per square millimeter of segmented stroma |
| histology | LUSC, LUAD or other |
| OS_weeks | overall survival in weeks from diagnosis |
| event | 1 dead, 0 alive at last follow-up, empty when the follow-up is missing |
| stage | 1 to 4 for stage I to IV |
| os_source | BOMI2: ENCIMA or Swedish_2015; empty for BOMI1 |

Lineages of this cohort: {lineages}. The lineage of a cell is the first marker of the priority
order it is positive for, and Other is the cell that is positive for none of them.

Every proportion has all stromal cells of the patient as its denominator, so the {n_lineages} of
them sum to one. A density is a summed count over a summed area, both taken over the patient's
cores, not the mean of the per-core densities.
"""


def per_case(cohort):
    """Counts, proportions, densities, the exposure and the clinical columns for one cohort."""
    lineages = rc.LINEAGES[cohort]

    cells_path = paths.intermediate(cohort, "cells.tsv")
    cells = pd.read_csv(cells_path, sep="\t", usecols=["core", "Case", "tissue", "lineage", "TrxR1"],
                        dtype={"core": str, "Case": str, "tissue": str, "lineage": str,
                               "TrxR1": "Int8"})
    areas = pd.read_csv(paths.intermediate(cohort, "areas.csv"),
                        dtype={"core": str, "tissue": str})
    print(f"[{cohort}] {len(cells)} cells, {cells.Case.nunique()} patients, "
          f"{len(areas)} area rows", flush=True)

    unknown = set(cells.lineage) - set(lineages)
    if unknown:
        raise ValueError(f"{cells_path} holds lineages the cohort does not have: {sorted(unknown)}")

    case_of_core = cells.drop_duplicates("core").set_index("core")["Case"]
    areas = areas.assign(Case=areas.core.map(case_of_core))
    orphan = int(areas.Case.isna().sum())
    if orphan:
        print(f"[{cohort}] {orphan} area rows name a core the cell table does not have, dropped",
              flush=True)
        areas = areas[areas.Case.notna()]

    # Segmented stroma area per patient, summed over the patient's cores.
    area_stroma = areas[areas.tissue == "Stroma"].groupby("Case").area_mm2.sum()

    stroma = cells[cells.tissue == "Stroma"]
    # One row per patient of the cell table, also for a patient without stromal cells (count 0).
    out = pd.DataFrame({"n_stroma": stroma.groupby("Case").size()},
                       index=cells.groupby("Case").size().index)
    with_call = stroma[stroma.TrxR1.notna()]
    out["n_stroma_trxr1"] = with_call.groupby("Case").size()
    out[["n_stroma", "n_stroma_trxr1"]] = out[["n_stroma", "n_stroma_trxr1"]].fillna(0).astype(int)

    # The exposure comes from the shared module, so that this table and every chapter that
    # recomputes it cannot drift apart.
    exposure = rc.stromal_trxr1_fraction(cells).set_index("Case")
    disagree = int((out.n_stroma.reindex(exposure.index).astype("Int64")
                    != exposure.n_stroma.astype("Int64")).sum())
    if disagree or len(exposure) != len(out):
        raise ValueError(f"the stromal cell count here and in redox.cells disagree for "
                         f"{disagree} patients, and the two tables hold {len(out)} and "
                         f"{len(exposure)} patients; redox.cells.stromal_trxr1_fraction drops "
                         f"patients below MIN_STROMAL_CELLS ({rc.MIN_STROMAL_CELLS} stromal cells), "
                         f"which is the most likely cause of a patient missing from one table")
    out["stromal_trxr1"] = exposure["stromal_trxr1"]

    # Proportions of the stromal cells, and stromal densities.
    counted = stroma.groupby(["Case", "lineage"]).size().unstack(fill_value=0)
    for name in lineages:
        out[f"prop_{name}"] = counted.get(name, 0) / out.n_stroma
    for name in lineages:
        n = stroma[stroma.lineage == name].groupby("Case").size()
        out[f"dens_{name}_S"] = n.reindex(out.index).fillna(0) / area_stroma

    out.index.name = "Case"
    out = out.reset_index()
    # event and stage are whole numbers with a patient or two missing, so they are read onto the
    # nullable integer type; the default parser would make them floats and write 1.0 for a death.
    clinical = pd.read_csv(paths.intermediate(cohort, "clinical.csv"),
                           dtype={"Case": str, "event": "Int8", "stage": "Int8"})
    out = out.merge(clinical[["Case"] + CLINICAL], on="Case", how="left")
    missing = int(out.histology.isna().sum())
    if missing:
        raise ValueError(f"{missing} patients of the cell table are not in clinical.csv")
    return out.sort_values("Case", key=lambda s: s.astype(str)).reset_index(drop=True)


def run(cohort):
    """Build percase.csv and percase.md for one cohort and print a one-line check."""
    out = per_case(cohort)
    target = paths.intermediate(cohort, "percase.csv")
    out.to_csv(target, index=False)
    prop_rows = "".join(
        f"| prop_{n} | this lineage as a proportion of all stromal cells of the patient |\n"
        for n in rc.LINEAGES[cohort])
    n_lineages = len(rc.LINEAGES[cohort])
    (target.parent / "percase.md").write_text(DOC.format(
        cohort=cohort, n_rows=len(out), prop_rows=prop_rows,
        lineages=", ".join(rc.LINEAGES[cohort]),
        n_lineages=COUNT_WORDS.get(n_lineages, str(n_lineages))))
    share = out[[f"prop_{n}" for n in rc.LINEAGES[cohort]]].sum(axis=1)
    print(f"[{cohort}] wrote {target.name} ({len(out)} patients, {out.shape[1]} columns); "
          f"proportions sum to {share.min():.12f} at least and {share.max():.12f} at most; "
          f"stromal TrxR1 fraction median {out.stromal_trxr1.median():.4f}", flush=True)


def main():
    """Command line: run one cohort or both."""
    ap = argparse.ArgumentParser(description=__doc__.split("RUN:")[0].strip(),
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cohort", choices=["BOMI2", "BOMI1", "both"], default="both")
    a = ap.parse_args()
    for cohort in (paths.COHORTS if a.cohort == "both" else (a.cohort,)):
        run(cohort)
    return 0


if __name__ == "__main__":
    sys.exit(main())
