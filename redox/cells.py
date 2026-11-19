"""
The per-cell table, the core qualifying rules and the paper's exposure.

WHAT: loads cells.tsv, decides which cores qualify for a chapter, and computes the lineage-free
      stromal TrxR1 fraction per patient.
FROM: work_root/intermediate/<cohort>/cells.tsv, one row per Tumor or Stroma cell of a primary core.
TO:   DataFrames, and a set of core identifiers.
WHY:  the exposure of the whole paper and the two core qualifying rules are each defined once, so
      that Chapter 1 and Chapter 4 cannot drift apart.

cells.tsv columns: core, Case, tissue, x_raw, y_raw, CD3, CD68, ASMA, CK, lineage, TrxR1
(BOMI1 carries CD3 and CK only, and its lineage is T, CK or Other). x_raw and y_raw are the export
positions (BOMI2 pixels, BOMI1 micrometers); load_cells returns them as x_um and y_um, converted in
Python with redox.coords.to_um, in double precision.
"""
import pandas as pd

from . import coords, paths

# Lineage names in priority order, last entry the fallback. BOMI2 assigns CD3 before CD68 before
# CK before alpha-SMA; BOMI1 has no CD68 and no alpha-SMA channel, so its "Other" holds the
# fibroblasts and macrophages that BOMI2 separates.
LINEAGES = {"BOMI2": ["T", "Mac", "CK", "Fib", "Other"],
            "BOMI1": ["T", "CK", "Other"]}

# Marker columns per cohort, in priority order, with the lineage each one names.
MARKERS = {"BOMI2": [("CD3", "T"), ("CD68", "Mac"), ("CK", "CK"), ("ASMA", "Fib")],
           "BOMI1": [("CD3", "T"), ("CK", "CK")]}

# Core qualifying rules (RULE_CH3 for chapter3.py, RULE_CH4 for chapter4.py): minimum stromal cells, TrxR1-positive stromal cells, TrxR1-negative ones.
RULE_CH3 = (50, 5, 5)
RULE_CH4 = (60, 25, 25)

# Minimum stromal cells for a patient to carry the Chapter 1 exposure.
MIN_STROMAL_CELLS = 20


def load_cells(cohort, tissue=None, columns=None):
    """cells.tsv for a cohort.

    tissue   "Stroma" or "Tumor" keeps only that inForm tissue category; none keeps both.
    columns  the columns to return; none returns all of them. The tissue column is read whenever
             it is needed for the filter and dropped again when it was not asked for.

    The file stores the export positions as x_raw and y_raw. Asking for x_um or y_um (or for every
    column) reads them and returns x_um, y_um = redox.coords.to_um(raw, cohort), computed here in
    double precision from the export values. With every column
    the raw pair is replaced by the micrometer pair; x_raw and y_raw are returned only when asked
    for by name.
    """
    path = paths.intermediate(cohort, "cells.tsv")
    if not path.exists():
        raise FileNotFoundError(
            f"{path} is missing. Run mif/02_cells.R and mif/03_percase.py first.")
    raw_of = {"x_um": "x_raw", "y_um": "y_raw"}
    if columns is None:
        want, um, keep_raw = None, ["x_um", "y_um"], []
    else:
        um = [c for c in columns if c in raw_of]
        keep_raw = [c for c in columns if c in raw_of.values()]
        want = [raw_of.get(c, c) for c in columns] + (["tissue"] if tissue is not None else [])
        want = list(dict.fromkeys(want))
    # float_precision "round_trip" makes pandas parse a position back to the exact double that was
    # written; the default parser can be one last bit off.
    d = pd.read_csv(path, sep="\t", usecols=want, float_precision="round_trip",
                    dtype={"core": str, "Case": str, "tissue": str, "lineage": str,
                           "CD3": "Int8", "CD68": "Int8", "ASMA": "Int8", "CK": "Int8",
                           "TrxR1": "Int8"})
    if any(raw_of[c] not in d.columns for c in um):
        raise KeyError(f"{path} has no x_raw and y_raw columns. It predates the raw-position "
                       "schema of the package; rerun mif/02_cells.R.")
    # The micrometer pair takes the place of the raw pair, so the columns keep the file's order.
    um_of = {v: k for k, v in raw_of.items()}
    order = []
    for c in d.columns:
        if c in um_of and um_of[c] in um:
            d[um_of[c]] = coords.to_um(d[c].to_numpy(dtype=float), cohort)
            order.append(um_of[c])
            if c in keep_raw:
                order.append(c)
        else:
            order.append(c)
    d = d[order]
    if tissue is not None:
        d = d[d.tissue == tissue]
        if columns is not None and "tissue" not in columns:
            d = d.drop(columns=["tissue"])
    return d.reset_index(drop=True)


def _stromal_with_call(cells):
    """Stromal cells that carry a TrxR1 call and a position, the set every core rule counts on."""
    d = cells
    if "tissue" in d.columns:
        d = d[d.tissue == "Stroma"]
    d = d[d.TrxR1.notna()]
    xy = [c for c in ("x_um", "y_um") if c in d.columns]
    if xy:
        d = d.dropna(subset=xy)
    return d


def qualifying_cores(cells, min_stromal, min_pos, min_neg):
    """The cores that pass a rule, as a set of core identifiers.

    A core qualifies when it holds at least min_stromal stromal cells with a TrxR1 call, of which
    at least min_pos are TrxR1-positive and at least min_neg are TrxR1-negative. Pass RULE_CH3 or
    RULE_CH4 with a star rather than retyping the three numbers.
    """
    d = _stromal_with_call(cells)
    g = d.groupby("core")["TrxR1"]
    n = g.size()
    k = g.sum().astype(int)
    ok = (n >= min_stromal) & (k >= min_pos) & ((n - k) >= min_neg)
    return set(n.index[ok])


def stromal_trxr1_fraction(cells, min_cells=MIN_STROMAL_CELLS):
    """The paper's exposure: the TrxR1-positive fraction of a patient's stromal cells.

    Mean of the TrxR1 call over the patient's stromal cells that carry one, with no lineage gating.
    n_stroma counts every stromal cell of the patient, including those whose TrxR1 intensity is
    missing. Patients with fewer than min_cells
    stromal cells are dropped; on both cohorts that rule removes nobody, and it is applied anyway
    so that the value matches the written definition.
    """
    d = cells
    if "tissue" in d.columns:
        d = d[d.tissue == "Stroma"]
    g = d.groupby("Case")["TrxR1"]
    t = pd.DataFrame({"n_stroma": g.size(), "stromal_trxr1": g.mean().astype(float)})
    t = t[t.n_stroma >= min_cells].reset_index()
    t["Case"] = t.Case.astype(str)
    return t[["Case", "n_stroma", "stromal_trxr1"]]
