#!/usr/bin/env python
"""
Atlas stage 02, step sdrf_libraries: build the authoritative SDRF library table for E-MTAB-13526.

Paper items: part of the atlas behind Figure 2, Figure 3, Supplementary Figures 2 to 4 and Supplementary Data 1.
Run: python atlas/steps/stage2_emtab13526_02_sdrf_libraries.py   (normally through atlas/02_arm_emtab13526.py).
No random numbers are drawn.

The SDRF is the authority for histology and for the sorted fraction. Library
name suffixes (_B1/_B2/_T1/_T2) do NOT encode the fraction, so nothing here is
inferred from a string suffix.

Writes: output/RNA/meta_dataset/stage2/sdrf_libraries.csv
"""
import os
import pandas as pd

import os as _os  # noqa: E402
import sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # atlas/
from common import ROOT, SCRNA, SCRATCH, DATA, RSCRIPT  # noqa: E402,F401  paths from config.json
SDRF = os.path.join(SCRNA, "E-MTAB-13526",
                    "metadata/E-MTAB-13526.sdrf.txt")
OUT = os.path.join(ROOT, "output/RNA/meta_dataset/stage2")
os.makedirs(OUT, exist_ok=True)

HIST = {
    "lung adenocarcinoma": "LUAD",
    "lung squamous cell carcinoma": "LUSC",
    "non-small cell lung cancer": "NSCLC_unspecified",
    "normal": "normal",
}
TISSUE = {
    "tumor": "primary_tumour",
    "normal tissue adjacent to tumor": "adjacent_normal",
    "upper lobe": "normal_lung",
    "lower lobe": "normal_lung",
}
ENRICH = {
    "CD45+": "CD45pos",
    "CD235a-": "CD235a_neg",
    "CD45+ + MDSC": "CD45pos_MDSC",
    "MDSC": "MDSC",
}


def one(g, col):
    """Return the single value of column col within one library (error if it varies)."""
    v = sorted(set(str(x) for x in g[col].dropna()))
    if len(v) != 1:
        raise RuntimeError(f"{col} not constant within library: {v}")
    return v[0]


def main():
    """Build one row per library from the SDRF and write sdrf_libraries.csv."""
    s = pd.read_csv(SDRF, sep="\t", dtype=str, low_memory=False)
    s.columns = [c.strip() for c in s.columns]
    s["assay_prefix"] = s["Assay Name"].str.replace(r"_S\d+$", "", regex=True)

    rows = []
    for lib, g in s.groupby("Source Name"):
        rows.append({
            "library_id": lib,
            "individual": one(g, "Characteristics[individual]"),
            "original_source_name": one(g, "Characteristics[original source name]"),
            "disease_sdrf": one(g, "Characteristics[disease]"),
            "facs_sdrf": one(g, "Characteristics[FACS]"),
            "sampling_site_sdrf": one(g, "Characteristics[sampling site]"),
            "sex_sdrf": one(g, "Characteristics[sex]"),
            "age_sdrf": one(g, "Characteristics[age]"),
            "grading_sdrf": "|".join(sorted(set(g["Characteristics[tumor grading]"].dropna()))),
            "library_construction_sdrf": one(g, "Comment[library construction]"),
            "ena_sample": one(g, "Comment[ENA_SAMPLE]"),
            "assay_prefixes": "|".join(sorted(set(g["assay_prefix"]))),
            "n_assay_prefixes": g["assay_prefix"].nunique(),
            "n_ena_runs": g["Comment[ENA_RUN]"].nunique(),
        })
    lib = pd.DataFrame(rows).sort_values("library_id").reset_index(drop=True)

    lib["histology"] = lib["disease_sdrf"].map(HIST)
    lib["tissue_class"] = lib["sampling_site_sdrf"].map(TISSUE)
    lib["enrichment"] = lib["facs_sdrf"].map(ENRICH)
    for c in ("histology", "tissue_class", "enrichment"):
        if lib[c].isna().any():
            raise RuntimeError(f"unmapped {c}: {lib.loc[lib[c].isna()].to_dict('records')}")

    # patient-level histology must be single-valued
    bad = lib.groupby("individual")["histology"].nunique()
    if (bad > 1).any():
        raise RuntimeError(f"histology varies within patient: {bad[bad>1]}")

    # assay prefixes must be unique to one library (the join key's precondition)
    exploded = lib.assign(p=lib["assay_prefixes"].str.split("|")).explode("p")
    dup = exploded.groupby("p")["library_id"].nunique()
    print("assay prefixes total:", len(dup), " shared by >1 library:", int((dup > 1).sum()))
    if (dup > 1).any():
        raise RuntimeError(f"assay prefix shared across libraries: {dup[dup>1]}")

    lib.to_csv(os.path.join(OUT, "sdrf_libraries.csv"), index=False)

    print("\nlibraries:", len(lib))
    print("\nhistology x enrichment (libraries)")
    print(pd.crosstab(lib["histology"], lib["enrichment"], margins=True).to_string())
    print("\nhistology x tissue_class (libraries)")
    print(pd.crosstab(lib["histology"], lib["tissue_class"], margins=True).to_string())
    print("\npatients per histology")
    print(lib.groupby("histology")["individual"].nunique().to_string())
    print("\nwrote", os.path.join(OUT, "sdrf_libraries.csv"))


if __name__ == "__main__":
    main()
