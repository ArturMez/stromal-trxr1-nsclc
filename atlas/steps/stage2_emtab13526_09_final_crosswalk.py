#!/usr/bin/env python
"""
Atlas stage 02, step crosswalk: final, verified crosswalk for E-MTAB-13526.

Paper items: part of the atlas behind Figure 2, Figure 3, Supplementary Figures 2 to 4 and Supplementary Data 1.
Run: python atlas/steps/stage2_emtab13526_09_final_crosswalk.py   (normally through atlas/02_arm_emtab13526.py).
No random numbers are drawn; every check (F1 to F3, G4 to G9) is printed as PASS or FAIL.

AUTHORITY FOR library_id: the (barcode, UMI) fingerprint of the match step. That route
reads the actual deposited count matrices and is independent of every naming
convention, of the author's free-text sangerID, and of the SDRF's own site
labels.

AUTHORITY FOR histology and sorted fraction: the SDRF.

Three routes are compared:
  A  ordinal name   'Patient 15 Tumor 2' -> P15_T2
  B  sequencing id  sangerID -> SDRF FASTQ file names -> library
  F  fingerprint    author (barcode, total_counts) -> raw library   [AUTHORITY]

Writes:
  output/RNA/meta_dataset/stage2/crosswalk.csv
  output/RNA/meta_dataset/stage2/crosswalk_conflicts.csv
"""
import os
import re
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stage2_emtab13526_05_routeB import build_route_b, parse_patient_sample  # noqa: E402

import os as _os  # noqa: E402
import sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # atlas/
from common import ROOT, SCRNA, SCRATCH, DATA, RSCRIPT  # noqa: E402,F401  paths from config.json
S2 = os.path.join(ROOT, "output/RNA/meta_dataset/stage2")
MD = os.path.join(ROOT, "output/RNA/meta_dataset")

pd.set_option("display.width", 270)
pd.set_option("display.max_rows", 400)
pd.set_option("display.max_colwidth", 70)

FAIL = []


def gate(name, ok, detail=""):
    """Print a PASS or FAIL line for one crosswalk check and record failures."""
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}  {detail}", flush=True)
    if not ok:
        FAIL.append(name)


SORTING_TO_FACS = {
    "Myeloid": "CD45+", "Myeloid ": "CD45+", "CD235a-": "CD235a-",
    "Myeloid + eMDSC": "CD45+ + MDSC", "Myeloid+ eMDSC": "CD45+ + MDSC",
    "MDSC": "MDSC", "eMDSC": "MDSC",
}
ENV_TO_SITECLASS = {"Tumour": "primary_tumour", "Background": "adjacent_normal",
                    "Healthy": "normal_lung"}


def main():
    """Compare the three routes, apply the checks and write crosswalk.csv and crosswalk_conflicts.csv."""
    obs = pd.read_parquet(os.path.join(S2, "author_obs.parquet"),
                          columns=["author_object", "patient", "sample", "sorting",
                                   "sangerID", "patient_sample", "environment"])
    lib = pd.read_csv(os.path.join(S2, "sdrf_libraries.csv"))
    fpa = pd.read_csv(os.path.join(S2, "fingerprint_assignment.csv"))
    p1 = pd.read_csv(os.path.join(MD, "emtab13526_pass1.csv")).rename(
        columns={"library": "library_id"})

    grp = (obs.groupby("patient_sample", observed=True)
              .agg(author_object=("author_object", "first"),
                   patient=("patient", "first"),
                   sample=("sample", "first"),
                   sorting=("sorting", "first"),
                   sangerID=("sangerID", "first"),
                   environment=("environment", "first"),
                   author_n_cells=("patient", "size"))
              .reset_index())

    print("=" * 78)
    print("ROUTE COMPARISON")
    print("=" * 78)

    grp["library_id_A"] = grp["patient_sample"].map(parse_patient_sample)
    rb = build_route_b()[["patient_sample", "library_id_B"]]
    grp = grp.merge(rb, on="patient_sample", how="left")
    grp = grp.merge(fpa[["grp", "greedy_library", "greedy_score", "runner_up_score"]]
                    .rename(columns={"grp": "patient_sample",
                                     "greedy_library": "library_id_F",
                                     "greedy_score": "fingerprint_score",
                                     "runner_up_score": "fingerprint_runner_up"}),
                    on="patient_sample", how="left")

    gate("F1 fingerprint resolved all 81 groups bijectively",
         grp["library_id_F"].notna().all() and grp["library_id_F"].nunique() == 81
         and set(grp["library_id_F"]) == set(lib["library_id"]),
         f"distinct={grp['library_id_F'].nunique()}")
    gate("F2 every fingerprint score >= 0.99 and beats runner-up 5-fold",
         bool((grp["fingerprint_score"] >= 0.99).all())
         and bool((grp["fingerprint_score"] > 5 * grp["fingerprint_runner_up"]).all()),
         f"min_score={grp['fingerprint_score'].min():.4f} "
         f"max_runner_up={grp['fingerprint_runner_up'].max():.4f}")

    covB = grp["library_id_B"].notna()
    agreeBF = grp.loc[covB, "library_id_B"] == grp.loc[covB, "library_id_F"]
    gate("F3 sequencing-id route agrees with fingerprint where it resolves",
         bool(agreeBF.all()), f"covered={int(covB.sum())}/81 agree={int(agreeBF.sum())}")

    disA = grp[grp["library_id_A"] != grp["library_id_F"]]
    print(f"\n  Route A (ordinal name) disagrees with the fingerprint for "
          f"{len(disA)} of 81 groups:")
    if len(disA):
        print(disA[["patient_sample", "sangerID", "library_id_A",
                    "library_id_B", "library_id_F", "author_n_cells"]].to_string(index=False))

    grp["library_id"] = grp["library_id_F"]
    x = grp.merge(lib, on="library_id", how="left", validate="one_to_one")

    print("\n" + "=" * 78)
    print("VERIFICATION AGAINST THE SDRF")
    print("=" * 78)

    bad = x[x["patient"] != x["individual"]]
    gate("G4 author patient == SDRF individual", len(bad) == 0, f"mismatches={len(bad)}")
    if len(bad):
        print(bad[["patient_sample", "patient", "library_id", "individual"]].to_string(index=False))

    x["facs_from_author"] = x["sorting"].map(SORTING_TO_FACS)
    fbad = x[x["facs_from_author"] != x["facs_sdrf"]]
    gate("G5 author sorting == SDRF FACS",
         len(fbad) == 0,
         f"unmapped={sorted(set(x.loc[x['facs_from_author'].isna(),'sorting']))} "
         f"mismatches={len(fbad)}")
    if len(fbad):
        print(fbad[["library_id", "patient_sample", "individual", "sorting",
                    "facs_from_author", "facs_sdrf"]].to_string(index=False))

    x["site_from_author"] = x["environment"].map(ENV_TO_SITECLASS)
    sbad = x[x["site_from_author"] != x["tissue_class"]]
    gate("G6 author environment == SDRF sampling site class", len(sbad) == 0,
         f"mismatches={len(sbad)}")
    if len(sbad):
        print(sbad[["library_id", "individual", "patient_sample", "environment",
                    "site_from_author", "sampling_site_sdrf", "tissue_class",
                    "author_n_cells"]].to_string(index=False))

    x = x.merge(p1[["library_id", "bc_ge_100umi", "bc_ge_500umi", "bc_ge_1000umi",
                    "bc_ge_2000umi", "nonzero_entries", "total_umi"]],
                on="library_id", how="left", validate="one_to_one")
    over = x[x["author_n_cells"] > x["bc_ge_100umi"]]
    gate("G7 author cells <= barcodes with >=100 UMI", len(over) == 0,
         f"libraries exceeding={len(over)}")
    if len(over):
        print(over[["library_id", "patient_sample", "author_n_cells",
                    "bc_ge_100umi"]].to_string(index=False))

    ph = x.groupby("individual")["histology"].nunique()
    gate("G8 histology single-valued per patient", bool((ph == 1).all()),
         f"patients={len(ph)}")

    gate("G9 every SDRF library has author cells",
         int((x["author_n_cells"] == 0).sum()) == 0 and len(x) == 81,
         f"libraries with 0 cells={int((x['author_n_cells']==0).sum())}")

    # ---- conflict register ----------------------------------------------
    conf = []
    for _, r in x.iterrows():
        if r["site_from_author"] != r["tissue_class"]:
            conf.append({"library_id": r["library_id"], "individual": r["individual"],
                         "field": "tissue_class",
                         "author_value": r["environment"],
                         "sdrf_value": r["sampling_site_sdrf"],
                         "n_cells": r["author_n_cells"],
                         "note": "author h5ad and SDRF disagree on tumour vs adjacent normal"})
        if r["facs_from_author"] != r["facs_sdrf"]:
            conf.append({"library_id": r["library_id"], "individual": r["individual"],
                         "field": "enrichment",
                         "author_value": r["sorting"],
                         "sdrf_value": r["facs_sdrf"],
                         "n_cells": r["author_n_cells"],
                         "note": "author h5ad and SDRF disagree on the sorted fraction"})
    conf = pd.DataFrame(conf)
    conf.to_csv(os.path.join(S2, "crosswalk_conflicts.csv"), index=False)

    x["tissue_class_disputed"] = x["site_from_author"] != x["tissue_class"]
    x["enrichment_disputed"] = x["facs_from_author"] != x["facs_sdrf"]
    x["crosswalk_route_A_agrees"] = x["library_id_A"] == x["library_id"]

    print("\n" + "=" * 78)
    print("FINAL CROSSWALK")
    print("=" * 78)
    show = ["library_id", "patient_sample", "individual", "histology", "facs_sdrf",
            "enrichment", "tissue_class", "author_object", "author_n_cells",
            "bc_ge_100umi", "fingerprint_score", "crosswalk_route_A_agrees",
            "tissue_class_disputed", "enrichment_disputed"]
    print(x.sort_values("library_id")[show].to_string(index=False))

    print("\ntotal author cells:", int(x["author_n_cells"].sum()))
    print("\ncells per histology x FACS")
    print(pd.pivot_table(x, index="histology", columns="facs_sdrf",
                         values="author_n_cells", aggfunc="sum",
                         margins=True, fill_value=0).to_string())
    print("\ncells per histology x tissue_class")
    print(pd.pivot_table(x, index="histology", columns="tissue_class",
                         values="author_n_cells", aggfunc="sum",
                         margins=True, fill_value=0).to_string())

    keep = ["library_id", "patient_sample", "sangerID", "author_object", "patient",
            "sample", "sorting", "environment", "individual", "original_source_name",
            "disease_sdrf", "facs_sdrf", "sampling_site_sdrf", "histology",
            "tissue_class", "enrichment", "sex_sdrf", "age_sdrf", "grading_sdrf",
            "library_construction_sdrf", "ena_sample", "author_n_cells",
            "bc_ge_100umi", "bc_ge_500umi", "bc_ge_1000umi", "nonzero_entries",
            "total_umi", "library_id_A", "library_id_B", "library_id_F",
            "fingerprint_score", "fingerprint_runner_up",
            "crosswalk_route_A_agrees", "tissue_class_disputed", "enrichment_disputed"]
    out = x[keep].sort_values("library_id").reset_index(drop=True)
    out.to_csv(os.path.join(S2, "crosswalk.csv"), index=False)
    print("\nwrote", os.path.join(S2, "crosswalk.csv"))
    print("wrote", os.path.join(S2, "crosswalk_conflicts.csv"))
    if len(conf):
        print("\nCONFLICT REGISTER")
        print(conf.to_string(index=False))

    print("\n" + "=" * 78)
    print("GATES FAILED: " + ", ".join(FAIL) if FAIL else "ALL GATES PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
