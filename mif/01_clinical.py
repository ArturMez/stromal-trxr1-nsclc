#!/usr/bin/env python3
"""
The per-patient clinical table of a cohort, in two versions.

WHAT: merges the raw clinical files of BOMI2 and BOMI1 into one table per cohort, restricts it to
      the patients the analysis keeps, and writes a private full table and a de-identified one.
      It produces no figure or table of the paper itself; every later step starts from its
      clinical.csv (overall survival, event, stage and histology behind every hazard ratio).
FROM: raw_root/Bomi2/{BOMI2_clinical_ENCIMA.xlsx, the Swedish clinical database of Uppsala II (file
      name from config.json input_files "bomi2_clinical_database"), cores to cases_ Bomi2-Redox1.xlsx,
      List of cores to exclude from Bomi2 20210602.xlsx} and raw_root/Bomi1/{clinicopathological data
      of Bomi1 ....xls, Core_to_patient_matching_Bomi1_Redox1reducedpanel.xlsx}.
TO:   work_root/intermediate/<cohort>/clinical_full.csv, which keeps every merged column and is
      never shared, and work_root/intermediate/<cohort>/clinical.csv, the six de-identified columns
      that every chapter script reads through redox.clinical.load_clinical (plus a .md sidecar for
      each).
WHY:  one script merges both cohorts, writes only under work_root, applies the stage grouping in
      one place, and keeps the identifying columns (dates, free-text comments) in a file that is
      marked private and stays out of the repository and out of every figure.

The de-identified table holds Case, histology, OS_weeks, event, stage and os_source, and nothing
else. Its rows are the patients of the analyzed cohort, which is derived here from the core key
files alone, so that the clinical step does not have to wait for the per-cell step:
  BOMI2  the 542 cores of the core map, less the 2 the map marks for exclusion and the 52 the
         annotation sheet flags, leaves 488 cores and 294 patients.
  BOMI1  every core of the core map is primary; the one patient whose clinical record says "Remove due
         to mix-up" is left out with its two cores, which leaves 338 cores and 216 patients.

RUN: python3 mif/01_clinical.py [--cohort BOMI2|BOMI1|both]
"""
import argparse
import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from redox import paths  # noqa: E402

warnings.filterwarnings("ignore")

# Raw file names, per cohort. They are the file names as they arrived from the two biobanks and
# are never changed here; the folder they sit in comes from the configuration file.
BOMI2_FILES = {
    "encima": "BOMI2_clinical_ENCIMA.xlsx",
    "swedish": paths.input_file("bomi2_clinical_database"),
    "cores": "cores to cases_ Bomi2-Redox1.xlsx",
    "exclusions": "List of cores to exclude from Bomi2 20210602.xlsx",
}
BOMI1_FILES = {
    "clinical": "clinicopathological data of Bomi1 for the Redox1 panel cases 20200713.xls",
    "cores": "Core_to_patient_matching_Bomi1_Redox1reducedpanel.xlsx",
}

# The Swedish clinical database has its header on the second row and no usable column names, so
# the columns are named by position.
SWEDISH_COLUMNS = ["Projektnr", "pTNM", "Stage_pTNM", "Histo_main", "Histo_sub", "HistoWHO2015",
                   "Grade", "Surgery_date", "Age", "Gender", "Vital_date", "Status", "Smoking",
                   "WHO_PS", "AdjTreatment", "AdjTreatment_type", "Comments"]
SWEDISH_EXTRA = ["pTNM", "Stage_pTNM", "Histo_main", "Histo_sub", "HistoWHO2015", "Grade",
                 "AdjTreatment", "AdjTreatment_type", "Surgery_date", "Vital_date", "Status",
                 "SE_OS_weeks_2015", "Comments"]

# Code books. The Swedish database and the BOMI1 file both code stage on the seven-step scale
# IA, IB, IIA, IIB, IIIA, IIIB, IV; every survival model in the paper uses the four-step scale.
STAGE_SEVEN_TO_FOUR = {1: 1, 2: 1, 3: 2, 4: 2, 5: 3, 6: 3, 7: 4}
GENDER = {1: "Male", 2: "Female"}
SMOKING = {1: "Yes", 2: "Yes", 3: "No"}          # current and former smoker to Yes, never to No
EVENT_TEXT = {0: "Alive", 1: "Dead"}             # Swedish status 2 and missing stay missing
HISTOLOGY = {1: ("Lung squamous cell carcinoma", "LUSC"), 2: ("Lung adenocarcinoma", "LUAD")}

# Column names of the merged BOMI2 table (the names of the ENCIMA table, which is the base).
OS_COLUMN = "Time_Diagnosis_Last_followup (OS)"
FRONT = ["ID or PAD_year", "in_BOMI2", "OS_source", "Tumor_type", "Tumor_type_code", OS_COLUMN,
         "Event_last_followup", "Age", "Gender", "clinicl_stge", "WHO_PS", "Smoking"]

# The core name of a BOMI2 image, as "glass,row,letter". The annotation sheet names the same core
# in that form already, which is what makes the two files joinable.
BOMI2_CORE = re.compile(r"BOMI2_Redox_(\d+)_Core\[(\d+),(\d+),([A-Za-z])\]")

CLINICAL_DOC = """# clinical.csv ({cohort})

One row per patient of the analyzed cohort, de-identified: nothing here but the six columns below.
{n_rows} patients. Written by mif/01_clinical.py from the raw clinical files.

| column | meaning |
|---|---|
| Case | patient identifier, as a string |
| histology | LUSC, LUAD or other |
| OS_weeks | overall survival in weeks from diagnosis |
| event | 1 dead, 0 alive at last follow-up, empty when the follow-up is missing |
| stage | 1 to 4 for stage I to IV |
| os_source | BOMI2: ENCIMA or Swedish_2015, the clinical source of the follow-up; empty for BOMI1 |

{note}
"""

FULL_DOC = """# clinical_full.csv ({cohort})

PRIVATE. One row per patient of the merged clinical table, {n_rows} rows and {n_cols} columns,
written by mif/01_clinical.py. It keeps every column of the raw clinical files, including
{private}. It is never shared and never committed. The chapter scripts read clinical.csv instead,
with two exceptions: tumor_vs_stroma.py s1 reads its histology code and smoking columns, and
redox.clinical.load_covariates reads age, sex and WHO performance status for the covariate block of
Supplementary Table 1a (tumor_vs_stroma.py s9). Nothing else of it reaches an output.

{note}
"""


# ---------------------------------------------------------------------------- BOMI2 core key files

def bomi2_core_table():
    """Core, patient, metastasis flag and exclusion flag for all 542 BOMI2 cores.

    The two key files are the same ones the per-cell step reads: the core map gives the patient of
    each core and marks two cores for exclusion, and the annotation sheet gives the tissue source
    (T1 and T2 primary, M1 and M2 metastasis) and flags 52 cores to exclude.
    """
    folder = paths.raw("Bomi2")

    def key(name):
        """The "glass,row,letter" key of a core name, or none when the name does not parse."""
        m = BOMI2_CORE.search(str(name))
        return f"{m.group(1)},{m.group(3)},{m.group(4).upper()}" if m else None

    cores = pd.read_excel(folder / BOMI2_FILES["cores"])
    cores["key"] = cores["core"].map(key)
    cores = cores.dropna(subset=["key"]).drop_duplicates("key")
    is_exclude = cores["Patient"].astype(str) == "exclude"

    sheet = pd.read_excel(folder / BOMI2_FILES["exclusions"], sheet_name="All annotations")
    sheet = sheet.rename(columns={sheet.columns[3]: "TO_EXCLUDE"})
    sheet["key"] = sheet["TMA"].astype(str).str.replace(" ", "", regex=False).str.upper()
    flagged = set(sheet.loc[sheet["TO_EXCLUDE"].astype(str).str.upper() == "YES", "key"])
    source = sheet.drop_duplicates("key").set_index("key")["Source"]

    t = pd.DataFrame({
        "key": cores["key"].to_numpy(),
        "Case": np.where(is_exclude, None, cores["Patient"].to_numpy()),
        "source": [source.get(k) for k in cores["key"]],
    })
    t["TM"] = t["source"].astype(str).str[0]
    t["excluded"] = is_exclude.to_numpy() | t["key"].isin(flagged).to_numpy()
    return t


def bomi2_analyzed_cases():
    """The 294 BOMI2 patients the paper analyzes, from the core key files.

    A core is kept when neither key file excludes it, and a patient enters the cohort when a kept
    core of that patient is primary tissue.

    One patient is reclassified: the annotation sheet calls both of that patient's surviving cores
    metastasis, a review of the tissue found them to be primary, and the correction takes the
    cohort from 293 to 294 patients. The rule derives the correction from the sheet instead of
    naming the cores: after the exclusions, exactly one patient is left with metastasis cores only,
    and it is that patient. This keeps core identifiers out of the code.
    """
    t = bomi2_core_table()
    kept = t[(~t.excluded) & t.Case.notna()].copy()
    by_case = kept.groupby("Case")["TM"].apply(lambda s: set(s))
    metastasis_only = sorted(by_case.index[by_case.apply(lambda s: s == {"M"})])
    if len(metastasis_only) != 1:
        raise ValueError(
            f"expected exactly one metastasis-only patient to reclassify, found "
            f"{len(metastasis_only)}. The annotation sheet has changed; check it before trusting "
            "the cohort count.")
    kept.loc[kept.Case.isin(metastasis_only), "TM"] = "T"
    cases = sorted({str(c) for c in kept.loc[kept.TM == "T", "Case"]})
    print(f"[BOMI2] {len(t)} cores, {int(t.excluded.sum())} excluded, "
          f"{int((~t.excluded).sum())} kept, {len(cases)} patients "
          f"(one metastasis-only patient reclassified as primary)", flush=True)
    return cases


# ---------------------------------------------------------------------------------- BOMI2 merge

def _first_int(x):
    """The first run of digits in a cell of the Swedish database, as an integer."""
    m = re.search(r"\d+", str(x))
    return int(m.group()) if m else np.nan


def merge_bomi2():
    """The merged BOMI2 clinical table: ENCIMA, plus the patients only the Swedish database has.

    ENCIMA covers 306 of the 356 patients and is the base. The 50 BOMI2 patients it does not cover
    are taken from the Swedish clinical database, where survival is computed from the surgery and
    vital dates and censored at the closing date of that extract, 17 June 2015. The Swedish columns
    that ENCIMA does not carry are appended with an SE_ prefix for every patient.
    """
    folder = paths.raw("Bomi2")
    cores = pd.read_excel(folder / BOMI2_FILES["cores"])
    in_cohort = set(pd.to_numeric(cores.loc[cores.Patient.astype(str) != "exclude", "Patient"],
                                  errors="coerce").dropna().astype(int))

    en = pd.read_excel(folder / BOMI2_FILES["encima"])
    empty = [c for c in en.columns if en[c].notna().sum() == 0]
    en = en.drop(columns=empty)
    en = en.drop(columns=[c for c in ["RFS", "RFS event", "TTR", "TTR event"] if c in en.columns])
    en["pid"] = pd.to_numeric(en["ID or PAD_year"], errors="coerce")
    en["OS_source"] = "ENCIMA"
    print(f"[BOMI2] ENCIMA {len(en)} patients, {en.shape[1] - 2} columns, "
          f"{len(empty)} empty columns dropped", flush=True)

    sw = pd.read_excel(folder / BOMI2_FILES["swedish"], sheet_name="Kliniska data", header=1)
    sw = sw.rename(columns={sw.columns[i]: n
                            for i, n in enumerate(SWEDISH_COLUMNS[:len(sw.columns)])})
    sw["pid"] = sw["Projektnr"].apply(_first_int)
    sw = sw.dropna(subset=["pid"])
    sw["pid"] = sw["pid"].astype(int)
    sw = sw.drop_duplicates("pid", keep="first")
    surgery = pd.to_datetime(sw["Surgery_date"], errors="coerce")
    vital = pd.to_datetime(sw["Vital_date"], errors="coerce")
    sw["SE_OS_weeks_2015"] = (vital - surgery).dt.days / 7.0

    missing = sorted(in_cohort - set(en["pid"].dropna().astype(int)))
    indexed = sw.set_index("pid")
    rows = []
    for pid in missing:
        if pid not in indexed.index:
            continue
        r = indexed.loc[pid]
        code = _first_int(r["Histo_main"]) if pd.notna(r["Histo_main"]) else None
        name, short = HISTOLOGY.get(code, (np.nan, np.nan))
        stage = pd.to_numeric(r["Stage_pTNM"], errors="coerce")
        rows.append({
            "ID or PAD_year": pid, "pid": pid, "OS_source": "Swedish_2015",
            "TMA_cohort": np.nan, "Tumor_type": name, "Tumor_type_code": short,
            OS_COLUMN: r["SE_OS_weeks_2015"],
            "Event_last_followup": EVENT_TEXT.get(r["Status"] if pd.notna(r["Status"]) else None,
                                                  np.nan),
            "Age": r["Age"],
            "Gender": GENDER.get(_first_int(r["Gender"]) if pd.notna(r["Gender"]) else None,
                                 np.nan),
            "clinicl_stge": STAGE_SEVEN_TO_FOUR.get(int(stage) if pd.notna(stage) else None,
                                                    np.nan),
            "WHO_PS": r["WHO_PS"],
            "Smoking": SMOKING.get(_first_int(r["Smoking"]) if pd.notna(r["Smoking"]) else None,
                                   np.nan),
        })
    print(f"[BOMI2] {len(rows)} of {len(missing)} patients outside ENCIMA taken from the Swedish "
          "database", flush=True)
    merged = pd.concat([en, pd.DataFrame(rows)], ignore_index=True)

    extra = sw.set_index("pid")[SWEDISH_EXTRA].add_prefix("SE_")
    extra.columns = [c.replace("SE_SE_", "SE_") for c in extra.columns]
    merged = merged.merge(extra, left_on="pid", right_index=True, how="left")
    merged["in_BOMI2"] = merged["pid"].isin(in_cohort)
    front = [c for c in FRONT if c in merged.columns]
    merged = merged[front + [c for c in merged.columns if c not in front and c != "pid"]]
    merged.insert(0, "Case", merged["ID or PAD_year"].astype(str))
    return merged


# ---------------------------------------------------------------------------------- BOMI1 merge

def bomi1_analyzed_cases():
    """The BOMI1 patients, from the core map: every core is primary. A patient whose raw clinical
    record says to remove them is left out (one patient, "Remove due to mix-up" in Comments), as
    mif/02_cells.R leaves out their cores. The rule reads the comment, so no identifier is written here."""
    cc = pd.read_excel(paths.raw("Bomi1") / BOMI1_FILES["cores"])
    cc.columns = [str(c).strip() for c in cc.columns]
    clin = pd.read_excel(paths.raw("Bomi1") / BOMI1_FILES["clinical"])
    removed = {str(int(c)) for c in clin.loc[clin["Comments"].astype(str).str.strip().str.lower()
                                             .str.startswith("remove"), "Project ID"]}
    cases = sorted({str(int(c)) for c in pd.to_numeric(cc["Case"], errors="coerce").dropna()} - removed)
    print(f"[BOMI1] {len(cc)} cores, {len(removed)} patient(s) removed by the clinical record, "
          f"{len(cases)} patients", flush=True)
    return cases


def merge_bomi1():
    """The BOMI1 clinical table, renamed into the BOMI2 column schema.

    The raw file codes histology as 1 squamous, 2 adenocarcinoma and 3, 4 and 9 as other kinds that
    the paper groups together; survival is in days and is converted to weeks so that the two
    cohorts share one unit; stage stays on its seven-step scale here and is grouped further down.
    """
    cl = pd.read_excel(paths.raw("Bomi1") / BOMI1_FILES["clinical"])
    cl.columns = [str(c).strip() for c in cl.columns]
    days = pd.to_numeric(cl["Time to death (360 dagar)"], errors="coerce")
    out = pd.DataFrame({
        "Case": cl["Project ID"].astype("Int64").astype(str),
        "Tumor_type_code": pd.to_numeric(cl["Histology 2011-4 PMI"], errors="coerce")
                             .map({k: v[1] for k, v in HISTOLOGY.items()}),
        OS_COLUMN: days / 7.0,
        "OS_days_raw": days,
        "Event_last_followup": pd.to_numeric(cl["Dead"], errors="coerce").map(EVENT_TEXT),
        "clinicl_stge": pd.to_numeric(cl["Stage"], errors="coerce"),
        "Age": pd.to_numeric(cl["Age"], errors="coerce"),
        "Gender": pd.to_numeric(cl["Gender"], errors="coerce").map(GENDER),
    })
    keep = [c for c in cl.columns if c not in
            ("Project ID", "Histology 2011-4 PMI", "Time to death (360 dagar)", "Dead", "Stage",
             "Age", "Gender")]
    out = pd.concat([out, cl[keep].reset_index(drop=True)], axis=1)
    print(f"[BOMI1] clinical file {len(out)} patients", flush=True)
    return out


# ------------------------------------------------------------------------------- de-identifying

def deidentify(full, cohort, cases):
    """The six shared columns, for the patients of the analyzed cohort, sorted by Case.

    histology  anything that is neither LUSC nor LUAD, including a missing code, becomes "other",
               which is the group the paper reports and the group every histology model excludes.
    stage      BOMI2 is already on the four-step scale; BOMI1 is grouped from the seven-step one.
    event      1 for a death, 0 for a patient alive at last contact, empty when the raw file has
               neither, which happens for one patient in each cohort.
    """
    d = full.set_index(full["Case"].astype(str))
    missing = [c for c in cases if c not in d.index]
    if missing:
        raise ValueError(f"{cohort}: {len(missing)} analyzed patients are not in the clinical "
                         f"table, the first of them at position {cases.index(missing[0])}")
    d = d.loc[cases]
    stage = pd.to_numeric(d["clinicl_stge"], errors="coerce")
    if cohort == "BOMI1":
        stage = stage.map(STAGE_SEVEN_TO_FOUR)
        source = pd.Series([pd.NA] * len(d), index=d.index, dtype="object")
    else:
        if stage.dropna().max() > 4:
            raise ValueError("BOMI2 stage is expected on the four-step scale already")
        source = d["OS_source"].astype("object")
    histology = d["Tumor_type_code"].where(d["Tumor_type_code"].isin(["LUSC", "LUAD"]), "other")
    # The two integer columns stay on the nullable Int8 type, so that a patient without a follow-up
    # or without a stage is written as an empty field and every other patient as a plain integer.
    # A float column would write 1.0 where the table means one death.
    out = pd.DataFrame({
        "histology": histology,
        "OS_weeks": pd.to_numeric(d[OS_COLUMN], errors="coerce"),
        "event": d["Event_last_followup"].map({"Dead": 1, "Alive": 0}).astype("Int8"),
        "stage": stage.astype("Int8"),
        "os_source": source,
    })
    out.insert(0, "Case", list(d.index))
    return out.reset_index(drop=True)


NOTE = {
    "BOMI2": ("BOMI2 stage is the clinical stage of the ENCIMA table for the patients it covers "
              "and the pathological pTNM of the Swedish database for the rest; the Swedish "
              "patients are censored at 17 June 2015, which is the closing date of that extract. "
              "os_source says which of the two a patient's follow-up came from."),
    "BOMI1": ("BOMI1 stage is grouped from the seven-step scale IA to IV by the fixed map "
              "{1,2} to 1, {3,4} to 2, {5,6} to 3, 7 to 4, which is the map every published BOMI1 "
              "survival model used. The cohort has one clinical source, so os_source is empty."),
}
PRIVATE = {
    "BOMI2": ("the surgery and vital dates of every patient, the pTNM string, the free-text "
              "comments of the Swedish database and its granular histology"),
    "BOMI1": ("the free-text comments, the treatment and recurrence columns and the block and row "
              "of the tissue array"),
}


def run(cohort):
    """Merge, de-identify and write both tables (and their .md sidecars) of one cohort."""
    if cohort == "BOMI2":
        full, cases = merge_bomi2(), bomi2_analyzed_cases()
    else:
        full, cases = merge_bomi1(), bomi1_analyzed_cases()

    target = paths.intermediate(cohort, "clinical_full.csv")
    full.to_csv(target, index=False)
    (target.parent / "clinical_full.md").write_text(FULL_DOC.format(
        cohort=cohort, n_rows=len(full), n_cols=full.shape[1], private=PRIVATE[cohort],
        note=NOTE[cohort]))
    print(f"[{cohort}] wrote {target.name} ({len(full)} rows, {full.shape[1]} columns, PRIVATE)",
          flush=True)

    out = deidentify(full, cohort, cases)
    target = paths.intermediate(cohort, "clinical.csv")
    out.to_csv(target, index=False)
    (target.parent / "clinical.md").write_text(CLINICAL_DOC.format(
        cohort=cohort, n_rows=len(out), note=NOTE[cohort]))
    counts = out.histology.value_counts().to_dict()
    print(f"[{cohort}] wrote {target.name} ({len(out)} patients, {counts}, "
          f"{int(out.event.isna().sum())} without follow-up, "
          f"{int(out.stage.isna().sum())} without stage)", flush=True)


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
