"""
The de-identified per-patient clinical table.

WHAT: loads the six clinical columns every chapter needs.
FROM: work_root/intermediate/<cohort>/clinical.csv, written by mif/01_clinical.py.
TO:   a DataFrame with Case (string), histology, OS_weeks, event, stage, os_source.
WHY:  one loader with one schema means a chapter script cannot pick up a different event flag or
      stage coding by accident.

Columns:
  Case       patient identifier as a string, so that a numeric-looking identifier never becomes an
             integer in one table and a string in another
  histology  LUSC, LUAD or other
  OS_weeks   overall survival in weeks from diagnosis
  event      1 dead, 0 alive at last follow-up, empty when the follow-up is missing
  stage      1 to 4 for stage I to IV, the linear covariate of every survival model in the paper
  os_source  BOMI2 only: ENCIMA or Swedish_2015, the clinical source of the follow-up
"""
import pandas as pd

from . import paths

COLUMNS = ["Case", "histology", "OS_weeks", "event", "stage", "os_source"]


def load_clinical(cohort):
    """The clinical table of a cohort, Case as string, in the order given by COLUMNS."""
    path = paths.intermediate(cohort, "clinical.csv")
    if not path.exists():
        raise FileNotFoundError(
            f"{path} is missing. Run mif/01_clinical.py first.")
    d = pd.read_csv(path, dtype={"Case": str, "histology": str, "os_source": str})
    missing = [c for c in COLUMNS if c not in d.columns]
    if missing:
        raise ValueError(f"{path} is missing the columns {missing}")
    return d[COLUMNS]


# The three covariates of the sensitivity models in Supplementary Table 1a (tumor_vs_stroma.py s9). They are not
# in clinical.csv, which keeps the six columns every chapter uses, so they are read here from the private
# clinical_full.csv by name, and nothing else of that table leaves this function.
COVARIATE_SOURCE = {"BOMI2": {"age": "Age", "sex": "Gender", "ps": "WHO_PS"},
                    "BOMI1": {"age": "Age", "sex": "Gender", "ps": "WHO"}}
SEX_CODE = {"Male": 1.0, "Female": 0.0}


def load_covariates(cohort):
    """Age (years), sex (male 1, female 0) and WHO performance status (0 to 4) of each patient.

    Read from clinical_full.csv, which mif/01_clinical.py writes beside clinical.csv. A value that is
    missing or not numeric stays missing, so that a model drops the patient rather than guessing; an
    unexpected sex code stops the run.
    """
    path = paths.intermediate(cohort, "clinical_full.csv")
    if not path.exists():
        raise FileNotFoundError(f"{path} is missing. Run mif/01_clinical.py first.")
    src = COVARIATE_SOURCE[cohort]
    f = pd.read_csv(path, dtype={"Case": str}, usecols=["Case"] + list(src.values()), low_memory=False)
    sex = f[src["sex"]].map(SEX_CODE)
    bad = f[src["sex"]].notna() & sex.isna()
    if bad.any():
        raise ValueError(f"{path}: unexpected sex codes {sorted(f.loc[bad, src['sex']].astype(str).unique())}")
    return pd.DataFrame({"Case": f.Case, "age": pd.to_numeric(f[src["age"]], errors="coerce"),
                         "male": sex, "ps": pd.to_numeric(f[src["ps"]], errors="coerce")})
