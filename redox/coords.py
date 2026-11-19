"""
One coordinate convention for both cohorts.

WHAT: converts the raw inForm cell positions to micrometers.
FROM: the "Cell X Position" and "Cell Y Position" columns of the raw exports.
TO:   an array of micrometers, whatever the cohort.
WHY:  the two cohorts do not store positions in the same unit. Uppsala II (BOMI2) is in pixels at
      0.4964 micrometers per pixel; Uppsala I (BOMI1) is already in micrometers. Treating both as
      micrometers would make a 30 micrometer radius 14.9 micrometers in BOMI2. One constant, used
      everywhere including the tissue areas.

      Evidence for 0.4964: on the BOMI2 cores and the BOMI1 cores (the export holds 340 cores, 338
      after the exclusion of one patient with 2 cores), the median core diameter is 1.19 mm in
      both cohorts once BOMI2 is scaled, and the median nearest-neighbor distance is 6.25 against
      6.51 micrometers.
"""
import numpy as np

MPP = {"BOMI2": 0.4964, "BOMI1": 1.0}


def scale(cohort):
    """Micrometers per raw coordinate unit for a cohort."""
    if cohort not in MPP:
        raise KeyError(f"no coordinate scale recorded for {cohort}")
    return MPP[cohort]


def to_um(xy, cohort):
    """Raw coordinates to micrometers. Accepts any array-like of shape (n, 2) or (n,)."""
    return np.asarray(xy, dtype=float) * scale(cohort)
