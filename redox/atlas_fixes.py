"""
Two corrections to the single-cell atlas object, applied wherever the object is read.

WHAT: duplicate(), the cells to leave out, and relabel(), the corrected study label of identifiers.
FROM: identifier arrays read from the obs table of the atlas object (patient_id, study, donor_id, block_id).
TO:   boolean masks and relabeled arrays for chapter2/02_working_sets.py, 07_smoking_atlas.py, 09_txnrd1_profile.py
      and atlas/05_obs_summary.py. The object itself is not changed.
WHY:  two source problems found when the sources of the atlas were checked against the published studies.

Duplicate patients. Zhang et al. 2022 (PMID 35027529) distributed, as their "batch 3", eight patients of Lambrechts
et al. 2018 (PMID 29988129). In the object they are patient_id Zhang_Chen_2022_DL1 to DL8 (technical_replicate_group
batch3, 50,651 cells), and their library names (BT1290 and so on) and cell barcodes are those of the E-MTAB-6149 and
E-MTAB-6653 deposits, which the object also holds as Lambrechts_Thienpont_2018_P1 to P8. The redistributed copies are
left out of every analysis, so that no cell enters twice; the original deposits are kept.

Study label. The study of GEO GSE325414 was labeled Wennerberg_2025 when the object was built. Its publication is
Jimenez et al., Clin Cancer Res 2026 (PMID 42307634), so the label is Jimenez_Bruno_2026 (first and last author, year,
as for the other studies). The build step (atlas/steps/wave1_arm_gse325414.py) writes the corrected label; objects
built earlier are relabeled here when read, in the study, donor, patient and block identifiers.
"""
import numpy as np

DUPLICATE_PATIENT_PREFIX = "Zhang_Chen_2022_DL"
OLD_LABEL, NEW_LABEL = "Wennerberg_2025", "Jimenez_Bruno_2026"


def duplicate(patient_ids):
    """True for the cells of the redistributed copies of the Lambrechts 2018 patients."""
    return np.char.startswith(np.asarray(patient_ids).astype(str), DUPLICATE_PATIENT_PREFIX)


def relabel(values):
    """Identifiers with the old study label at their start replaced by the corrected one."""
    v = np.asarray(values).astype(str)
    hit = np.char.startswith(v, OLD_LABEL)
    if not hit.any():
        return v
    # np.where widens the string type to the longer label; assigning into v would truncate it
    return np.where(hit, np.char.replace(v, OLD_LABEL, NEW_LABEL, count=1), v)
