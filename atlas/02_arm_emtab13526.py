#!/usr/bin/env python3
"""
Atlas stage 02: the E-MTAB-13526 arm (De Zuani et al.) and its merge into the panel object.

WHAT: runs, in order, the steps below (build scripts in atlas/steps/, paths from config.json).
      Run time: about one hour (author obs over HTTP about 5 minutes).
      Paper items: part of the atlas behind Figure 2, Figure 3, Supplementary Figures 2 to 4 and
      Supplementary Data 1.
FROM: raw_root/external_scRNA/E-MTAB-13526 (81 Matrix Market triplets and the SDRF), HTTP range reads of
      the two author h5ad files on the EBI FTP server (97 GiB, never downloaded; needs network), and the
      stage 01 panel object.
TO:   work_root/atlas/output/RNA/meta_dataset/: emtab13526_pass1.csv; stage2/{author_obs*.parquet,
      sdrf_libraries.csv, fingerprint/, fingerprint_scores.csv, fingerprint_assignment.csv, crosswalk.csv,
      crosswalk_conflicts.csv, compartment_map_emtab13526.csv, arm_emtab13526.h5ad and tables,
      stage2_blockqc_*}; stage1a/panel_object.h5ad replaced by the merged object (2,043,367 cells,
      21 blocks; the three-arm object is kept as stage1a/panel_object_prestage2.h5ad);
      scratch/e13526/ (the per-library extracts).
WHY:  adds the E-MTAB-13526 arm to the panel object. The author "to remove" flag of this arm is
      carried in the arm but not into the full-gene object (64 cells).

Steps:
  pass1              emtab13526_pass1.py (UMI floors per library; reads 9.9 GB)
  author_obs         stage2_emtab13526_01_read_author_obs.py (imports rangefile.py)
  sdrf_libraries     stage2_emtab13526_02_sdrf_libraries.py
  fingerprint        stage2_emtab13526_06_fingerprint.py
  match              stage2_emtab13526_08_match2.py (Hungarian assignment)
  crosswalk          stage2_emtab13526_09_final_crosswalk.py (imports stage2_emtab13526_05_routeB.py)
  compartment_map    stage2_emtab13526_11_compartment_map.py (61 labels, explicit)
  extract            stage2_emtab13526_13_extract.py (4 worker processes; asserts 33,538 features)
  build_arm          stage2_emtab13526_14_build_arm.py
  verify_arm         stage2_emtab13526_15_verify_arm.py (reads the obs schema from atlas/data/)
  merge_add          stage2_merge_01_add_emtab13526.py
  merge_blockqc      stage2_merge_02_blockqc.py (block QC of the merged draft)
  merge_verify       stage2_merge_03_verify.py (moves the draft into stage1a/panel_object.h5ad)

RUN: python atlas/02_arm_emtab13526.py [--list] [--dry-run] [--only step ...] [--from step]
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import run_steps  # noqa: E402

STEPS = [
    ("pass1", ("emtab13526_pass1.py",), None),
    ("author_obs", ("stage2_emtab13526_01_read_author_obs.py",), None),
    ("sdrf_libraries", ("stage2_emtab13526_02_sdrf_libraries.py",), None),
    ("fingerprint", ("stage2_emtab13526_06_fingerprint.py",), None),
    ("match", ("stage2_emtab13526_08_match2.py",), None),
    ("crosswalk", ("stage2_emtab13526_09_final_crosswalk.py",), None),
    ("compartment_map", ("stage2_emtab13526_11_compartment_map.py",), None),
    ("extract", ("stage2_emtab13526_13_extract.py", "{SCRATCH}/e13526", "4"), None),
    ("build_arm", ("stage2_emtab13526_14_build_arm.py", "{SCRATCH}/e13526"), None),
    ("verify_arm", ("stage2_emtab13526_15_verify_arm.py",), None),
    ("merge_add", ("stage2_merge_01_add_emtab13526.py",), None),
    ("merge_blockqc", ("stage2_merge_02_blockqc.py",), None),
    ("merge_verify", ("stage2_merge_03_verify.py",), None),
]

if __name__ == "__main__":
    sys.exit(run_steps("02", STEPS))
