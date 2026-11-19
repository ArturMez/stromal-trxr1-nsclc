#!/usr/bin/env python3
"""
Atlas stage 04: the full-transcriptome arms, fullgene_object.h5ad, panel_object_v2.h5ad, the block
QC on 24 blocks, its transfer, and the obs repair.

WHAT: runs, in order, the steps below (build scripts in atlas/steps/, paths from config.json).
      Run time: several hours; the LuCA and E-MTAB-13526 extractions dominate. Disk: about 24 GB of
      full arms and 25 GB for the object.
      Paper items: fullgene_object.h5ad is the atlas that chapter2/ reads for Figure 2, Figure 3,
      Supplementary Figures 2 to 4 and Supplementary Data 1.
FROM: the raw deposits, the stage 01 to 03 arms and stage1a/panel_object.h5ad (harmonized obs of the
      four original arms).
TO:   work_root/atlas/output/RNA/meta_dataset/fullbuild/{arm_*_full.h5ad, e6149/,
      fullgene_object.h5ad (1,888,910 cells x 50,552 genes), fullgene_gene_presence.csv,
      verify_report.csv, panel_object_v2.h5ad (2,500,510 x 1,721, 24 blocks), blockqc_v2/*}.
WHY:  builds the full-gene atlas object. Two obs repairs run in the stage that creates the defects,
      right after the steps they repair:
        transfer_qc  fullbuild_06_transfer_qc.py writes block_qc_pass and block_qc_assessed onto
                     the full-gene object and sets n_genes to -1 (with n_genes_is_full) on
                     panel_object_v2 where no full-transcriptome count exists
        repair_obs   fullbuild_10_repair_obs.py retypes n_counts_panel, n_genes_panel (float, from
                     panel_object_v2 by cell id), malignancy_score (all missing) and is_malignant
                     (nullable, all missing); its backup goes to scratch/repair_obs_backup
      Other obs and var columns present in earlier objects are not rebuilt; none of them is read by
      chapter2/ or by the obs summary of stage 05.

Steps:
  extract_luca, extract_emtab13526, extract_emtab6653,
  dump6149_rows (R), export6149_full (R, the stage 01 exporter with the full gene list),
  extract_emtab6149, extract_wave1, merge, verify (SEED 31), merge_panel,
  blockqc_v2 (stage2_merge_02_blockqc.py with BLOCKQC_DRAFT and BLOCKQC_OUT on panel_object_v2),
  transfer_qc, repair_obs, state (fullbuild_state.py against the registry in atlas/data/)

RUN: python atlas/04_build_objects.py [--list] [--dry-run] [--only step ...] [--from step]
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import run_steps  # noqa: E402

STEPS = [
    ("extract_luca", ("fullbuild_01_extract_luca.py",), None),
    ("extract_emtab13526", ("fullbuild_01b_extract_emtab13526.py",), None),
    ("extract_emtab6653", ("fullbuild_01c_extract_emtab6653.py",), None),
    ("dump6149_rows", ("fullbuild_01d_dump6149_rows.R", "{MD}/fullbuild/e6149", "{SCRNA}/E-MTAB-6149/processed"), None),
    ("export6149_full", ("stage1a_export6149.R", "{MD}/fullbuild/e6149", "{MD}/fullbuild/e6149/all_gene_ids.txt",
                         "{SCRNA}/E-MTAB-6149/processed", "{MD}/stage1a/panel_genes.csv"), None),
    ("extract_emtab6149", ("fullbuild_01d_extract_emtab6149.py",), None),
    ("extract_wave1", ("fullbuild_01e_extract_wave1.py", "all"), None),
    ("merge", ("fullbuild_02_merge.py",), None),
    ("verify", ("fullbuild_04_verify.py",), None),
    ("merge_panel", ("fullbuild_05_merge_panel.py",), None),
    ("blockqc_v2", ("stage2_merge_02_blockqc.py",),
     {"BLOCKQC_DRAFT": "{MD}/fullbuild/panel_object_v2.h5ad", "BLOCKQC_OUT": "{MD}/fullbuild/blockqc_v2/"}),
    ("transfer_qc", ("fullbuild_06_transfer_qc.py",), None),
    ("repair_obs", ("fullbuild_10_repair_obs.py", "--backup-dir", "{SCRATCH}/repair_obs_backup"), None),
    ("state", ("fullbuild_state.py",), None),
]

if __name__ == "__main__":
    sys.exit(run_steps("04", STEPS))
