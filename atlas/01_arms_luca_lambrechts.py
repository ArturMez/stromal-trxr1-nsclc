#!/usr/bin/env python3
"""
Atlas stage 01: the gene panel, the three original arms (LuCA, E-MTAB-6149, E-MTAB-6653), the
three-arm panel object, its stage-1 block QC and the E-MTAB-6653 label transfer.

WHAT: runs, in order, the steps below (build scripts in atlas/steps/, paths from config.json; see
      atlas/common.py). Run time: about 30 minutes for the panel and arms, 4 minutes for the block QC
      and about 45 minutes for the label transfer.
      Paper items: this stage builds the first part of the atlas behind Figure 2, Figure 3,
      Supplementary Figures 2 to 4 and Supplementary Data 1 (chapter2/ reads the final object).
FROM: raw_root/external_scRNA/{LuCA, E-MTAB-6149, E-MTAB-6653} (00_download.sh).
TO:   work_root/atlas/output/RNA/meta_dataset/: build-plan_compartment_map.csv, lambrechts.csv,
      stage1a/{_features6149.csv, panel_genes.csv, panel_genes_source_presence.csv, arm_luca.h5ad and
      its tables, arm_emtab6149.h5ad and tables, arm_emtab6653.h5ad, panel_object.h5ad (1,146,594 cells
      x 1,906 genes, 19 blocks) and its tables, blockqc_*, labeltransfer_*}; scratch/e6149/ (the R
      export and panel_resolution.csv) and scratch/labeltransfer_emtab6653/.
WHY:  builds the gene-panel object of the three original arms. Order matters: the label transfer
      asserts the 1,146,594-cell three-arm object, so it must run before stage 02 adds E-MTAB-13526.

Steps:
  build_plan           build-plan.py; without the Census obs cache (acquisition.py, not packaged) it
                       writes the compartment map only, the one output a later step reads
  lambrechts           lambrechts.py (R through config "rscript" or PATH; scratch under work_root);
                       writes lambrechts.csv
  dump6149_features    stage1a_dump6149_features.R (paths are arguments)
  panel                stage1a_panel.py (SEED 0 for the HVG cell sample)
  arm_luca             stage1a_arm_luca.py; asserts 66,296 Lambrechts cells dropped, 1,061,555 kept
  arm_luca_verify      stage1a_arm_luca_verify.py
  resolve6149          stage1a_resolve6149_panel.py: panel_resolution.csv (resolution of the panel
                       genes in E-MTAB-6149)
  export6149           stage1a_export6149.R (paths are arguments 3 and 4)
  arm_emtab6149        stage1a_arm_emtab6149.py; asserts 52,698 cells
  arm_emtab6653        stage1a_arm_emtab6653.py; reads the obs schema from atlas/data/
  merge_panel_object   stage1a_merge_panel_object.py
  verify_panel_object  stage1a_verify_panel_object.py
  blockqc_01..03       the stage-1 block QC (stage 02 repeats it on the larger object, but the
                       label transfer reads this flag); asserts 1,146,594 cells and 1,906 genes
  lt_build, lt_validate, lt_apply, lt_write
                       stage1b_label_transfer_emtab6653_v2.py: build, validate, apply and write,
                       with the variant and floor recorded in the object (perbatch, 0.5); SEED 0

RUN: python atlas/01_arms_luca_lambrechts.py [--list] [--dry-run] [--only step ...] [--from step]
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import run_steps  # noqa: E402

LT = "stage1b_label_transfer_emtab6653_v2.py"
STEPS = [
    ("build_plan", ("build-plan.py",), None),
    ("lambrechts", ("lambrechts.py",), None),
    ("dump6149_features", ("stage1a_dump6149_features.R", "{SCRNA}/E-MTAB-6149/processed",
                           "{MD}/stage1a/_features6149.csv"), None),
    ("panel", ("stage1a_panel.py",), None),
    ("arm_luca", ("stage1a_arm_luca.py",), None),
    ("arm_luca_verify", ("stage1a_arm_luca_verify.py",), None),
    ("resolve6149", ("stage1a_resolve6149_panel.py", "{SCRATCH}/e6149"), None),
    ("export6149", ("stage1a_export6149.R", "{SCRATCH}/e6149", "{SCRATCH}/e6149/panel_rows.txt",
                    "{SCRNA}/E-MTAB-6149/processed", "{MD}/stage1a/panel_genes.csv"), None),
    ("arm_emtab6149", ("stage1a_arm_emtab6149.py", "{SCRATCH}/e6149"), None),
    ("arm_emtab6653", ("stage1a_arm_emtab6653.py",), None),
    ("merge_panel_object", ("stage1a_merge_panel_object.py",), None),
    ("verify_panel_object", ("stage1a_verify_panel_object.py",), None),
    ("blockqc_01", ("blockqc_01_detection_by_block.py",), None),
    ("blockqc_02", ("blockqc_02_rule_and_flag.py",), None),
    ("blockqc_03", ("blockqc_03_write_flag.py",), None),
    ("lt_build", (LT, "build"), None),
    ("lt_validate", (LT, "validate"), None),
    ("lt_apply", (LT, "apply", "--variant", "perbatch", "--floor", "0.5"), None),
    ("lt_write", (LT, "write", "--variant", "perbatch", "--floor", "0.5"), None),
]

if __name__ == "__main__":
    sys.exit(run_steps("01", STEPS))
