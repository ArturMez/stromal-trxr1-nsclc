#!/usr/bin/env python3
"""
Atlas stage 03: the three wave-1 arms (STTT, GSE325414, GSE246215) at the panel gene axis.

WHAT: runs, in order, the steps below (build scripts in atlas/steps/, paths from config.json).
      Run time: GSE246215 and GSE325414 take about a minute each, STTT about 80 minutes (a dense
      10 GB csv streamed line by line).
      Paper items: part of the atlas behind Figure 2, Figure 3, Supplementary Figures 2 to 4 and
      Supplementary Data 1.
FROM: raw_root/external_scRNA/{GSE325414, GSE246215, STTT_Zhang2022} and the panel gene list of
      stage1a/panel_object.h5ad (stage 02).
TO:   work_root/atlas/output/RNA/meta_dataset/wave1/{gse325414_donor_meta.csv, arm_gse325414.h5ad,
      arm_gse246215.h5ad, arm_sttt.h5ad}; scratch/recon/{g325414_rowmap.json, g325414_panel.tsv.gz}.
WHY:  brings three more public studies to the panel gene axis. The GSE325414 arm needs a donor table,
      a row map and panel triplets; gse325414_prepare writes all three from the deposit.

Steps:
  gse325414_prepare  wave1_gse325414_prepare.py: donor table from the series matrix, row map and
                     panel triplets from the deposit matrix
  arm_gse246215      wave1_arm_gse246215.py
  arm_gse325414      wave1_arm_gse325414.py (reads its inputs from scratch/recon)
  arm_sttt           wave1_arm_sttt.py

RUN: python atlas/03_arms_wave1.py [--list] [--dry-run] [--only step ...] [--from step]
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import run_steps  # noqa: E402

STEPS = [
    ("gse325414_prepare", ("wave1_gse325414_prepare.py", "all"), None),
    ("arm_gse246215", ("wave1_arm_gse246215.py",), None),
    ("arm_gse325414", ("wave1_arm_gse325414.py",), None),
    ("arm_sttt", ("wave1_arm_sttt.py",), None),
]

if __name__ == "__main__":
    sys.exit(run_steps("03", STEPS))
