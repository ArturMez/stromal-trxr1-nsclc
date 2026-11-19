#!/usr/bin/env python3
"""
Atlas stage 06: the object of the epithelial cells of primary LUSC and LUAD tumors at all genes (the cancer cells).

PAPER: the cancer cells used as a reference in Figure 2c (fourth column), Supplementary Figures 2a, 2c and 2d and the
       cancer-cell columns of the sheets Thioredoxin programs, Thioredoxin cycle genes and TXN and TXNIP of
       Supplementary Data 1. chapter2/02_working_sets.py --cells can reads this object (config.json key "epi_object")
       and keeps the cells that the original authors labeled as cancer cells.

WHAT: the atlas object of stage 04 holds fibroblasts, myeloid cells and T and NK cells at all genes. This stage builds
      the same kind of object for the epithelial cells. It does not write its own extraction code: it copies the
      full-transcriptome extraction and merge steps of stage 04 (steps/fullbuild_01*_extract_*.py and
      steps/fullbuild_02_merge.py) into the scratch folder with exact text replacements, each checked to occur as
      often as expected, and runs the copies in order:
        - the selection: level1_compartment epithelial_malignant, tissue_class primary_tumour, histology LUSC or LUAD
          (instead of the fibroblast, myeloid and T/NK compartments of all samples);
        - the outputs: work_root/atlas/output/RNA/meta_dataset/epibuild/ (arm_<study>_epi.h5ad, then
          epi_fullgene_object.h5ad and epi_gene_presence.csv); nothing of stage 04 is overwritten;
        - GSE246215 is skipped (a fibroblast-only deposit).
      The inputs are those of stage 04: the panel object and arm tables of stages 01 to 03 and the deposits under
      scrna_raw, so stages 00 to 04 must have run first.
FROM: work_root/atlas/output/RNA/meta_dataset/{stage1a,wave1,fullbuild/e6149}/ and raw_root/external_scRNA (or the
      config key scrna_raw), through atlas/common.py.
TO:   work_root/atlas/output/RNA/meta_dataset/epibuild/epi_fullgene_object.h5ad (set config.json "epi_object" to this
      path), the arms and epi_gene_presence.csv beside it; the generated step copies in work_root/atlas/scratch/
      epibuild_steps/; one log per step in work_root/atlas/logs/06_<step>.log.
RUN:  python atlas/06_cancer_cell_object.py [--list] [--dry-run] [--only step ...] [--from step]
      Several hours and tens of GB, like stage 04; no random numbers.
"""
import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

ATLAS = Path(__file__).resolve().parent
sys.path.insert(0, str(ATLAS))
from common import LOGS, SCRATCH  # noqa: E402

SRC = ATLAS / "steps"
DST = Path(SCRATCH) / "epibuild_steps"
MD = "{ROOT}/output/RNA/meta_dataset"           # as written inside the step scripts (an f-string of theirs)

# the cell selection of the extraction steps, and its epithelial version
SEL_OLD = "(P.level1_compartment.isin(TARGET))"
SEL_NEW = ('(P.level1_compartment.isin(TARGET)) & (P.tissue_class == "primary_tumour") '
           '& (P.histology.isin(["LUSC", "LUAD"]))')
# every copy imports atlas/common.py from the package, wherever the copy sits
COMMON = [("_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # atlas/",
           f'_sys.path.insert(0, "{ATLAS}")  # the package atlas/ (common.py)', 1)]
TARGET = [('TARGET = ["fibroblast", "myeloid", "T_NK"]', 'TARGET = ["epithelial_malignant"]', 1)]
OUT_DIR = [(f'OUT = f"{MD}/fullbuild"', f'OUT = f"{MD}/epibuild"', 1)]
# (old text, new text, expected count) per step file
PLAN = {
    "fullbuild_01_extract_luca.py": TARGET + OUT_DIR + [
        ("arm_luca_full.h5ad", "arm_luca_epi.h5ad", 2), (SEL_OLD, SEL_NEW, 1)],
    "fullbuild_01b_extract_emtab13526.py": TARGET + OUT_DIR + [
        ("arm_emtab13526_full.h5ad", "arm_emtab13526_epi.h5ad", 2), (SEL_OLD, SEL_NEW, 1)],
    "fullbuild_01c_extract_emtab6653.py": TARGET + OUT_DIR + [
        ("arm_emtab6653_full.h5ad", "arm_emtab6653_epi.h5ad", 2), (SEL_OLD, SEL_NEW, 1)],
    "fullbuild_01d_extract_emtab6149.py": TARGET + [
        (f'DST = f"{MD}/fullbuild/arm_emtab6149_full.h5ad"', f'DST = f"{MD}/epibuild/arm_emtab6149_epi.h5ad"', 1),
        (SEL_OLD, SEL_NEW, 1)],
    "fullbuild_01e_extract_wave1.py": TARGET + OUT_DIR + [
        ("arm_sttt_full.h5ad", "arm_sttt_epi.h5ad", 1), ("arm_gse325414_full.h5ad", "arm_gse325414_epi.h5ad", 1),
        ("arm_gse246215_full.h5ad", "arm_gse246215_epi.h5ad", 1),
        ("return D[D.level1_compartment.isin(TARGET)].copy()",
         'return D[D.level1_compartment.isin(TARGET) & (D.tissue_class == "primary_tumour") '
         '& D.histology.isin(["LUSC", "LUAD"])].copy()', 1),
        ('if which in ("all", "gse246215"):', 'if which == "gse246215":   # no epithelial cells in this deposit', 1)],
    "fullbuild_02_merge.py": [
        (f'D = f"{MD}/fullbuild"', f'D = f"{MD}/epibuild"', 1),
        ('os.environ.get("FULLGENE_DST", f"{D}/fullgene_object.h5ad")',
         'os.environ.get("EPI_DST", f"{D}/epi_fullgene_object.h5ad")', 1),
        ('_full.h5ad"),', '_epi.h5ad"),', 6), ('_full.h5ad")]', '_epi.h5ad")]', 1),
        ("fullgene_gene_presence.csv", "epi_gene_presence.csv", 2)],
}
# the steps in run order: name, generated script, arguments
STEPS = [("extract_luca", "epibuild_01_extract_luca.py", []),
         ("extract_emtab6149", "epibuild_01d_extract_emtab6149.py", []),
         ("extract_emtab6653", "epibuild_01c_extract_emtab6653.py", []),
         ("extract_wave1", "epibuild_01e_extract_wave1.py", ["all"]),
         ("extract_emtab13526", "epibuild_01b_extract_emtab13526.py", []),
         ("merge", "epibuild_02_merge.py", [])]


def generate():
    """Write the epithelial copies of the stage 04 steps to DST, every replacement checked by its count."""
    DST.mkdir(parents=True, exist_ok=True)
    for name, reps in PLAN.items():
        text = (SRC / name).read_text()
        for old, new, n in COMMON + reps:
            k = text.count(old)
            if k != n:
                raise SystemExit(f"{name}: expected {n} of {old!r}, found {k}")
            text = text.replace(old, new)
        if f'"{MD}/fullbuild"' in text.replace(f'{MD}/fullbuild/e6149', ""):
            raise SystemExit(f"{name}: an output path still points at the stage 04 folder")
        out = DST / name.replace("fullbuild_", "epibuild_")
        out.write_text(f"# Generated by atlas/06_cancer_cell_object.py from atlas/steps/{name}; do not edit.\n" + text)
    print(f"[06] wrote {len(PLAN)} step copies to {DST}", flush=True)


def main(argv=None):
    """Generate the step copies, then run the selected steps in order, each in its own process with a log."""
    ap = argparse.ArgumentParser(description="atlas stage 06: the epithelial-cell object")
    ap.add_argument("--list", action="store_true", help="print the steps and exit")
    ap.add_argument("--dry-run", action="store_true", help="print the commands without running them")
    ap.add_argument("--only", nargs="+", default=None, help="run only these steps, in stage order")
    ap.add_argument("--from", dest="start", default=None, help="run from this step to the end")
    a = ap.parse_args(argv)
    names = [s[0] for s in STEPS]
    for n in (a.only or []) + ([a.start] if a.start else []):
        if n not in names:
            ap.error(f"unknown step {n}; steps are {names}")
    todo = [s for s in STEPS if s[0] in a.only] if a.only else STEPS[names.index(a.start):] if a.start else STEPS
    cmds = [(name, [sys.executable, str(DST / script), *args]) for name, script, args in todo]
    if a.list or a.dry_run:
        for name, cmd in cmds:
            print(f"{name:20s} {' '.join(cmd)}")
        return 0
    generate()
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join([str(SRC), str(ATLAS), env.get("PYTHONPATH", "")])
    t_stage = time.time()
    for name, cmd in cmds:
        log = f"{LOGS}/06_{name}.log"
        print(f"[06] {name}: {' '.join(cmd)}  (log {log})", flush=True)
        t0 = time.time()
        with open(log, "w") as fh:
            rc = subprocess.call(cmd, stdout=fh, stderr=subprocess.STDOUT, env=env, cwd=str(SRC))
        print(f"[06] {name}: exit {rc} after {time.time() - t0:.0f} s", flush=True)
        if rc != 0:
            print(f"[06] stopped at {name}; see {log}", flush=True)
            return rc
    print(f"[06] done in {time.time() - t_stage:.0f} s", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
