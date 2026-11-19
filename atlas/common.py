"""
Paths and the step runner shared by every atlas stage.

WHAT: resolves the two folders the atlas build uses and runs the numbered steps of a stage in order.
FROM: config.json (or the file named by REDOX_CONFIG), the same file the rest of the package reads:
        raw_root      the downloaded public deposits live in raw_root/external_scRNA, laid out as
                      00_download.sh writes them (override with the optional key "scrna_raw")
        work_root     every file the build writes goes under work_root/atlas (override with the
                      optional key "atlas_root")
        rscript       optional; the R interpreter, found on PATH when not given
TO:   the constants below, imported by every step in atlas/steps/.
WHY:  the build steps take the project root from the configuration, not from their own text. The
      steps share one folder layout below the atlas root
      (output/RNA/meta_dataset/{stage1a,stage2,wave1,fullbuild}).

Layout under ROOT (= work_root/atlas):
  output/RNA/meta_dataset/          every table and object the steps write
  scratch/                          working files of the steps
  logs/                             one log per step, written by run_steps
"""
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

ATLAS = Path(__file__).resolve().parent
REPO = ATLAS.parent
STEPS = ATLAS / "steps"
DATA = str(ATLAS / "data")


def _load_config():
    """Read config.json (or the file named by REDOX_CONFIG) and check that raw_root and work_root are set."""
    named = os.environ.get("REDOX_CONFIG")
    path = Path(named).expanduser() if named else REPO / "config.json"
    if not path.exists():
        raise FileNotFoundError(f"no configuration at {path}; copy config.example.json to config.json")
    cfg = json.loads(path.read_text())
    for key in ("raw_root", "work_root"):
        if not cfg.get(key):
            raise KeyError(f"{path} has no {key}")
    return cfg


CFG = _load_config()
SCRNA = str(Path(CFG.get("scrna_raw") or Path(CFG["raw_root"]) / "external_scRNA").expanduser())
ROOT = str(Path(CFG.get("atlas_root") or Path(CFG["work_root"]) / "atlas").expanduser())
MD = f"{ROOT}/output/RNA/meta_dataset"
SCRATCH = f"{ROOT}/scratch"
LOGS = f"{ROOT}/logs"
RSCRIPT = CFG.get("rscript") or shutil.which("Rscript") or "Rscript"

for _d in (MD, f"{MD}/stage1a", f"{MD}/stage2", f"{MD}/wave1", f"{MD}/fullbuild",
           f"{MD}/fullbuild/blockqc_v2", SCRATCH, LOGS):
    os.makedirs(_d, exist_ok=True)


def step_cmd(step):
    """The command line of one step: (script, arg, ...) with {SCRATCH} and friends filled in."""
    script, *args = step
    fill = {"SCRATCH": SCRATCH, "MD": MD, "ROOT": ROOT, "SCRNA": SCRNA}
    args = [a.format(**fill) for a in args]
    path = STEPS / script
    if script.endswith(".R"):
        return [RSCRIPT, str(path), *args]
    return [sys.executable, str(path), *args]


def run_steps(stage, steps, argv=None):
    """Run the steps of a stage in order, each in its own process, stopping at the first failure.

    steps is a list of (name, (script, arg, ...), env) with env a dict or None. Command line:
      --list          print the steps and exit
      --only a b ...  run only the named steps (in stage order)
      --from name     run from this step to the end
      --dry-run       print the commands without running them
    """
    import argparse
    ap = argparse.ArgumentParser(description=f"atlas stage {stage}")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--only", nargs="+", default=None)
    ap.add_argument("--from", dest="start", default=None)
    a = ap.parse_args(argv)
    names = [s[0] for s in steps]
    for n in (a.only or []) + ([a.start] if a.start else []):
        if n not in names:
            ap.error(f"unknown step {n}; steps are {names}")
    todo = steps
    if a.only:
        todo = [s for s in steps if s[0] in a.only]
    elif a.start:
        todo = steps[names.index(a.start):]
    if a.list or a.dry_run:
        for name, step, env in todo:
            extra = " ".join(f"{k}={v.format(MD=MD, SCRATCH=SCRATCH)}" for k, v in (env or {}).items())
            print(f"{name:28s} {extra + ' ' if extra else ''}{' '.join(step_cmd(step))}")
        return 0
    t_stage = time.time()
    for name, step, env in todo:
        cmd = step_cmd(step)
        e = dict(os.environ)
        e["PYTHONPATH"] = os.pathsep.join([str(STEPS), str(ATLAS), e.get("PYTHONPATH", "")])
        for k, v in (env or {}).items():
            e[k] = v.format(MD=MD, SCRATCH=SCRATCH)
        log = f"{LOGS}/{stage}_{name}.log"
        print(f"[{stage}] {name}: {' '.join(cmd)}  (log {log})", flush=True)
        t0 = time.time()
        with open(log, "w") as fh:
            rc = subprocess.call(cmd, stdout=fh, stderr=subprocess.STDOUT, env=e, cwd=str(STEPS))
        print(f"[{stage}] {name}: exit {rc} after {time.time() - t0:.0f} s", flush=True)
        if rc != 0:
            print(f"[{stage}] stopped at {name}; see {log}", flush=True)
            return rc
    print(f"[{stage}] done in {time.time() - t_stage:.0f} s", flush=True)
    return 0
