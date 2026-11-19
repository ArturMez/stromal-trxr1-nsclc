"""
Every path the package uses, resolved in one place.

WHAT: turns the two roots of the configuration file into the folders the scripts read and write.
FROM: the JSON file named by the REDOX_CONFIG environment variable, or config.json beside the
      repository root when that variable is not set.
TO:   accessors raw(), intermediate(), output(), input_file(), examples() and pellinen(), plus the
      parsed configuration in CFG.
WHY:  no script holds an absolute path. A script that needs a file asks for it here, so that
      moving the data means editing one JSON file and nothing else. Folders under
      work_root are created on demand, so a fresh clone with a config file and the data can run
      any script without a separate setup step.
"""
import json
import os
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

COHORTS = ("BOMI2", "BOMI1")


def _load_config():
    """Read the configuration file and check that the two required roots are present."""
    named = os.environ.get("REDOX_CONFIG")
    path = Path(named).expanduser() if named else REPO / "config.json"
    if not path.exists():
        raise FileNotFoundError(
            f"no configuration at {path}. Copy config.example.json to config.json and set the roots, "
            "or point REDOX_CONFIG at your own copy.")
    cfg = json.loads(path.read_text())
    for key in ("raw_root", "work_root"):
        if not cfg.get(key):
            raise KeyError(f"{path} has no {key}")
    cfg["_config_path"] = str(path)
    return cfg


CFG = _load_config()
RAW_ROOT = Path(CFG["raw_root"]).expanduser()
WORK_ROOT = Path(CFG["work_root"]).expanduser()


def raw(*parts):
    """A file under raw_root. Read only: nothing in the package ever writes here."""
    return RAW_ROOT.joinpath(*parts)


def intermediate(cohort, name=None):
    """work_root/intermediate/<cohort>/<name>, creating the folder. Without a name, the folder."""
    folder = WORK_ROOT / "intermediate" / str(cohort)
    folder.mkdir(parents=True, exist_ok=True)
    return folder if name is None else folder / name


def output(chapter, name=None):
    """work_root/output/chapter<N>/<name>, creating the folder. Without a name, the folder.

    The chapter may be given as a number (4) or as a folder name ("chapter4"); both land in the
    same place, so that a script may pass whichever reads better at the call site.
    """
    tag = str(chapter)
    if not tag.startswith("chapter"):
        tag = f"chapter{tag}"
    folder = WORK_ROOT / "output" / tag
    folder.mkdir(parents=True, exist_ok=True)
    return folder if name is None else folder / name


def input_file(key):
    """The file name of a raw input that config.json names under "input_files" (key -> file name).

    For example the Swedish clinical database of Uppsala II is key "bomi2_clinical_database".
    """
    named = (CFG.get("input_files") or {}).get(key)
    if not named:
        raise KeyError(f'{CFG["_config_path"]} has no "input_files" entry "{key}"')
    return named


def examples():
    """The example-core file (config.json "examples"). It ships with the data, never with the code."""
    named = CFG.get("examples")
    if not named:
        raise KeyError(f'{CFG["_config_path"]} has no "examples" entry')
    return Path(named).expanduser()


def pellinen(*parts):
    """The optional folder (config.json "pellinen") with the published fibroblast data of Pellinen and
    colleagues, or none when the configuration leaves it empty."""
    named = CFG.get("pellinen")
    return Path(named).expanduser().joinpath(*parts) if named else None
