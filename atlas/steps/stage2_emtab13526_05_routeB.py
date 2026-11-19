#!/usr/bin/env python
"""
Atlas stage 02 (imported by the crosswalk step): build Route B for ALL 81 author sample groups.

Paper items: part of the atlas behind Figure 2, Figure 3, Supplementary Figures 2 to 4 and Supplementary Data 1.
Provides build_route_b and parse_patient_sample, used by stage2_emtab13526_09_final_crosswalk.py;
running it directly only prints the resolution. No random numbers are drawn.

Route B maps obs['sangerID'] to an SDRF library by matching the sangerID's
tokens against the FASTQ file names the SDRF records for each library
(Assay Name / Scan Name / SUBMITTED_FILE_NAME). It is completely independent of
the 'Patient N Tumor K' -> Pn_Tk name convention used by Route A.

Prints the full sangerID -> library resolution and every disagreement with
Route A. Writes nothing; the crosswalk step consumes the same logic.
"""
import os
import re

import pandas as pd

import os as _os  # noqa: E402
import sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # atlas/
from common import ROOT, SCRNA, SCRATCH, DATA, RSCRIPT  # noqa: E402,F401  paths from config.json
S2 = os.path.join(ROOT, "output/RNA/meta_dataset/stage2")
SDRF = os.path.join(SCRNA, "E-MTAB-13526",
                    "metadata/E-MTAB-13526.sdrf.txt")

pd.set_option("display.width", 260)
pd.set_option("display.max_rows", 400)
pd.set_option("display.max_colwidth", 70)


def sanger_tokens(sid):
    """Split a free-text sangerID into candidate sequencing-run tokens."""
    parts = [p for p in re.split(r"\s*-\s(?=\S)|\s*\|\s*", sid)]
    toks = []
    for p in parts:
        p = p.strip()
        if p:
            toks.append(p)
    return toks


def build_route_b():
    """Resolve each author sangerID to an SDRF library through the FASTQ file names the SDRF records."""
    s = pd.read_csv(SDRF, sep="\t", dtype=str, low_memory=False)
    s.columns = [c.strip() for c in s.columns]

    # every file-name string the SDRF associates with each library
    namecols = [c for c in ["Assay Name", "Scan Name", "Comment[SUBMITTED_FILE_NAME]"]
                if c in s.columns]
    lib_names = {}
    for lib, g in s.groupby("Source Name"):
        vals = set()
        for c in namecols:
            vals.update(str(v) for v in g[c].dropna().unique())
        lib_names[lib] = vals

    obs = pd.read_parquet(os.path.join(S2, "author_obs.parquet"),
                          columns=["patient_sample", "sangerID"])
    sid_of = (obs.groupby("patient_sample", observed=True)["sangerID"]
                 .agg(lambda x: sorted(set(map(str, x)))[0]).to_dict())

    rows = []
    for ps, sid in sorted(sid_of.items()):
        toks = sanger_tokens(sid)
        hits = {}
        for lib, names in lib_names.items():
            matched = set()
            for t in toks:
                # token must appear as a whole filename prefix, delimited by _ or .
                pat = re.compile(r"(?:^|[/_.])" + re.escape(t) + r"(?:[_.]|$)")
                if any(pat.search(n) for n in names):
                    matched.add(t)
            if matched:
                hits[lib] = matched
        rows.append({
            "patient_sample": ps,
            "sangerID": sid,
            "n_tokens": len(toks),
            "n_libs_matched": len(hits),
            "libs_matched": "|".join(sorted(hits)),
            "library_id_B": sorted(hits)[0] if len(hits) == 1 else None,
            "tokens_matched": "|".join(sorted({t for v in hits.values() for t in v})),
        })
    return pd.DataFrame(rows)


def parse_patient_sample(ps):
    """Parse an author 'Patient N Tumor K' (or Background / Healthy) name into the library naming convention."""
    m = re.fullmatch(r"Patient (\d+) Tumour (\d+)", ps)
    if m:
        return f"P{m.group(1)}_T{m.group(2)}"
    m = re.fullmatch(r"Patient (\d+) Background (\d+)", ps)
    if m:
        return f"P{m.group(1)}_B{m.group(2)}"
    m = re.fullmatch(r"Donor (\d+) Healthy (\d+)", ps)
    if m:
        return f"D{m.group(1)}_{m.group(2)}"
    return None


if __name__ == "__main__":
    rb = build_route_b()
    rb["library_id_A"] = rb["patient_sample"].map(parse_patient_sample)
    rb["agree"] = rb["library_id_A"] == rb["library_id_B"]

    print("Route B coverage")
    print("  groups:", len(rb))
    print("  resolved to exactly one library:", int(rb["library_id_B"].notna().sum()))
    print("  ambiguous (>1 library):", int((rb["n_libs_matched"] > 1).sum()))
    print("  unresolved (0 libraries):", int((rb["n_libs_matched"] == 0).sum()))
    print("  Route B library ids distinct:", rb["library_id_B"].nunique())

    print("\nUnresolved / ambiguous:")
    bad = rb[rb["library_id_B"].isna()]
    print(bad[["patient_sample", "sangerID", "n_libs_matched", "libs_matched"]].to_string(index=False)
          if len(bad) else "  none")

    print("\nDisagreements between Route A (name ordinal) and Route B (sequencing id):")
    d = rb[rb["library_id_B"].notna() & ~rb["agree"]]
    print(d[["patient_sample", "sangerID", "library_id_A", "library_id_B"]].to_string(index=False)
          if len(d) else "  none")

    print("\nFull resolution:")
    print(rb[["patient_sample", "sangerID", "library_id_A", "library_id_B",
              "n_libs_matched", "agree"]].to_string(index=False))
