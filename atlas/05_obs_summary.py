#!/usr/bin/env python3
"""
Atlas stage 05: the obs summary of the full-gene object, checked against the stored reference.

WHAT: counts cells, patients and QC-passing cells of fullgene_object.h5ad by block, study, source,
      histology, tissue class, level-1 compartment, level-2 label (within level 1), block_qc_pass,
      and for the chapter2/ selection (primary tumor, block_qc_pass, LUSC or LUAD) by histology and
      compartment; then compares the table with atlas/reference/fullgene_obs_summary.tsv. Both corrections of
      redox/atlas_fixes.py are applied first, as in every analysis: the study label is corrected, and the cells of
      the redistributed copies of the Lambrechts 2018 patients (the Lambrechts copies in the Zhang data, which
      are included once) are left out of every section and counted once, in the section left_out_copies, so
      that the summary still accounts for every cell of the object.
      Paper items: the atlas counts quoted in the Methods ("Atlas") and Supplementary Methods
      ("scRNA Atlas"): 13 public studies, 24 batches (blocks) of which one fails quality control.
FROM: output/RNA/meta_dataset/fullbuild/fullgene_object.h5ad under the atlas root (obs only, h5py;
      X is never opened), or the file given with --object.
TO:   output/RNA/meta_dataset/fullbuild/fullgene_obs_summary.tsv under the atlas root (or --out),
      and a verdict on stdout; exit code 1 when any count differs from the reference.
WHY:  a full rebuild of the atlas is tens of GB; this table is the target a rebuild must hit, and it
      carries every quantity chapter2/ takes from the atlas (blocks, patients, compartments, QC).
      The reference was written with --write-reference; its checksum is in the atlas README.

Columns of the summary (tab separated):
  section     which grouping: total, block_id, study, source_key, histology, tissue_class,
              level1_compartment, level2_celltype (key level1|level2), block_qc_pass,
              chapter2_selection (key histology|level1, cells of primary tumor that pass block QC),
              chapter2_selection_level2 (key histology|level1|level2, the same cells),
              left_out_copies (key all, the cells atlas_fixes.duplicate leaves out)
  key         the group (the obs values, with "tumour" written "tumor" as the package is in US English)
  n_cells     cells in the group
  n_patients  distinct obs patient_id in the group
  n_qc_pass   cells of the group with block_qc_pass True

RUN: python atlas/05_obs_summary.py [--object PATH] [--out PATH] [--reference PATH] [--write-reference]
"""
import argparse
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import common  # noqa: E402
from redox import atlas_fixes  # noqa: E402

REFERENCE = Path(__file__).resolve().parent / "reference" / "fullgene_obs_summary.tsv"
COLS = ["block_id", "study", "source_key", "histology", "tissue_class", "level1_compartment",
        "level2_celltype", "patient_id", "block_qc_pass"]
PRIMARY = "primary_tumour"   # the obs vocabulary of the object (British spelling)


def dec(a):
    """Decode an array of bytes or strings to an object array of str."""
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in a], dtype=object)


def col(h, name):
    """Read one obs column from the open h5ad file (categorical, nullable or plain) as an array."""
    g = h["obs"][name]
    if isinstance(g, h5py.Group):
        if "categories" in g:
            cats = dec(g["categories"][:])
            return cats[g["codes"][:]]
        v = g["values"][:].astype(object)
        v[g["mask"][:]] = None
        return v
    v = g[:]
    return dec(v) if v.dtype.kind in "SO" else v


def summarize(path):
    """Count cells, patients and QC-passing cells per group; returns the summary table (no random numbers)."""
    with h5py.File(path, "r") as h:
        o = pd.DataFrame({c: col(h, c) for c in COLS})
    for c in ("block_id", "study", "patient_id"):    # the corrected study label (redox.atlas_fixes)
        if c in o.columns:
            o[c] = atlas_fixes.relabel(o[c].to_numpy())
    o["block_qc_pass"] = o["block_qc_pass"].astype(str).isin(["True", "true", "1"])
    o["level2_celltype"] = o["level1_compartment"].astype(str) + "|" + o["level2_celltype"].astype(str)
    copies = atlas_fixes.duplicate(o["patient_id"].to_numpy())   # the Lambrechts copies (redox.atlas_fixes)
    left_out, o = o[copies], o[~copies]
    rows = []

    def add(section, frame, by):
        """Append the counts of frame grouped by the column by, as rows of the given section."""
        g = frame.groupby(by, dropna=False, observed=True)
        t = pd.DataFrame({"n_cells": g.size(), "n_patients": g.patient_id.nunique(),
                          "n_qc_pass": g.block_qc_pass.sum()}).reset_index()
        t = t.rename(columns={by: "key"})
        t.insert(0, "section", section)
        rows.append(t)

    tot = o.assign(all="all")
    add("total", tot, "all")
    for c in ("block_id", "study", "source_key", "histology", "tissue_class", "level1_compartment",
              "level2_celltype", "block_qc_pass"):
        add(c, o, c)
    sel = o[(o.tissue_class == PRIMARY) & o.block_qc_pass & o.histology.isin(["LUSC", "LUAD"])]
    sel = sel.assign(hist_comp=sel.histology.astype(str) + "|" + sel.level1_compartment.astype(str))
    add("chapter2_selection", sel, "hist_comp")
    sel_all = sel.assign(hist_comp=sel.histology.astype(str) + "|all")
    add("chapter2_selection", sel_all, "hist_comp")
    sel_l2 = sel.assign(hist_l2=sel.histology.astype(str) + "|" + sel.level2_celltype.astype(str))
    add("chapter2_selection_level2", sel_l2, "hist_l2")
    add("left_out_copies", left_out.assign(all="all"), "all")
    out = pd.concat(rows, ignore_index=True)
    # the object spells two tissue classes the British way; the summary writes them in US English
    out["key"] = out["key"].astype(str).str.replace("tumour", "tumor", regex=False)
    for c in ("n_cells", "n_patients", "n_qc_pass"):
        out[c] = out[c].astype(np.int64)
    return out.sort_values(["section", "key"]).reset_index(drop=True)


def main():
    """Build the summary, write it, and compare it with the reference (exit code 1 on a difference)."""
    ap = argparse.ArgumentParser(description=__doc__.split("Columns")[0].strip(),
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--object", default=f"{common.MD}/fullbuild/fullgene_object.h5ad")
    ap.add_argument("--out", default=f"{common.MD}/fullbuild/fullgene_obs_summary.tsv")
    ap.add_argument("--reference", default=str(REFERENCE))
    ap.add_argument("--write-reference", action="store_true",
                    help="write the summary as the new reference instead of checking against it")
    a = ap.parse_args()
    s = summarize(a.object)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    s.to_csv(a.out, sep="\t", index=False)
    tot = s[s.section == "total"].iloc[0]
    print(f"{a.object}: {tot.n_cells:,} cells, {tot.n_patients:,} patients, {tot.n_qc_pass:,} pass block QC; "
          f"{(s.section == 'block_id').sum()} blocks, {(s.section == 'study').sum()} studies")
    print(f"wrote {a.out}")
    if a.write_reference:
        s.to_csv(a.reference, sep="\t", index=False)
        print(f"wrote the reference {a.reference}")
        return 0
    ref = pd.read_csv(a.reference, sep="\t", dtype={"key": str})
    m = ref.merge(s, on=["section", "key"], how="outer", suffixes=("_ref", "_new"), indicator=True)
    bad = m[(m._merge != "both")]
    for c in ("n_cells", "n_patients", "n_qc_pass"):
        bad = pd.concat([bad, m[(m._merge == "both") & (m[f"{c}_ref"] != m[f"{c}_new"])]])
    bad = bad.drop_duplicates(["section", "key"])
    if len(bad):
        print(f"FAIL: {len(bad)} of {len(m)} rows differ from the reference")
        print(bad.to_string(index=False, max_rows=60))
        return 1
    print(f"PASS: all {len(m)} rows equal the reference {a.reference}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
