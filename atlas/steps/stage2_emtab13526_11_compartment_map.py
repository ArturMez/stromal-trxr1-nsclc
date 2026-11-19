#!/usr/bin/env python
"""
Atlas stage 02, step compartment_map: explicit author-label -> level1_compartment map for E-MTAB-13526.

Paper items: part of the atlas behind Figure 2, Figure 3, Supplementary Figures 2 to 4 and Supplementary Data 1.
Run: python atlas/steps/stage2_emtab13526_11_compartment_map.py   (normally through atlas/02_arm_emtab13526.py).
No random numbers are drawn.

EXPLICIT, never substring: "malignant cell" and "mast cell" both contain "t cell", and a substring
mapper would misfile them. Every one of the author's labels is listed literally below and the
script asserts that the map covers the observed vocabulary exactly, with no
label matched twice and none left over.

Labels the author themselves marked "(to remove)", and the unnamed leftover
cluster "8,1", are mapped but flagged so the build can drop them.

Writes: output/RNA/meta_dataset/stage2/compartment_map_emtab13526.csv
"""
import os

import pandas as pd

import os as _os  # noqa: E402
import sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # atlas/
from common import ROOT, SCRNA, SCRATCH, DATA, RSCRIPT  # noqa: E402,F401  paths from config.json
S2 = os.path.join(ROOT, "output/RNA/meta_dataset/stage2")

pd.set_option("display.width", 270)
pd.set_option("display.max_rows", 500)
pd.set_option("display.max_colwidth", 70)

# label : (level1_compartment, level2_celltype, confidence, note)
MAP = {
    # ---------------- T / NK ----------------
    "Naive T cells":                    ("T_NK", "CD4_T", "medium", "naive T, CD4/CD8 not resolved by the author label"),
    "CD8+ T cells":                     ("T_NK", "CD8_T", "high", ""),
    "CD8+ (low) T cells":               ("T_NK", "CD8_T", "high", ""),
    "Cycling CD8+ T cells":             ("T_NK", "CD8_T", "high", ""),
    "Cytotoxic T cells (low)":          ("T_NK", "CD8_T", "medium", "cytotoxic, lineage not stated"),
    "Cycling cytotoxic (low) T cells":  ("T_NK", "CD8_T", "medium", "cytotoxic, lineage not stated"),
    "Exhausted T cells":                ("T_NK", "CD8_T", "medium", "exhausted, lineage not stated"),
    "Exhausted cytotoxic T cells":      ("T_NK", "CD8_T", "medium", ""),
    "Exhausted cytotoxic T cells 2":    ("T_NK", "CD8_T", "medium", ""),
    "Cycling exhausted cytotoxic T cells": ("T_NK", "CD8_T", "medium", ""),
    "Downregulated T cells":            ("T_NK", "T_unspecified", "medium", "author QC-ish label, low-count T cells"),
    "Tregs":                            ("T_NK", "Treg", "high", ""),
    "Gamma delta T cells":              ("T_NK", "T_unspecified", "high", "gamma delta, no ct_high slot"),
    "NK cells":                         ("T_NK", "NK", "high", ""),
    "NK cells (higher cytotoxicity)":   ("T_NK", "NK", "high", ""),
    "NK cells (lower cytotoxicity)":    ("T_NK", "NK", "high", ""),

    # ---------------- B / plasma ----------------
    "B cells":                          ("B_plasma", "B", "high", ""),
    "LYZ+ B cells":                     ("B_plasma", "B", "medium", "LYZ is myeloid; possible myeloid doublet or ambient"),
    "Plasma B cells":                   ("B_plasma", "plasma", "high", ""),
    "Immature plasma B cells":          ("B_plasma", "plasma", "high", ""),
    "Cycling plasma B cells":           ("B_plasma", "plasma", "high", ""),

    # ---------------- myeloid ----------------
    "cDC2":                             ("myeloid", "DC", "high", ""),
    "pDCs":                             ("myeloid", "pDC", "high", ""),
    "Anti-inflammatory (high) alveolar macrophages":   ("myeloid", "macrophage", "high", "alveolar"),
    "Anti-inflammatory (medium) alveolar macrophages": ("myeloid", "macrophage", "high", "alveolar"),
    "Anti-inflammatory (low) alveolar macrophages":    ("myeloid", "macrophage", "high", "alveolar"),
    # NOTE the trailing space. It is in the data, it is the only form of this
    # label that occurs, and stripping it here would silently drop 1,929 cells.
    "Cycling anti-inflammatory alveolar macrophages ": ("myeloid", "macrophage", "high",
                                                        "author label carries a trailing space"),
    "Anti-inflammatory (high) macrophages":   ("myeloid", "macrophage", "high", ""),
    "Anti-inflammatory (medium) macrophages": ("myeloid", "macrophage", "high", ""),
    "Anti-inflammatory (low) macrophages":    ("myeloid", "macrophage", "high", ""),
    "Cycling anti-inflammatory macrophages":  ("myeloid", "macrophage", "high", ""),
    "Pro-inflammatory (high) macrophages":    ("myeloid", "macrophage", "high", ""),
    "Pro-inflammatory (medium) macrophages":  ("myeloid", "macrophage", "high", ""),
    "Pro-inflammatory (low) macrophages":     ("myeloid", "macrophage", "high", ""),
    "SPP1+ macrophages":                     ("myeloid", "macrophage", "high", ""),
    "SPP1+ cycling anti-inflammatory macrophages": ("myeloid", "macrophage", "high", ""),
    "Collagen+ macrophages":                 ("myeloid", "macrophage", "low",
                                              "UNSURE: collagen in a macrophage is a classic fibroblast doublet or ambient signature; 93.8 pct sit in CD235a- libraries"),
    "Immature myeloid cells":                ("myeloid", "monocyte", "medium", "immature, MDSC-gated material"),
    "Immature myeloid cells 1":              ("myeloid", "monocyte", "medium", ""),
    "Immature myeloid cells 2":              ("myeloid", "monocyte", "medium", ""),
    "Immature anti-inflammatory myeloid cells": ("myeloid", "monocyte", "medium", ""),
    "Cycling immature anti-inflammatory myeloid cells": ("myeloid", "monocyte", "medium", ""),
    "STAB+ anti-inflammatory cells":         ("myeloid", "macrophage", "low",
                                              "UNSURE: STAB1 is a scavenger-receptor macrophage marker, but the label omits the lineage word entirely"),
    "CAMLs":                                 ("myeloid", "macrophage", "low",
                                              "UNSURE: cancer-associated macrophage-like cells; myeloid lineage in the literature but a fused giant-cell entity, and 95.9 pct sit in CD235a- libraries"),

    # ---------------- epithelial ----------------
    "AT2 cells":                        ("epithelial_malignant", "AT2", "high", ""),
    "Club cells":                       ("epithelial_malignant", "club", "high", ""),
    "Ciliated epithelial cells":        ("epithelial_malignant", "ciliated", "high", ""),
    "Epithelial cells":                 ("epithelial_malignant", "epithelial_unspecified", "high", ""),
    "Transitioning epithelial cells":   ("epithelial_malignant", "epithelial_unspecified", "medium",
                                         "author term, not a standard state; candidate malignant or EMT-like"),
    "Cycling epithelial cells":         ("epithelial_malignant", "epithelial_unspecified", "high", ""),
    "Cycling epithelial cells 2":       ("epithelial_malignant", "epithelial_unspecified", "high", ""),

    # ---------------- fibroblast ----------------
    "Activated adventitial fibroblasts": ("fibroblast", "fibroblast_unspecified", "high",
                                          "the ONLY tumour fibroblast label in the atlas"),
    "Adventitial fibroblasts":           ("fibroblast", "fibroblast_unspecified", "high",
                                          "background object only"),

    # ---------------- endothelial ----------------
    "Lymphatic endothelial cells":      ("endothelial", "lymphatic_EC", "high",
                                         "no blood or capillary endothelium exists in this annotation"),

    # ---------------- other ----------------
    "Mast cells":                       ("other", "mast", "high", ""),
    "Cycling mast cells":               ("other", "mast", "high", ""),

    # ---------------- author-flagged junk ----------------
    "Mast cells,3 (to remove)":         ("other", "mast", "drop", "author marked 'to remove'"),
    "Mast cells,7 (to remove)":         ("other", "mast", "drop", "author marked 'to remove'"),
    "Mast cells,9 (to remove)":         ("other", "mast", "drop", "author marked 'to remove'"),
    "Cycling myeloid cells,2,0 (to remove)": ("myeloid", "monocyte", "drop", "author marked 'to remove'"),
    "8,1":                              ("other", "unassigned", "drop",
                                         "UNSURE: unnamed leftover cluster, a bare leiden id"),
}

# mural is deliberately absent: the annotation contains no pericyte or smooth
# muscle label at all.


def main():
    """Check that the map covers the observed author labels exactly and write compartment_map_emtab13526.csv."""
    obs = pd.read_parquet(os.path.join(S2, "author_obs.parquet"),
                          columns=["Cell types"])
    seen = obs["Cell types"].astype(str).value_counts()

    missing = sorted(set(seen.index) - set(MAP))
    extra = sorted(set(MAP) - set(seen.index))
    if missing:
        raise RuntimeError(f"author labels with no mapping: {missing}")
    if extra:
        raise RuntimeError(f"mapped labels not present in the data: {extra}")

    rows = []
    for lbl, (c1, c2, conf, note) in MAP.items():
        rows.append({"author_cell_type": lbl, "level1_compartment": c1,
                     "level2_celltype": c2, "confidence": conf,
                     "n_cells": int(seen[lbl]), "note": note})
    m = pd.DataFrame(rows).sort_values(["level1_compartment", "n_cells"],
                                       ascending=[True, False])
    m.to_csv(os.path.join(S2, "compartment_map_emtab13526.csv"), index=False)

    print(f"labels mapped: {len(m)}   cells: {m['n_cells'].sum():,}   "
          f"(author total {len(obs):,})")
    assert m["n_cells"].sum() == len(obs)
    assert m["author_cell_type"].nunique() == len(m)

    print("\ncells per level1_compartment")
    print(m.groupby("level1_compartment")["n_cells"].agg(["sum", "count"])
           .rename(columns={"sum": "cells", "count": "labels"})
           .sort_values("cells", ascending=False).to_string())

    print("\nlabels flagged UNSURE or drop")
    print(m[m["confidence"].isin(["low", "drop"])]
          [["author_cell_type", "level1_compartment", "n_cells", "confidence", "note"]]
          .to_string(index=False))

    print("\nfull map")
    print(m[["author_cell_type", "level1_compartment", "level2_celltype",
             "confidence", "n_cells"]].to_string(index=False))
    print("\nwrote", os.path.join(S2, "compartment_map_emtab13526.csv"))


if __name__ == "__main__":
    main()
