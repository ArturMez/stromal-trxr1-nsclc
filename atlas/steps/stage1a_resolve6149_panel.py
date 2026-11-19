#!/usr/bin/env python
"""
Resolve the panel genes onto the Ensembl ids of the 2018 E-MTAB-6149 deposit.

Atlas stage 01, step resolve6149: resolve the panel genes in E-MTAB-6149.

Writes `panel_resolution.csv`, which stage1a_arm_emtab6149.py reads from its scratch folder. 26 panel
genes carry a different Ensembl id in the 2018 E-MTAB-6149 featuredata build than in the LuCA
reference (TXNIP and MRC1 among them); they are resolved by unique gene symbol, and 65 panel genes are
absent from the deposit.
Paper items: none directly; part of the atlas behind Figure 2, Figure 3, Supplementary Figures 2 to 4
and Supplementary Data 1. No random numbers are drawn.

RULE, per panel gene, in panel order:
  ensembl_id       the panel Ensembl id is a row of the 6149 matrix
  symbol_fallback  otherwise, exactly one matrix row carries the panel gene symbol in the
                   deposit featuredata; that row's Ensembl id is used
  absent           otherwise; written later as a structural zero

FROM: stage1a/panel_genes.csv (stage1a_panel.py) and stage1a/_features6149.csv
      (stage1a_dump6149_features.R: ensembl_id, gene_symbol, in_matrix for the deposit featuredata)
TO:   <scratch>/panel_resolution.csv  gene_symbol, ensembl_id, matrix_ensembl_id, resolution,
                                      source_set, is_hvg (one row per panel gene)
      <scratch>/panel_rows.txt        the matrix Ensembl ids to export, one per line, handed to
                                      stage1a_export6149.R as its second argument
CHECKS: 1,968 ensembl_id, 26 symbol_fallback and 65 absent on the 2,059-gene panel (the counts
      recorded in the arm's var gene_id_resolution); a different deposit or panel stops here.

RUN: python atlas/steps/stage1a_resolve6149_panel.py <scratch>   (normally through atlas/01_arms_luca_lambrechts.py)
"""
import os
import sys

import pandas as pd

import os as _os  # noqa: E402
import sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # atlas/
from common import ROOT, SCRATCH  # noqa: E402  paths from config.json

OUT = f"{ROOT}/output/RNA/meta_dataset/stage1a"
SCR = sys.argv[1] if len(sys.argv) > 1 else f"{SCRATCH}/e6149"
EXPECT = {"ensembl_id": 1968, "symbol_fallback": 26, "absent": 65}


def resolve(panel, feat):
    """One row per panel gene: how its expression is found in the deposit matrix."""
    inm = feat[feat.in_matrix.astype(str).isin(["True", "TRUE", "1"])]
    ids = set(inm.ensembl_id)
    by_sym = inm.groupby("gene_symbol").ensembl_id.agg(list)
    rows = []
    for r in panel.itertuples():
        if r.ensembl_id in ids:
            mid, how = r.ensembl_id, "ensembl_id"
        elif len(by_sym.get(r.gene_symbol, [])) == 1:
            mid, how = by_sym[r.gene_symbol][0], "symbol_fallback"
        else:
            mid, how = "", "absent"
        rows.append((r.gene_symbol, r.ensembl_id, mid, how, r.source_set, r.is_hvg))
    return pd.DataFrame(rows, columns=["gene_symbol", "ensembl_id", "matrix_ensembl_id", "resolution",
                                       "source_set", "is_hvg"])


def main():
    """Resolve each panel gene to a row of the E-MTAB-6149 matrix and write panel_resolution.csv and panel_rows.txt."""
    os.makedirs(SCR, exist_ok=True)
    panel = pd.read_csv(f"{OUT}/panel_genes.csv")
    feat = pd.read_csv(f"{OUT}/_features6149.csv")
    res = resolve(panel, feat)
    n = res.resolution.value_counts().to_dict()
    print(f"panel genes {len(res)}: {n}")
    assert n == EXPECT, f"resolution counts {n} differ from the recorded {EXPECT}"
    used = res.loc[res.resolution != "absent", "matrix_ensembl_id"]
    assert used.is_unique, "two panel genes resolve onto one matrix row"
    res.to_csv(f"{SCR}/panel_resolution.csv", index=False)
    with open(f"{SCR}/panel_rows.txt", "w") as fh:
        fh.write("\n".join(used) + "\n")
    print(f"wrote {SCR}/panel_resolution.csv and panel_rows.txt ({len(used)} ids)")
    print(res[res.resolution == "symbol_fallback"][["gene_symbol", "ensembl_id",
                                                    "matrix_ensembl_id"]].to_string(index=False))


if __name__ == "__main__":
    main()
