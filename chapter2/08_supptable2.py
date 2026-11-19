#!/usr/bin/env python3
"""
Chapter 2, step 08: the tables behind the sheets of Supplementary Data 1 that come from the single-cell analysis.

WHAT: s1, the sheet Thioredoxin programs: one row per program (26, in the printed order and groups), its
          genes, and for each cell set the permutation z and the comparison-gene percentile of 06, whether
          the entry is bold (at or above the 97.5th percentile, z of at least 1.96 and, for programs of
          several genes, robust to leaving out any one gene) or n.t. (fewer than 3 genes detected, or a
          single gene that fails the detection rule).
      s2, the sheets Gene programs and Gene lists (the table the Supplementary Figure 3 footnote cites):
          one row per program of the panels of Figure 3 and Supplementary Figure 3 (every primary,
          contrast and technical-control program of a cell set's panel, figure_panels.json), with its
          source, reference, the genes of the list, the genes tested in each cell set whose panel holds
          it, and whether it is drawn (a program is drawn when, in at least one of its cell sets, the
          library-size split does not reproduce it, the rule of 05_plots.py; with fewer than 3 tested genes
          it is not scored); and every gene of every list, one row per gene.
      s3, the cancer cells and the additional tests as the workbook prints them: the columns of the cancer
          cells in the sheet Thioredoxin programs (all, LUSC, LUAD), the sheet Thioredoxin cycle genes and
          the sheet TXN and TXNIP. Written when 06 --cells can, 06 --only s4 and 11 have run.
      10_supplementary_data1.py writes the Excel workbook from these tables (the file names keep the numbers
      of the earlier table layout).
ITEMS: Supplementary Data 1, sheets Thioredoxin programs, Gene programs, Gene lists, Thioredoxin cycle genes
      and TXN and TXNIP (the sheet Read me is made by 10; the sheets Aldehyde enzyme genes and T and NK subtypes
      come from the tables of 04 and 12).
FROM: work_root/output/chapter2/trx_results.csv (06); gene_programs_definitions.csv and
      supptable_gene_programs.csv (01), program_results.csv and fig2_rows.csv (03).
TO:   work_root/output/chapter2/supptable2_programs.csv (s1), supptable2b_programs.csv and
      supptable2b_genes.csv (s2), supptable2_cancer.csv, supptable2_cycle_genes.csv and
      supptable2_txn_txnip.csv (s3), each with a .md sidecar.
WHY:  the sheets of the workbook need one tidy table each, with the values already rounded as printed;
      this step makes them from the test tables, so the workbook holds no values of its own.
RUN:  python3 chapter2/08_supptable2.py [--only s1 s2 s3]
      No random numbers.
"""
import argparse
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ch2lib as L  # noqa: E402


# ====================================================================== s1
def s1_table(a):
    """supptable2_programs.csv: z, percentile and bold flag of the 26 programs per cell set, as printed."""
    R = pd.read_csv(L.out("trx_results.csv"), float_precision="round_trip")
    ix = R.set_index(["cell", "set_id"])
    rows = []
    # 06 writes the sets of every cell set in the printed order, with their group and row label
    for first in R[R.cell == next(iter(L.CELL_SETS))].itertuples():
        k = first.set_id
        rec = dict(group=first.group, set_id=k, program=first.label, genes=first.listed)
        for c in L.CELL_SETS:
            r = ix.loc[(c, k)]
            nt = r.state == "nt"
            rec[f"{c}_z"] = np.nan if nt else r.z
            rec[f"{c}_percentile"] = np.nan if nt else r.percentile
            rec[f"{c}_bold"] = r.state == "hi"
            rec[f"{c}_z_printed"] = "n.t." if nt else f"{r.z:.1f}"
            rec[f"{c}_percentile_printed"] = "n.t." if nt else f"{r.percentile:.1f}"
        rows.append(rec)
    T = pd.DataFrame(rows)
    o = L.out("supptable2_programs.csv")
    T.to_csv(o, index=False)
    cols = [("group", "row group"), ("set_id", "set (06_thioredoxin.py)"), ("program", "row label"),
            ("genes", "the Genes column: the listed genes, or the source program")]
    for c, spec in L.CELL_SETS.items():
        cols += [(f"{c}_z", f"{spec['label']}: permutation z (empty when n.t.)"),
                 (f"{c}_percentile", f"{spec['label']}: percentile among genes detected as often as TXNRD1"),
                 (f"{c}_bold", f"{spec['label']}: True when printed in bold"),
                 (f"{c}_z_printed", "z as printed (1 dp, or n.t.)"),
                 (f"{c}_percentile_printed", "percentile as printed (1 dp, or n.t.)")]
    L.sidecar(o, "Supplementary Data 1, sheet Thioredoxin programs: TXNRD1 "
              "expression against the thioredoxin system, the NRF2 program and "
              "candidate upstream signals, per program and cell set (values from trx_results.csv).", cols)
    print(f"wrote supptable2_programs.csv, {len(T)} programs, "
          f"{int(sum(T[f'{c}_bold'].sum() for c in L.CELL_SETS))} bold entries")


# ====================================================================== s2
DASH = "-"                                     # printed for a cell set whose panel does not hold the program


def s2_programs(a):
    """supptable2b_programs.csv (one row per program) and supptable2b_genes.csv (one row per gene of each list)."""
    cfg = L.panels()
    G = pd.read_csv(L.out("supptable_gene_programs.csv"), dtype=str).fillna("")
    G["stem"] = G.program_id.map(L.stem)
    on_panel = {c: set(cfg["panel"][c]) for c in L.CELL_SETS}
    G = G[G.role.isin(cfg["suppfig3_roles"]) & G.stem.isin(set().union(*on_panel.values()))].copy()
    RS = pd.read_csv(L.out("program_results.csv"), dtype={"set_id": str})
    RS["stem"] = RS.set_id.map(L.stem)
    lab2cell = {s["label"]: c for c, s in L.CELL_SETS.items()}
    RS["cell"] = RS.cell_type.map(lab2cell)
    RS["rna"] = RS.verdict.str.startswith(("withdrawn", "not specific"))
    fig2 = set(pd.read_csv(L.out("fig2_rows.csv"), dtype={"set_id": str}).set_id.map(L.stem))
    group_order = {g: i for i, g in enumerate(cfg["groups"])}
    rows, genes = [], []
    for r in G.itertuples():
        cells = [c for c in L.CELL_SETS if r.stem in on_panel[c]]
        res = RS[(RS.stem == r.stem) & RS.cell.isin(cells) & (RS.role == r.role)]
        if res.empty:
            shown = "not drawn: fewer than 3 genes tested"
        elif (~res.rna).any():
            # the gene-program figure is Figure 3; the table of its rows keeps the name fig2_rows.csv
            shown = "Figure 3 and Supplementary Figure 3" if r.stem in fig2 else "Supplementary Figure 3"
        else:
            shown = "not drawn: the library-size split reproduces it in every cell set"
        m = re.match(r"^(.*?\d{4}, [^.]+)\.", r.reference)        # first author, year and journal
        ref = m.group(1) if m else r.reference
        rec = dict(group=r.group, program_id=r.program_id, program=r.label_as_printed, role=r.role,
                   source=r.source_type, identifier=r.source_identifier, reference=ref, pmid=r.pmid,
                   n_genes=int(r.n_genes_used), shown=shown)
        for c in L.CELL_SETS:
            n = int(getattr(r, f"n_tested_{c}"))
            rec[f"n_tested_{c}"] = n if c in cells else np.nan
            rec[f"{c}_printed"] = (f"{n}" if n >= L.MINSET else f"{n}, not scored") if c in cells else DASH
        rows.append(rec)
        tested = {c: set(getattr(r, f"not_testable_{c}").split()) for c in L.CELL_SETS}
        for g in r.genes.split():
            genes.append(dict(program_id=r.program_id, program=r.label_as_printed, gene=g,
                              **{f"tested_{c}": (int(g not in tested[c] and g != L.SPLITTER) if c in cells
                                                 else np.nan) for c in L.CELL_SETS}))
    T = pd.DataFrame(rows)
    T = T.assign(_g=T.group.map(group_order).fillna(len(group_order))).sort_values(["_g", "program"]).drop(columns="_g")
    o = L.out("supptable2b_programs.csv")
    T.to_csv(o, index=False)
    cols = [("group", "row group, as in Supplementary Figure 3"), ("program_id", "program id, as in every Chapter 2 table"),
            ("program", "label as printed"), ("role", "primary, contrast or technical control"),
            ("source", "kind of source"), ("identifier", "database identifier or PMC id"),
            ("reference", "first author, year and journal of the source"), ("pmid", "PubMed id"),
            ("n_genes", "genes in the list used (TXNRD1 included when listed; it is removed at scoring)"),
            ("shown", "where the program is drawn, or why it is not")]
    for c, spec in L.CELL_SETS.items():
        cols += [(f"n_tested_{c}", f"genes tested in {spec['short']}, empty when the program is not on their panel"),
                 (f"{c}_printed", f"the {spec['short']} cell as printed: the count, with ', not scored' below "
                                  f"{L.MINSET}, or '{DASH}' when the program is not on their panel")]
    L.sidecar(o, "Supplementary Data 1, sheet Gene programs: the gene "
              "programs of the panels of Figure 3 and Supplementary Figure 3, one row per program.", cols)
    Gn = pd.DataFrame(genes)
    og = L.out("supptable2b_genes.csv")
    Gn.to_csv(og, index=False)
    gcols = [("program_id", "program id"), ("program", "label as printed"), ("gene", "gene symbol of the list")]
    gcols += [(f"tested_{c}", f"1 when the gene is tested in {spec['short']}, 0 when not; empty when the program "
                              f"is not on their panel") for c, spec in L.CELL_SETS.items()]
    L.sidecar(og, "Supplementary Data 1, sheet Gene lists: the genes of the programs of the sheet Gene programs, one "
              "row per gene of each list (10_supplementary_data1.py adds the genes removed by a rule).", gcols)
    print(f"wrote supptable2b_programs.csv ({len(T)} programs; {T.shown.value_counts().to_dict()}) and "
          f"supptable2b_genes.csv ({len(Gn)} genes)")


# ====================================================================== s3
BEYOND = ("hi", "lo")       # printed bold: beyond the range in either direction, |z| >= 1.96 and robust


def printed(r):
    """z, percentile, bold and the two printed forms of one tested entry (a row with state, z and percentile)."""
    nt = r.state == "nt"
    return dict(z=np.nan if nt else r.z, percentile=np.nan if nt else r.percentile, bold=r.state in BEYOND,
                z_printed="n.t." if nt else f"{r.z:.1f}", percentile_printed="n.t." if nt else f"{r.percentile:.1f}")


def s3_additional(a):
    """The cancer cells and the additional tests as the workbook prints them. Written only when 06 and 11 have produced
    their tables (06 --cells can, 06 --only s4, 11); supptable2_programs.csv of s1 is not touched.
      supptable2_cancer.csv       the 26 programs in the cancer cells, all and by histology (Figure 2c, column
                                  Cancer cells, and Supplementary Figure 2c)
      supptable2_cycle_genes.csv  TXN, PRDX1, PRDX2 and SRXN1 one by one, the four cell sets (Supplementary Figure 2a)
      supptable2_txn_txnip.csv    TXNIP under the TXN split, the four cell sets (Supplementary Figure 2d)"""
    need = ["trx_results_with_extra.csv", "trx_cancer_histology.csv", "trx_cycle_genes.csv", "txn_txnip.csv"]
    missing = [n for n in need if not L.out(n).exists()]
    if missing:
        print(f"s3: not written, missing {missing} (run 06 --cells can, 06 --only s4 and 11_txn_txnip.py first)")
        return
    R = pd.read_csv(L.out("trx_results_with_extra.csv"), float_precision="round_trip")
    H = pd.read_csv(L.out("trx_cancer_histology.csv"), float_precision="round_trip")
    cols = {"can": R[R.cell == "can"].set_index("set_id"),
            "can_LUSC": H[H.histology == "LUSC"].set_index("set_id"),
            "can_LUAD": H[H.histology == "LUAD"].set_index("set_id")}
    head = {"can": "Cancer cells", "can_LUSC": "Cancer cells, LUSC", "can_LUAD": "Cancer cells, LUAD"}
    first = pd.read_csv(L.out("supptable2_programs.csv"), dtype=str, keep_default_na=False)
    rows = []
    for r0 in first.itertuples():
        rec = dict(group=r0.group, set_id=r0.set_id, program=r0.program, genes=r0.genes)
        for c, T in cols.items():
            rec.update({f"{c}_{k}": v for k, v in printed(T.loc[r0.set_id]).items()})
        rows.append(rec)
    C = pd.DataFrame(rows)
    o = L.out("supptable2_cancer.csv")
    C.to_csv(o, index=False)
    cc = [("group", "row group"), ("set_id", "set (06_thioredoxin.py)"), ("program", "row label"),
          ("genes", "the Genes column")]
    for c in cols:
        cc += [(f"{c}_z", f"{head[c]}: permutation z (empty when n.t.)"),
               (f"{c}_percentile", f"{head[c]}: percentile among genes detected as often as TXNRD1"),
               (f"{c}_bold", f"{head[c]}: True when printed in bold (beyond the 2.5th to 97.5th range, |z| >= 1.96, robust)"),
               (f"{c}_z_printed", "z as printed (1 dp, or n.t.)"),
               (f"{c}_percentile_printed", "percentile as printed (1 dp, or n.t.)")]
    L.sidecar(o, "Supplementary Data 1, sheet Thioredoxin programs, the columns of the cancer cells (the cells the "
              "original authors labeled as cancer cells): all, and those of LUSC and of LUAD (values from "
              "trx_results_with_extra.csv and trx_cancer_histology.csv). The rows are those of supptable2_programs.csv.",
              cc)
    G = pd.read_csv(L.out("trx_cycle_genes.csv"), float_precision="round_trip")
    cells = [c for c in list(L.CELL_SETS) + list(L.EXTRA_CELL_SETS) if c in set(G.cell)]
    gx = G.set_index(["cell", "set_id"])
    rows = []
    for k in dict.fromkeys(G.set_id):
        rec = dict(set_id=k, gene=k[2:])
        for c in cells:
            rec.update({f"{c}_{kk}": v for kk, v in printed(gx.loc[(c, k)]).items()})
        rows.append(rec)
    Gt = pd.DataFrame(rows)
    o = L.out("supptable2_cycle_genes.csv")
    Gt.to_csv(o, index=False)
    L.sidecar(o, "Supplementary Data 1, sheet Thioredoxin cycle genes: the four genes of the thioredoxin cycle scored one "
              "by one under the TXNRD1 split, per cell set (values from trx_cycle_genes.csv).",
              [("set_id", "set"), ("gene", "gene")] + [(f"{c}_...", f"{L.cell_spec(c)['label']}: z, percentile, bold and "
                                                       "the printed forms, as in supptable2_cancer.csv") for c in cells])
    X = pd.read_csv(L.out("txn_txnip.csv"), float_precision="round_trip").dropna(subset=["effect"])
    Xt = pd.DataFrame(dict(cell=X.cell, cell_type=X.cell_type, patients=X.patients.astype(int),
                           matched_pairs=X.matched_pairs.astype(int), effect=X.effect, z=X.z, percentile=X.percentile,
                           n_controls=X.n_controls.astype(int), tolerance_percent=(100 * X.tolerance).round().astype(int),
                           effect_printed=X.effect.map(lambda v: f"{v:.3f}"), z_printed=X.z.map(lambda v: f"{v:.1f}"),
                           percentile_printed=X.percentile.map(lambda v: f"{v:.1f}")))
    o = L.out("supptable2_txn_txnip.csv")
    Xt.to_csv(o, index=False)
    L.sidecar(o, "Supplementary Data 1, sheet TXN and TXNIP: the TXNIP score under the TXN split, per cell set (values "
              "from txn_txnip.csv).",
              [("cell", "cell set"), ("cell_type", "as printed"), ("patients", "patients"),
               ("matched_pairs", "pairs of a cell with high TXN and a partner"),
               ("effect", "mean paired difference in the TXNIP score, cells with high TXN minus partners"),
               ("z", "permutation z"), ("percentile", "percentile among the differences of the control genes"),
               ("n_controls", "control genes (mean count fraction close to TXN's)"),
               ("tolerance_percent", "window around TXN's mean count fraction that held them, percent"),
               ("effect_printed", "3 dp"), ("z_printed", "1 dp"), ("percentile_printed", "1 dp")])
    nb = int(sum(C[f"{c}_bold"].sum() for c in cols))
    print(f"wrote supptable2_cancer.csv ({len(C)} programs, {nb} bold entries), supptable2_cycle_genes.csv "
          f"({len(Gt)} genes, {len(cells)} cell sets) and supptable2_txn_txnip.csv ({len(Xt)} cell sets)")


def main():
    """Run the steps named with --only (default s1 to s3)."""
    ap = argparse.ArgumentParser(description=__doc__.split("RUN:")[0].strip(),
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", nargs="*", default=["s1", "s2", "s3"], choices=["s1", "s2", "s3"])
    a = ap.parse_args()
    t0 = time.time()
    if "s1" in a.only:
        s1_table(a)
    if "s2" in a.only:
        s2_programs(a)
    if "s3" in a.only:
        s3_additional(a)
    print(f"08_supptable2.py done in {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
