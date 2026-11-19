#!/usr/bin/env python3
"""
Chapter 2, step 10: Supplementary Data 1, the Excel workbook of the gene programs.

WHAT: one workbook, written from the tables of 08, 01, 04 and 12 without any new computation. Sheets:
        Read me               what each sheet holds and what its columns mean
        Thioredoxin programs  TXNRD1 expression against the thioredoxin system, the NRF2 program and
                              candidate upstream signals: the 26 programs of supptable2_programs.csv in
                              their order, with the values as printed (z and percentile to one decimal,
                              n.t. when not testable) and, per cell set, whether the entry is printed bold;
                              also the cancer cells (all, LUSC, LUAD; supptable2_cancer.csv)
        Thioredoxin cycle genes, TXN and TXNIP
                              the additional tests of Supplementary Figure 2a and 2d
                              (supptable2_cycle_genes.csv, supptable2_txn_txnip.csv). These two sheets and the
                              cancer columns are written when 08 s3 has produced their tables.
        Aldehyde enzyme genes the five genes of the aldehyde-metabolizing enzymes (a program of the sheet
                              Thioredoxin programs) tested one by one in fibroblasts: detection, pairs, effect,
                              permutation mean and SD and z, as printed, from fibroblast_genes_scan.csv of
                              04_fibroblast_genes.py, written when that table exists
        Gene programs         the 30 programs of the panels of Figure 3 and Supplementary Figure 3, in their
                              order: the label with "(contrast)" for the contrast program, the source as
                              the database identifier or the location in the publication, the reference
                              with its PMID, the genes, the genes tested per cell set and where the program
                              is drawn, each built from supptable2b_programs.csv and the where_in_source
                              column of supptable_gene_programs.csv
        T and NK subtypes     the four programs of T and NK cells in Figure 3, in all T and NK cells and within
                              conventional CD4 T cells, regulatory T cells, CD8 T cells and NK cells, with all cells
                              and without the cells of one study (PMID 29942094): tnk_subtypes.csv of
                              12_tnk_subtypes.py, written when that table exists
        Gene lists            every gene of every program of the Gene programs sheet: the genes used, from
                              supptable2b_genes.csv, with the genes tested per cell set, and the edits to
                              the published lists from supptable_gene_programs.csv (genes removed by a rule,
                              with the reason, and symbols replaced by the atlas symbol)
      The program labels of the Gene lists sheet are those of the Gene programs sheet.
ITEMS: Supplementary Data 1 (the file Supplementary_Data_1.xlsx).
FROM: work_root/output/chapter2/supptable2_programs.csv, supptable2b_programs.csv, supptable2b_genes.csv
      (08_supptable2.py) and the where_in_source, edits and symbols_renamed columns of
      supptable_gene_programs.csv (01_programs.py s2), fibroblast_genes_scan.csv (04_fibroblast_genes.py s1) and
      tnk_subtypes.csv (12_tnk_subtypes.py).
TO:   work_root/output/supplementary/Supplementary_Data_1.xlsx and its sibling Supplementary_Data_1.md.
      The cell values are reproducible; the file bytes are not, because openpyxl stamps the save time
      into the workbook properties.
WHY:  the program tables are reference data, read in a spreadsheet rather than on paper, so they are
      supplied as one workbook.
RUN:  python3 chapter2/10_supplementary_data1.py
      Run after 04_fibroblast_genes.py, 08_supptable2.py and 12_tnk_subtypes.py. No random numbers.
"""
import re
import sys
import time
from pathlib import Path

import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Alignment, Font
from openpyxl.utils import get_column_letter

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ch2lib as L  # noqa: E402
from redox import paths  # noqa: E402

NAME = "Supplementary_Data_1"
NOT_TESTED = "n.t."
# One reason and its genes per edit, joined by "; " in the edits column of 01_programs.py
# ("<reason>: <gene> <gene> ..."). A reason may itself hold ": " or "; ", so the genes are found as the
# run of gene symbols at the end of each part; the counts are checked against the list sizes below.
EDIT = re.compile(r"(.+?): ((?:[A-Z0-9][A-Za-z0-9_.-]*)(?: [A-Z0-9][A-Za-z0-9_.-]*)*)(?:; |$)")


def target(name):
    """work_root/output/supplementary/<name>, creating the folder."""
    folder = paths.WORK_ROOT / "output" / "supplementary"
    folder.mkdir(parents=True, exist_ok=True)
    return folder / name


def read(name):
    """A table of the chapter 2 output folder, every value as text."""
    return pd.read_csv(L.out(name), dtype=str, keep_default_na=False)


def printed_number(text):
    """A value as the table printed it: the number when it is one, else the text (n.t., a dash)."""
    try:
        return float(text)
    except ValueError:
        return text


def printed_count(text):
    """A count as the table printed it: a whole number when it is one, else the text (a dash, "2, not scored")."""
    return int(text) if text.isdigit() else text


# The Source column of the sheet Gene programs: the database identifier when there is one, else the location in the
# publication, the first part of where_in_source ("Table 1", "Supplementary Table S5", "Methods section").
DATABASE = re.compile(r"(MSigDB|Reactome|GO:|KEGG)")
# The Drawn in column of the sheet Gene programs, from the shown column of 08_supptable2.py.
DRAWN_IN = {"Figure 3 and Supplementary Figure 3": "Figure 3, Supplementary Figure 3",
            "Supplementary Figure 3": "Supplementary Figure 3",
            "not drawn: the library-size split reproduces it in every cell set": "not drawn (RNA amount)",
            "not drawn: fewer than 3 genes tested": "not drawn (not scored)"}


def program_label(program, role):
    """The program label of the sheet Gene programs: the contrast program carries "(contrast)"."""
    return program + (" (contrast)" if role == "contrast" else "")


def source_location(identifier, where_in_source):
    """The Source column: the database identifier, else the location in the publication."""
    if DATABASE.match(identifier):
        return identifier
    where = re.split(r" \(|, |;", where_in_source)[0].strip()
    if not where:
        raise ValueError(f"no location in the publication: {where_in_source!r}")
    return where


def parse_edits(text):
    """The rule-based edits of one program as (reason, genes) pairs, in the order of the edits column."""
    out, pos = [], 0
    while pos < len(text):
        m = EDIT.match(text, pos)
        if m is None:
            raise ValueError(f"cannot read the edits column: {text[pos:]!r}")
        out.append((m.group(1), m.group(2).split()))
        pos = m.end()
    return out


# ---------------------------------------------------------------- the sheets
CANCER = {"can": "Cancer cells", "can_LUSC": "Cancer cells, LUSC", "can_LUAD": "Cancer cells, LUAD"}


def extra():
    """True when 08 s3 has written the tables of the cancer cells and the additional tests."""
    return all(L.out(n).exists() for n in ("supptable2_cancer.csv", "supptable2_cycle_genes.csv",
                                             "supptable2_txn_txnip.csv"))


def triple(r, c):
    """z, percentile and the bold mark of one entry as printed, and which of the three cells are bold."""
    is_bold = getattr(r, f"{c}_bold") == "True"
    return ([printed_number(getattr(r, f"{c}_z_printed")), printed_number(getattr(r, f"{c}_percentile_printed")),
             "bold" if is_bold else ""], [is_bold, is_bold, False])


def sheet_thioredoxin():
    """The sheet Thioredoxin programs: header, rows and the bold marks, with the cancer columns when 08 s3 has run."""
    T = read("supptable2_programs.csv")
    cols = {c: spec["label"] for c, spec in L.CELL_SETS.items()}
    if extra():
        C = read("supptable2_cancer.csv")
        if C.set_id.tolist() != T.set_id.tolist():
            raise ValueError("supptable2_cancer.csv does not hold the programs of supptable2_programs.csv in their order")
        T = T.merge(C.drop(columns=["group", "program", "genes"]), on="set_id", validate="one_to_one")
        cols.update(CANCER)
    header = ["Group", "Program", "Genes"]
    for label in cols.values():
        header += [f"{label}, z", f"{label}, percentile", f"{label}, bold"]
    rows, bold = [], []
    for r in T.itertuples():
        row, marks = [r.group, r.program, r.genes], []
        for c in cols:
            cells, m = triple(r, c)
            row += cells
            marks += m
        rows.append(row)
        bold.append([False] * 3 + marks)
    return header, rows, bold


def sheet_cycle_genes():
    """The four genes of the thioredoxin cycle one by one, the columns of the sheet Thioredoxin programs."""
    G = read("supptable2_cycle_genes.csv")
    cells = [c for c in list(L.CELL_SETS) + list(L.EXTRA_CELL_SETS) if f"{c}_z_printed" in G.columns]
    header = ["Gene"]
    for c in cells:
        label = L.cell_spec(c)["label"]
        header += [f"{label}, z", f"{label}, percentile", f"{label}, bold"]
    rows, bold = [], []
    for r in G.itertuples():
        row, marks = [r.gene], []
        for c in cells:
            cc, m = triple(r, c)
            row += cc
            marks += m
        rows.append(row)
        bold.append([False] + marks)
    return header, rows, bold


def sheet_txn_txnip():
    """The sheet TXN and TXNIP: header and one row per cell type."""
    X = read("supptable2_txn_txnip.csv")
    header = ["Cell type", "Patients", "Pairs", "Difference in TXNIP score", "z", "Percentile", "Comparison genes",
              "Window, %"]
    rows = [[r.cell_type, int(r.patients), int(r.matched_pairs), float(r.effect_printed), float(r.z_printed),
             float(r.percentile_printed), int(r.n_controls), int(r.tolerance_percent)] for r in X.itertuples()]
    return header, rows


# The sheet T and NK subtypes: blocks in this order, the four programs in the order of Figure 3.
SUBTYPE_ORDER = ["all", "CD4_T", "Treg", "CD8_T", "NK"]
PROGRAM_ORDER = ["myc_targets", "mtorc1", "cell_cycle", "treg"]
WITHOUT = "without_Guo_Zhang_2018"


def tnk():
    """True when 12_tnk_subtypes.py has written tnk_subtypes.csv."""
    return L.out("tnk_subtypes.csv").exists()


def three(v):
    """A value as printed with three decimals."""
    return float(f"{v:.3f}")


def p_printed(p):
    """A p value as printed: three decimals, and "< 0.001" below that."""
    return "< 0.001" if p < 0.001 else three(p)


def sheet_tnk_subtypes():
    """The four programs of T and NK cells in Figure 3, in all T and NK cells and within each subtype (20 rows), with
    the same test without the cells of one study (Guo_Zhang_2018) in the last five columns."""
    T = pd.read_csv(L.out("tnk_subtypes.csv"), float_precision="round_trip")
    header = ["Subtype", "Program", "Patients", "Pairs", "Batches", "Effect", "Permutation SD", "p", "Same direction",
              "Without one study: pairs", "Without one study: effect", "Without one study: permutation SD",
              "Without one study: p", "Without one study: same direction"]
    rows = []
    for st in SUBTYPE_ORDER:
        for pr in PROGRAM_ORDER:
            a = T[(T.run == "all") & (T.subtype == st) & (T.program == pr)]
            b = T[(T.run == WITHOUT) & (T.subtype == st) & (T.program == pr)]
            if len(a) != 1 or len(b) != 1:
                raise ValueError(f"tnk_subtypes.csv: no single row for {st}, {pr}")
            a, b = a.iloc[0], b.iloc[0]
            rows.append([a.subtype_label, a.program_label, int(a.patients), int(a.pairs), int(a.batches),
                         three(a.effect), three(a.perm_sd), p_printed(a.p_gauss),
                         f"{int(a.batches_same_direction)} of {int(a.batches_with_20_pairs)}",
                         int(b.pairs), three(b.effect), three(b.perm_sd), p_printed(b.p_gauss),
                         f"{int(b.batches_same_direction)} of {int(b.batches_with_20_pairs)}"])
    return header, rows


# The sheet Aldehyde enzyme genes: the genes of this program of the sheet Thioredoxin programs, in that sheet's order.
ALDEHYDE_SET = "C1_aldehyde"


def scan():
    """True when 04_fibroblast_genes.py has written fibroblast_genes_scan.csv."""
    return L.out("fibroblast_genes_scan.csv").exists()


def one(v):
    """A value as printed with one decimal; a value that rounds to zero is 0.0, without a minus sign."""
    v = float(f"{v:.1f}")
    return 0.0 if v == 0 else v


def sheet_aldehyde_genes():
    """The genes of the aldehyde-metabolizing enzymes tested one by one in fibroblasts (fibroblast_genes_scan.csv), in
    the order of the sheet Thioredoxin programs: detection in percent, pairs, effect, permutation mean and SD, and z."""
    T = read("supptable2_programs.csv")
    row = T[T.set_id == ALDEHYDE_SET]
    if len(row) != 1:
        raise ValueError(f"supptable2_programs.csv: no single row for {ALDEHYDE_SET}")
    genes = row.iloc[0].genes.split(", ")
    G = pd.read_csv(L.out("fibroblast_genes_scan.csv"), float_precision="round_trip")
    if G.gene.isna().any() or G.gene.duplicated().any():
        raise ValueError("fibroblast_genes_scan.csv: a gene is missing or listed twice")
    G = G.set_index("gene")
    missing = [g for g in genes if g not in G.index]
    if missing:
        raise ValueError(f"fibroblast_genes_scan.csv has no row for {', '.join(missing)}")
    header = ["Gene", "Detected, %", "Pairs", "Effect", "Permutation mean", "Permutation SD", "z"]
    rows = [[g, float(f"{100 * G.loc[g, 'detection']:.1f}"), int(G.loc[g, "n_pairs"]), three(G.loc[g, "eff"]),
             three(G.loc[g, "perm_mean"]), three(G.loc[g, "perm_sd"]), one(G.loc[g, "z"])] for g in genes]
    return header, rows


def sheet_programs():
    """The sheet Gene programs, one row per program, its row group in the first column."""
    T = read("supptable2b_programs.csv")
    where = read("supptable_gene_programs.csv").set_index("program_id").where_in_source
    header = ["Group", "Program", "Source", "Reference", "Genes"] + [spec["label"] for spec in L.CELL_SETS.values()]
    header += ["Drawn in"]
    rows = []
    for r in T.itertuples():
        reference = f"{r.reference} (PMID {r.pmid})" if r.pmid else r.reference
        rows.append([r.group, program_label(r.program, r.role), source_location(r.identifier, where[r.program_id]),
                     reference, int(r.n_genes)]
                    + [printed_count(getattr(r, f"{c}_printed")) for c in L.CELL_SETS] + [DRAWN_IN[r.shown]])
    return header, rows, T.program_id.tolist(), dict(zip(T.program_id, rows))


def sheet_genes(order, label):
    """label: program id to the program label of the Gene programs sheet."""
    G = read("supptable2b_genes.csv")
    D = read("supptable_gene_programs.csv").set_index("program_id")
    header = ["Program", "Gene", "Edit to the published list"] + [f"Tested, {spec['short']}"
                                                                  for spec in L.CELL_SETS.values()]
    word = {"1.0": "yes", "0.0": "no", "": "not on panel"}
    rows, n_used, n_removed = [], 0, 0
    for pid in order:
        d = D.loc[pid]
        renamed = dict(p.split(" -> ") for p in d.symbols_renamed.split("; ") if p)
        source_of = {new: old for old, new in renamed.items()}
        g = G[G.program_id == pid]
        if len(g) != int(d.n_genes_used):
            raise ValueError(f"{pid}: {len(g)} genes in supptable2b_genes.csv, {d.n_genes_used} in the definitions")
        for r in g.itertuples():
            edit = (f"the source names it {source_of[r.gene]}; replaced by the atlas symbol (HGNC)"
                    if r.gene in source_of else "")
            rows.append([label[pid], r.gene, edit] + [word[getattr(r, f"tested_{c}")] for c in L.CELL_SETS])
        removed = 0
        for reason, genes in parse_edits(d.edits):
            for gene in genes:
                rows.append([label[pid], gene, f"removed: {reason}"] + [""] * len(L.CELL_SETS))
            removed += len(genes)
        # every gene of the published list is either used (once, after the symbol replacement) or removed
        if int(d.n_genes_used) + removed != int(d.n_genes_in_source):
            raise ValueError(f"{pid}: {d.n_genes_used} used + {removed} removed != {d.n_genes_in_source} in the source")
        n_used += len(g)
        n_removed += removed
    return header, rows, n_used, n_removed


# The lines that describe the cancer cells and the additional tests are written only when their tables exist
# ("+" marks them); "t" marks the lines of the sheet T and NK subtypes and "a" those of the sheet Aldehyde enzyme genes.
README = [
    ("Supplementary Data 1", None),
    ("Gene programs of the single-cell analysis of TXNRD1-positive cells (Figures 2 and 3, Supplementary "
     "Figures 2 and 3).", None),
    ("", None),
    ("Sheet Thioredoxin programs", "bold"),
    ("TXNRD1 expression against the thioredoxin system, the NRF2 program and candidate upstream signals, one row "
     "per program (Figure 2c).", None),
    ("Group: the row group. Program: the row label. Genes: the genes of the list, or the published program it "
     "takes (Hallmark gene sets; the S-phase and G2/M genes of Tirosh and colleagues).", None),
    ("For each cell type (fibroblasts, macrophages, T and NK cells): z, the permutation z of the difference "
     "between TXNRD1-positive cells and their matched partners; percentile, the rank of this difference among the "
     "differences obtained when the same cells were split by genes detected as often as TXNRD1 (Methods); bold, "
     "marked when the entry is at or above the 97.5th percentile or at or below the 2.5th percentile, with an "
     "absolute z of at least 1.96, and robust to leaving "
     "out any one gene. n.t., not testable: fewer than three genes detected or, for a single gene, measured "
     "in fewer than 60 % or detected in fewer than 2 % of the cells (GPX1 is missing from the count tables of the "
     "largest studies). SELENOK, SELENOW, SELENOH and GPX1 include the counts of their former symbols in the "
     "studies that use them (SELK, SEPW1, C11orf31 and an Ensembl identifier).", None),
    ("Cancer cells: the cells labeled as cancer cells by the original authors, analyzed as a reference in the same "
     "way, all together (Figure 2c) and those of LUSC and of LUAD separately (Supplementary Figure 2c). In cancer "
     "cells the thioredoxin cycle includes SRXN1, which was not testable in the other cell types.", "+"),
    ("", "+"),
    ("Sheet Thioredoxin cycle genes", "+bold"),
    ("The genes of the thioredoxin cycle tested one by one (Supplementary Figure 2a), with z, percentile and bold "
     "as in the sheet Thioredoxin programs.", "+"),
    ("", "+"),
    ("Sheet TXN and TXNIP", "+bold"),
    ("TXNIP in cells with high TXN (Supplementary Figure 2d), one row per cell type. Cells with a TXN count "
     "fraction above the median of their batch and histology were paired with the other cells in the same way as "
     "for TXNRD1. Patients, Pairs: the patients and the matched pairs. Difference in TXNIP score: the mean paired "
     "difference, cells with high TXN minus their partners. z: its permutation z. Percentile: its rank among the "
     "differences obtained when the same cells were split in the same way by genes expressed as strongly as TXN. "
     "Comparison genes, Window: the number of these genes and the window around the mean count fraction of TXN "
     "that held them.", "+"),
    ("", "a"),
    ("Sheet Aldehyde enzyme genes", "abold"),
    ("The five genes of the aldehyde-metabolizing enzymes (sheet Thioredoxin programs, Figure 2c) tested one by one in "
     "fibroblasts, one row per gene. The values come from the test of single genes in fibroblasts (Supplementary "
     "Methods), with 100 permutations. A program score is the mean of the scores of its genes, so the effects of the "
     "genes show which of them carry the difference of the set. None of the five genes passed the Benjamini-Hochberg "
     "correction across all genes of that test (q above 0.20).", "a"),
    ("Detected, %: the fibroblasts in which the gene was detected, as a percentage of those in which it was measured. "
     "Pairs: the matched pairs with a value for the gene in both cells. Effect: the mean paired difference in the "
     "within-batch z score of the gene's count fraction, TXNRD1-positive cells minus their partners. Permutation mean "
     "and Permutation SD: the mean and the SD of this difference over 100 permutations of TXNRD1 status among the "
     "cells of each patient. z: the effect minus the permutation mean, divided by the permutation SD.", "a"),
    ("", None),
    ("Sheet Gene programs", "bold"),
    ("The programs of the panel of each cell type in Figure 3 and Supplementary Figure 3: published cell-state "
     "signatures and pathway gene sets, one contrast (MHC class I genes) and one technical control "
     "(dissociation stress score), one row per program.", None),
    ("Group: as in Supplementary Figure 3. Program: the program label. Source: the database identifier or, for a "
     "list taken from a publication, where it is in that publication. Reference: the publication, with its PubMed "
     "ID. Genes: the genes in the list as used, after the edits of the Gene lists sheet (TXNRD1 included when "
     "listed; it is removed at scoring).", None),
    ("Fibroblasts, Macrophages, T and NK cells: the genes of the list tested in that cell type (Methods); a dash "
     "marks a program that is not in the panel of that cell type, and not scored marks fewer than three tested "
     "genes. Drawn in: where the program is shown; programs whose difference the library-size split reproduced in "
     "every cell type of their panel (RNA amount) are listed but not drawn.", None),
    ("", None),
    ("Sheet T and NK subtypes", "tbold"),
    ("The four programs drawn for T and NK cells in Figure 3, tested within conventional CD4 T cells, regulatory T "
     "cells, CD8 T cells and NK cells, one row per subtype and program. The first four rows repeat the test in all T "
     "and NK cells (Figure 3). The cells of each subtype were paired and tested in the same way as all T and NK "
     "cells. Cells without one of the four subtype labels were not included, among them the cells of one sequencing "
     "batch that has no subtype labels.", "t"),
    ("Patients, Pairs, Batches: the patients, the matched pairs and the sequencing batches. Effect: the mean paired "
     "difference in program score, TXNRD1-positive cells minus their partners. Permutation SD and p: from permutation "
     "of TXNRD1 status among the cells of each patient, p not corrected for the number of tests. Same direction: the "
     "batches whose difference has the sign of the effect, among the batches with at least 20 pairs. Without one "
     "study: the same test without the cells of Guo and colleagues (PMID 29942094).", "t"),
    ("", "t"),
    ("Sheet Gene lists", "bold"),
    ("Every gene of every program of the Gene programs sheet, one row per gene and program, in the order of the "
     "published list; then the genes of the published list that a rule removed.", None),
    ("Program: as in the Gene programs sheet. Gene: the gene symbol. Edit to the published list: empty when "
     "the gene is used as published; for a used gene whose published symbol is not the symbol of the atlas, "
     "the published symbol; for a removed gene, the rule and its reason.", None),
    ("Tested, per cell type: yes when the gene is tested in that cell type, no when it is not (detected in "
     "too few cells, measured in too few cells, or TXNRD1 itself), not on panel when the program is not in the "
     "panel of that cell type; empty for a removed gene.", None),
]


def write_sheet(ws, header, rows, bold=None, widths=None):
    """Fill a worksheet: bold header, rows, optional bold cells, column widths, floats to one decimal, frozen header row."""
    ws.append(header)
    for cell in ws[1]:
        cell.font = Font(bold=True)
        cell.alignment = Alignment(wrap_text=True, vertical="top")
    for i, row in enumerate(rows):
        ws.append(row)
        if bold is not None:
            for j, b in enumerate(bold[i]):
                if b:
                    ws.cell(row=i + 2, column=j + 1).font = Font(bold=True)
    for j, w in enumerate(widths or []):
        ws.column_dimensions[get_column_letter(j + 1)].width = w
    for row in ws.iter_rows(min_row=2):
        for cell in row:
            if isinstance(cell.value, float):
                cell.number_format = "0.0"
    ws.freeze_panes = "A2"


def main():
    """Build the workbook and its .md description."""
    t0 = time.time()
    wb = Workbook()
    ws = wb.active
    ws.title = "Read me"
    more = extra()
    has_tnk = tnk()
    has_scan = scan()
    for text, style in README:
        if style and style.startswith("+"):
            if not more:
                continue
            style = style[1:]
        if style and style.startswith("t"):            # the lines of the sheet T and NK subtypes
            if not has_tnk:
                continue
            style = style[1:]
        if style and style.startswith("a"):            # the lines of the sheet Aldehyde enzyme genes
            if not has_scan:
                continue
            style = style[1:]
        ws.append([text])
        c = ws.cell(row=ws.max_row, column=1)
        c.alignment = Alignment(wrap_text=True, vertical="top")
        if style == "bold" or text == "Supplementary Data 1":
            c.font = Font(bold=True)
    ws.column_dimensions["A"].width = 120

    h1, r1, b1 = sheet_thioredoxin()
    write_sheet(wb.create_sheet("Thioredoxin programs"), h1, r1, b1, [26, 34, 40] + [11, 11, 7] * ((len(h1) - 3) // 3))
    n_extra = ""
    if more:
        hg, rg, bg = sheet_cycle_genes()
        write_sheet(wb.create_sheet("Thioredoxin cycle genes"), hg, rg, bg, [12] + [11, 11, 7] * ((len(hg) - 1) // 3))
        hx, rx = sheet_txn_txnip()
        write_sheet(wb.create_sheet("TXN and TXNIP"), hx, rx, None, [18, 10, 10, 14, 8, 11, 12, 10])
        for row in wb["TXN and TXNIP"].iter_rows(min_row=2, min_col=4, max_col=4):
            row[0].number_format = "0.000"            # the difference in score to three decimals, as printed
        n_extra = (f"| Thioredoxin cycle genes | {len(rg)} genes | supptable2_cycle_genes.csv |\n"
                   f"| TXN and TXNIP | {len(rx)} cell types | supptable2_txn_txnip.csv |\n")
    n_aldehyde = ""
    if has_scan:
        h5, r5 = sheet_aldehyde_genes()
        write_sheet(wb.create_sheet("Aldehyde enzyme genes"), h5, r5, None, [12, 11, 8, 9, 13, 13, 8])
        for row in wb["Aldehyde enzyme genes"].iter_rows(min_row=2, min_col=4, max_col=6):
            for cell in row:
                cell.number_format = "0.000"          # effect, permutation mean and SD to three decimals, as printed
        n_aldehyde = f"| Aldehyde enzyme genes | {len(r5)} genes | fibroblast_genes_scan.csv |"
    else:
        print("fibroblast_genes_scan.csv not found: the sheet Aldehyde enzyme genes is not written "
              "(run 04_fibroblast_genes.py first)")
    h2, r2, order, by_id = sheet_programs()
    write_sheet(wb.create_sheet("Gene programs"), h2, r2, None,
                [24, 40, 40, 44, 8] + [13] * len(L.CELL_SETS) + [32])
    n_tnk = ""
    if has_tnk:
        h4, r4 = sheet_tnk_subtypes()
        write_sheet(wb.create_sheet("T and NK subtypes"), h4, r4, None,
                    [26, 34, 9, 8, 8, 8, 13, 8, 13, 13, 13, 16, 12, 16])
        for row in wb["T and NK subtypes"].iter_rows(min_row=2):
            for j in (5, 6, 7, 10, 11, 12):           # effect, SD and p to three decimals, as printed
                if isinstance(row[j].value, float):
                    row[j].number_format = "0.000"
        n_tnk = f"| T and NK subtypes | {len(r4)} rows | tnk_subtypes.csv |"
    else:
        print("tnk_subtypes.csv not found: the sheet T and NK subtypes is not written (run 12_tnk_subtypes.py first)")
    h3, r3, n_used, n_removed = sheet_genes(order, {pid: row[1] for pid, row in by_id.items()})
    write_sheet(wb.create_sheet("Gene lists"), h3, r3, None, [40, 14, 60] + [14] * len(L.CELL_SETS))

    path = target(f"{NAME}.xlsx")
    wb.save(path)
    lines = [f"# {NAME}.xlsx", "",
             "Supplementary Data 1 of the manuscript, written by "
             "chapter2/10_supplementary_data1.py from supptable2_programs.csv, supptable2b_programs.csv, "
             "supptable2b_genes.csv and supptable_gene_programs.csv. The sheet Read me describes every column.", "",
             "| sheet | rows | source |", "|---|---|---|",
             "| Read me | text | this script |",
             f"| Thioredoxin programs | {len(r1)} programs | supptable2_programs.csv (printed z and percentile, "
             "bold marks)" + (", with the cancer cells of supptable2_cancer.csv" if more else "") + " |",
             *([n_extra.rstrip("\n")] if more else []),
             *([n_aldehyde] if n_aldehyde else []),
             f"| Gene programs | {len(r2)} programs | supptable2b_programs.csv and the where_in_source column of "
             "supptable_gene_programs.csv, by the rules of the printed table (label, source, reference, drawn in) |",
             *([n_tnk] if n_tnk else []),
             f"| Gene lists | {len(r3)} genes: {n_used} used, {n_removed} removed by a rule | supptable2b_genes.csv "
             "and the edits and symbols_renamed columns of supptable_gene_programs.csv |", ""]
    path.with_suffix(".md").write_text("\n".join(lines))
    print(f"wrote {path} (Thioredoxin programs {len(r1)} rows, Gene programs {len(r2)}, Gene lists {len(r3)}: "
          f"{n_used} used, {n_removed} removed) in {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
