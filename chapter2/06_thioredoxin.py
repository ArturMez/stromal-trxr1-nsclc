#!/usr/bin/env python3
"""
Chapter 2, step 06: the thioredoxin system, the NRF2 program and candidate upstream signals under the
TXNRD1 split, calibrated against genes detected as often as TXNRD1.

ITEMS: Figure 2c (s3: every set as its TXNRD1 effect against the range of the comparison genes, in
      fibroblasts, macrophages, T and NK cells and cancer cells); the numbers of the sheet Thioredoxin programs
      of Supplementary Data 1 (built from trx_results.csv by 10_supplementary_data1.py) and of the thioredoxin
      section, "TXNRD1 expression in stromal cells is linked to the thioredoxin cycle, but not to hypoxia,
      inflammation or smoking"; Supplementary Figure 2a (the four genes of the cycle one by one, --split-cycle
      and s4) and Supplementary Figure 2c (cancer cells by histology, --cells can --histology and s4).
      Figure 2b is a schematic with no data and is not drawn here.

WHAT: s1, per cell set (fibroblasts, macrophages, T and NK cells; and cancer cells with --cells can), for
          the 26 sets of the sheet Thioredoxin programs of Supplementary Data 1 (modules of hand-made
          gene lists, programs of the package gmt and single genes; TXNRD1 removed from every list;
          the other TrxR1 substrates TXNDC17 and TXNL1 and the group Selenium status and ferroptosis are
          among them):
        - effect: mean paired difference in the set score (mean within-batch z of the count fractions
          of its tested genes), TXNRD1-detecting minus a non-detecting cell of the same batch and
          histology matched on genes detected (the pairing of 03, SEED = 6);
        - z, normal p and empirical p from 1,000 permutations of TXNRD1 detection within patient (the
          detection labels are shuffled among the cells of each patient) with the matcher re-run in each
          (the permutation stream of 03's REAL arm);
        - leave one gene out for every module of two or more tested genes (z of each leave-one-out
          score against the same permutations);
        - the comparison-gene rule: the same cells split by each gene detected within 20 % of TXNRD1's
          detection rate (the window is widened in steps of 5 % until it holds at least 200 genes; it
          never needed widening here), the pairing repeated with a fresh generator (SEED = 6) for each;
          the percentile of the TXNRD1 effect among these effects, the genes that are tested members of
          the set left out for that set. A set is unusual when TXNRD1's effect is at or above the 97.5th
          percentile (or at or below the 2.5th) of the comparison genes; it is called robust when
          |z| >= 1.96 and the effect holds with each gene left out in turn;
        - a module needs 3 tested genes (the package floor); one with fewer is not testable (n.t.):
          its effect and percentile are left empty, its z is written but not shown. Single genes are
          scored alone; one that fails the package detection rule (measured in fewer than 60 % or
          detected in fewer than 2 % of the cells) is not testable, GPX1 in every cell set because the
          count tables of the largest studies have no GPX1 row;
        - former gene symbols are merged into the current ones before scoring (MERGE: SELK into
          SELENOK, SEPW1 into SELENOW, C11orf31 into SELENOH, ENSG00000233276 into GPX1; no study
          measures both names, which every run checks). No other set holds these genes.
      s2, the three cell sets together (trx_results.csv, with the drawing class of every entry) and the
          counts that the text and the legends quote (trx_numbers.csv); with --cells can run before, the
          cancer cells are added in a separate table (trx_results_with_extra.csv).
      s3, the panel of Figure 2c (file suppfig2c_programs_calibrated).
      s4, the additional tests: trx_cycle_genes.csv, trx_cancer_histology.csv, trx_additional_numbers.csv and
          the panels panel_trx_cycle_genes (Supplementary Figure 2a) and panel_trx_cancer_histology
          (Supplementary Figure 2c). The variant runs that feed s4 are
        --cells can      the cells the original authors labeled as cancer cells, a fourth column of Figure 2c
        --split-cycle    TXN, PRDX1, PRDX2 and SRXN1 scored one by one (s1 only, *_split tables)
        --histology      the cells of LUSC or of LUAD within the sample of the main run (s1 only)
FROM: work_root/intermediate/chapter2/work_<cell set>/ (02_working_sets.py) and
      work_root/output/chapter2/gene_programs.gmt (01_programs.py: the programs ifn_gamma,
      matrix.tgfb, cell_cycle, e2f_targets, myc_targets and the NRF2 target list).
TO:   work_root/output/chapter2/, every table with a sibling .md:
        s1  trx_sets.gmt                    the 26 sets as scored (written by the first cell set)
            trx_tests_<cell>.csv            one row per set: effect, z, p, leave one gene out, comparison genes
        s2  trx_results.csv, trx_numbers.csv
        s3  suppfig2c_programs_calibrated.png/.pdf   (Figure 2c)
        s4  trx_cycle_genes.csv, trx_cancer_histology.csv, trx_additional_numbers.csv,
            panel_trx_cycle_genes and panel_trx_cancer_histology (.png/.pdf)
      and the logs work_root/logs/chapter2/06_<cell>.log.
WHY:  the thioredoxin cycle, its NADPH supply, the NRF2 targets and the candidate upstream signals are tested
      with one permutation test and one comparison-gene calibration, on one gmt, so that the effect of
      TXNRD1 on a set is judged against genes detected as often as TXNRD1. Every set's score depends on its
      own genes only, so scoring the 26 sets together gives the values of separate runs. The leave-one-out
      effects and the effect of every comparison split are summarized in trx_tests_<cell>.csv and not
      written gene by gene.
RUN:  python3 chapter2/06_thioredoxin.py [--only s1 s2 s3 s4] [--cells fib mac tnk can] [--nperm 1000]
      [--split-cycle] [--histology LUSC|LUAD] [--leave-out-study PREFIX [PREFIX ...]]
      Without --only the steps s1 to s3 run; s4 is named explicitly, after the variant runs.
      --leave-out-study removes the cells of one study or several (batch-name prefixes) before the cap, as
      03_program_tests.py does, and writes every table, log and panel with the suffix _without_<prefix> beside the
      ordinary ones, which it leaves untouched (trx_sets.gmt is read, not rewritten). The run that backs the paper
      is --leave-out-study Jimenez_Bruno_2026 (the ablation study GSE325414; Methods sentence on it).
      Seeds: SEED = 6 (cap, observed pairing, permutation stream, comparison-gene splits); 1,000 permutations;
      cap 2,000 cells per batch and histology; no other random numbers (the panels have no jitter).
"""
import argparse
import math
import os
import sys
import time
import warnings
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402
from scipy.sparse import csr_matrix  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ch2lib as L  # noqa: E402
from redox import paths  # noqa: E402
from redox import plotting as RP  # noqa: E402  (chooses the Agg backend)
import matplotlib.pyplot as plt  # noqa: E402

SEED = L.SEED       # 6
CAP = 2000          # cells per batch x histology, as in 03
SHEET = "the sheet Thioredoxin programs of Supplementary Data 1"
MINCTRL = 200       # comparison-gene rule: widen the detection window until at least this many comparison genes
MODULE_FLOOR = 3    # a module needs this many tested genes (the package MINSET)
LOGS = paths.WORK_ROOT / "logs" / "chapter2"

# ---------------------------------------------------------------- the sets of Supplementary Data 1
# kind: "module" (hand-made list, needs 3 tested genes), "program" (a program of the package gmt, the
# same floor), "gene" (one gene, scored alone). genes: the list, or the package set id it is taken from.
# label: the row label of the sheet Thioredoxin programs of Supplementary Data 1; panel: the row label of
# Figure 2c; listed: the genes column of that sheet.
TARGET_SET = "nrf2_targets"
DIRECT_PARTNERS = ("TXN", "TXNRD1", "PRDX1")                      # left out of the other-targets set
SETS = [
    dict(set_id="M1_trx_prx_cycle", group="Thioredoxin system", kind="module",
         genes=["TXN", "PRDX1", "PRDX2", "SRXN1"], label="Thioredoxin cycle",
         # the panel label names SRXN1: it passes the detection rule in the cancer cells only, so the
         # score is TXN, PRDX1 and PRDX2 in the three stromal cell sets and all four genes in the cancer cells
         panel="TXN, PRDX1, PRDX2, SRXN1 (thioredoxin cycle)"),
    dict(set_id="G_TXNDC17", group="Thioredoxin system", kind="gene", genes=["TXNDC17"], label="TXNDC17 (TRP14)",
         panel="TXNDC17 (TRP14, TrxR1 substrate)"),
    dict(set_id="G_TXNL1", group="Thioredoxin system", kind="gene", genes=["TXNL1"], label="TXNL1 (TRP32)",
         panel="TXNL1 (TRP32, TrxR1 substrate)"),
    dict(set_id="M3_nadph_supply", group="Thioredoxin system", kind="module",
         genes=["G6PD", "PGD", "TKT", "TALDO1", "ME1", "IDH1"], label="NADPH supply", panel="NADPH supply"),
    dict(set_id="M2_mito_trx", group="Thioredoxin system", kind="module",
         genes=["TXN2", "TXNRD2", "PRDX3", "PRDX5"], label="Mitochondrial thioredoxin",
         panel="Mitochondrial thioredoxin"),
    dict(set_id="M4_selenocysteine", group="Thioredoxin system", kind="module",
         genes=["SEPHS2", "SEPSECS", "PSTK", "EEFSEC", "SECISBP2"], label="Selenocysteine incorporation",
         panel="Selenocysteine incorporation"),
    # panel label "TXNIP" without "(thioredoxin inhibitor)": TXN and TXNIP regulate each other (PMID 41736403)
    dict(set_id="G_TXNIP", group="Thioredoxin system", kind="gene", genes=["TXNIP"], label="TXNIP",
         panel="TXNIP"),
    dict(set_id="M6_selenium", group="Selenium status and ferroptosis", kind="module",
         genes=["SELENOW", "SELENOH", "SELENOK"], label="Selenium-sensitive selenoproteins",
         panel="Selenium-sensitive selenoproteins"),
    dict(set_id="G_GPX1", group="Selenium status and ferroptosis", kind="gene", genes=["GPX1"], label="GPX1",
         panel="GPX1"),
    dict(set_id="G_GPX4", group="Selenium status and ferroptosis", kind="gene", genes=["GPX4"], label="GPX4",
         panel="GPX4"),
    dict(set_id="G_PRDX6", group="Selenium status and ferroptosis", kind="gene", genes=["PRDX6"], label="PRDX6",
         panel="PRDX6"),
    dict(set_id="G_AIFM2", group="Selenium status and ferroptosis", kind="gene", genes=["AIFM2"],
         label="AIFM2 (FSP1)", panel="AIFM2 (FSP1)"),
    dict(set_id="REF_NRF2_OTHER", group="NRF2 program", kind="module",
         genes=("minus", TARGET_SET, DIRECT_PARTNERS), label="NRF2 targets other than TXN, PRDX1",
         panel="NRF2 targets other than TXN, PRDX1"),
    dict(set_id="G_HMOX1", group="NRF2 program", kind="gene", genes=["HMOX1"], label="HMOX1",
         panel="HMOX1"),
    dict(set_id="G_SQSTM1", group="NRF2 program", kind="gene", genes=["SQSTM1"], label="SQSTM1",
         panel="SQSTM1"),
    dict(set_id="M5_glutathione", group="NRF2 program", kind="module",
         genes=["GCLC", "GCLM", "GSR", "GSS", "SLC7A11"], label="Glutathione synthesis",
         panel="Glutathione synthesis"),
    dict(set_id="C5_hypoxia", group="Candidate upstream signals", kind="module",
         genes=["VEGFA", "ADM", "BNIP3", "PGK1", "SLC2A1", "CA9", "NDRG1", "P4HA1"], label="Hypoxia",
         panel="Hypoxia"),
    dict(set_id="C6a_nfkb_tnf", group="Candidate upstream signals", kind="module",
         genes=["NFKBIA", "TNFAIP3", "CXCL8", "IL1B"], label="NF-κB and TNF", panel="NF-κB and TNF"),
    dict(set_id="C6b_ifn_gamma", group="Candidate upstream signals", kind="program",
         genes=("package", "ifn_gamma"), label="Interferon gamma", panel="Interferon gamma",
         listed="Hallmark interferon gamma response"),
    dict(set_id="C6c_tgfb", group="Candidate upstream signals", kind="program",
         genes=("package", "matrix.tgfb"), label="TGF-β", panel="TGF-β", listed="Hallmark TGF-β signaling"),
    dict(set_id="C6d_il6_stat3", group="Candidate upstream signals", kind="module",
         genes=["SOCS3", "JUNB", "BCL3", "IL6"], label="IL-6 and STAT3", panel="IL-6 and STAT3"),
    dict(set_id="C1_aldehyde", group="Candidate upstream signals", kind="module",
         genes=["AKR1C1", "AKR1C2", "AKR1C3", "ALDH1A1", "ALDH3A1"], label="Aldehyde-metabolizing enzymes",
         panel="Aldehyde-metabolizing enzymes"),
    dict(set_id="C2_ahr", group="Candidate upstream signals", kind="module",
         genes=["CYP1A1", "CYP1B1", "AHRR", "TIPARP"], label="Aryl hydrocarbon receptor",
         panel="Aryl hydrocarbon receptor"),
    dict(set_id="C7a_cell_cycle", group="Proliferation", kind="program", genes=("package", "cell_cycle"),
         label="Cell cycle", panel="Cell cycle", listed="S-phase and G2/M genes of Tirosh and colleagues"),
    dict(set_id="C7b_e2f", group="Proliferation", kind="program", genes=("package", "e2f_targets"),
         label="E2F targets", panel="E2F targets", listed="Hallmark E2F targets"),
    dict(set_id="C7c_myc", group="Proliferation", kind="program", genes=("package", "myc_targets"),
         label="MYC targets", panel="MYC targets", listed="Hallmark MYC targets"),
]
GROUPS = list(dict.fromkeys(s["group"] for s in SETS))
SET_BY_ID = {s["set_id"]: s for s in SETS}
M1 = "M1_trx_prx_cycle"
# the four genes of M1 scored one by one (which gene carries the thioredoxin cycle),
# added only with --split-cycle, which writes *_split files and trx_sets_split.gmt and leaves every other file alone
SPLIT_CYCLE = [dict(set_id=f"G_{g}", group="Thioredoxin system", kind="gene", genes=[g], label=g, panel=g)
               for g in ("TXN", "PRDX1", "PRDX2", "SRXN1")]
NT_READING = f"not testable (fewer than {MODULE_FLOOR} genes detected)"
NT_GENE = (f"not testable (measured in fewer than {round(100 * L.PRESENCE)} % or detected in fewer than "
           f"{round(100 * L.DETFLOOR)} % of the cells)")
# Former symbols merged into the current ones before scoring, new name first. The atlas studies use
# one name or the other, never both (checked at every run), so the counts of the former name move into the current one
# and the current name counts as measured wherever the former one was. No other set holds either name.
MERGE = [("SELENOK", "SELK"), ("SELENOW", "SEPW1"), ("SELENOH", "C11orf31"), ("GPX1", "ENSG00000233276")]

# colors of the panels (blue = higher, red = lower; not the histology colors)
BLUE, RED, AXIS, BAND = "#2a78d6", "#e34948", "#c3c2b7", "#efeee9"
INK, INK2, MUTED, GRID = RP.INK, RP.INK2, RP.MUTED, RP.GRID


def universe(S):
    """nmeas, det, keepg and UNI as ch2lib.load_set computes them, for the counts now in S."""
    nmeas = np.bincount(S["acode"], minlength=S["PRES"].shape[0]).astype(float) @ S["PRES"]
    colnnz = np.diff(S["X"].indptr)
    det = np.where(nmeas > 0, colnnz / np.maximum(nmeas, 1), 0.0)
    keepg = np.where((nmeas >= L.PRESENCE * S["N"]) & (det >= L.DETFLOOR))[0]
    S.update(nmeas=nmeas, det=det, keepg=keepg, UNI=set(S["genes"][keepg]))


def merge_symbols(S, log):
    """Move the counts of each former symbol of MERGE into its current name and mark the current name measured
    wherever the former one was; stops when one study measures both names or one cell has counts under both."""
    genes, GI, X = S["genes"], S["GI"], S["X"]
    PRES = S["PRES"].copy()
    G = len(genes)
    cols = np.arange(G)
    for new, old in MERGE:
        if new not in GI or old not in GI:
            log(f"  merge {old} into {new}: not both in the gene list, nothing to merge")
            continue
        i, j = GI[new], GI[old]
        if np.any(PRES[:, i] & PRES[:, j]):
            raise SystemExit(f"{new} and {old} are both measured in one study: the merge would not be exact")
        both = (X[:, i].toarray().ravel() > 0) & (X[:, j].toarray().ravel() > 0)
        if both.any():
            raise SystemExit(f"{int(both.sum())} cells have counts under both {new} and {old}")
        log(f"  merge {old} into {new}: {X[:, j].nnz:,} cells with counts moved")
        cols[j] = i
        PRES[:, i] |= PRES[:, j]
        PRES[:, j] = False
    M = csr_matrix((np.ones(G, X.dtype), (np.arange(G), cols)), shape=(G, G))
    Xm = (X @ M).tocsc()
    Xm.sort_indices()
    if not np.allclose(np.asarray(Xm.sum(1)).ravel(), S["tot"], rtol=1e-6):
        raise SystemExit("the merge changed the total counts of some cells")
    S["X"], S["PRES"] = Xm, PRES
    universe(S)


def set_genes():
    """set id -> gene list, in the order the gmt holds them (package programs from gene_programs.gmt)."""
    pkg = L.read_gmt(L.out("gene_programs.gmt"))
    res = {}
    for s in SETS:
        g = s["genes"]
        if isinstance(g, tuple) and g[0] == "package":
            res[s["set_id"]] = list(pkg[g[1]][1])
        elif isinstance(g, tuple) and g[0] == "minus":
            res[s["set_id"]] = [x for x in pkg[g[1]][1] if x not in g[2]]
        else:
            res[s["set_id"]] = list(g)
    return res


def listed_text(s, genes):
    """The genes column of the sheet Thioredoxin programs: the list, or the name of the source program."""
    return s.get("listed") or ", ".join(genes)


def gmt_text():
    """The text of trx_sets.gmt: one tab-separated line per set (id, row label, genes)."""
    return "\n".join("\t".join([k, SET_BY_ID[k]["label"]] + g) for k, g in set_genes().items()) + "\n"


def read_gmt_unchanged():
    """A leave-out run reads trx_sets.gmt of the ordinary run instead of rewriting it, after checking that it holds
    exactly the sets this script would write."""
    path = L.out("trx_sets.gmt")
    if not path.exists() or path.read_text() != gmt_text():
        raise SystemExit("trx_sets.gmt is missing or differs from the sets of this script: run the ordinary s1 first")
    return path


def write_gmt(path=None):
    """Write the gmt of the sets (default trx_sets.gmt) with its sidecar; returns the path."""
    path = path or L.out("trx_sets.gmt")
    path.write_text(gmt_text())
    L.sidecar(path, f"The {len(SETS)} sets of {SHEET} as 06 scores them, one line per set, tab separated "
              "(package programs copied from gene_programs.gmt; TXNRD1 is removed at scoring).", [
                  ("field 1", "set id"), ("field 2", f"row label of {SHEET}"),
                  ("fields 3 on", "the genes of the list, in list order")])
    return path


# ====================================================================== s1
def s1_tests(cell, a):
    """Tests of the 26 sets in one cell set: observed pairing, within-patient permutations (SEED = 6), leave one gene
    out, and the comparison-gene percentile; writes trx_tests_<cell>.csv."""
    spec = L.cell_spec(cell)
    log = L.Log()
    if a.split_cycle:                    # the extra single-gene sets go to their own gmt
        gmt = write_gmt(L.out("trx_sets_split.gmt"))
    elif a.leave_out_study or a.histology:
        gmt = read_gmt_unchanged()
    else:
        gmt = write_gmt()
    S = L.load_set(cell, CAP, leave_out=a.leave_out_study, histology=a.histology)
    merged = []
    merge_symbols(S, merged.append)
    sfx = (L.leave_out_suffix(a.leave_out_study) + (f"_{a.histology}" if a.histology else "")
           + ("_split" if a.split_cycle else ""))
    P = L.score_programs(S, gmt, minset=1)
    PROG, names = P["PROG"], P["names"]
    N, TX, depth, strat, pidx = S["N"], S["TX"], S["depth"], S["strat"], S["pidx"]
    pat, studies = S["pat"], S["studies"]
    strata = L.Strata(strat)
    log("=" * 100)
    log(f"THIOREDOXIN SYSTEM AND UPSTREAM SIGNALS.  {spec['label']}.  nperm {a.nperm}, cap {CAP}, seed {SEED}"
        + (f", without {' and '.join(L.leave_out_prefixes(a.leave_out_study))}" if sfx else ""))
    log(f"  {N:,} cells, {len(np.unique(pat))} patients, {len(studies)} batches, TXNRD1 detected in {TX.mean():.4f}")
    for m in merged:
        log(m)
    for s in SETS:
        p = PROG[s["set_id"]]
        log(f"  {s['set_id']:20s} tested {len(p['tested'])}/{len(p['listed'])}   "
            f"not testable: {','.join(p['not_testable'])[:100]}")
    missing = [s["set_id"] for s in SETS if s["set_id"] not in names]
    if missing:
        log(f"  no tested gene at all, not scored: {missing}")

    L.check_matcher(depth, TX, strat)
    log("  match_fast returns the pairs of the verbatim matcher on 3 successive draws (checked)")
    scorer = L.Scorer.build(P, S)

    # ------------------------------------------------ observed pairing and permutations (03's REAL stream)
    aa, bb = L.match_fast(depth, TX, strata, np.random.default_rng(SEED))
    obs, _, _, cnt = scorer.base(aa, bb)
    lg_obs = scorer.logo(aa, bb)
    r = np.random.default_rng(SEED)
    nul, lg_nul = [], []
    for i in range(a.nperm):
        pa, pb = L.match_fast(depth, L.perm_expo(TX, pidx, r), strata, r)
        if len(pa) < L.MINPAIR:
            continue
        nul.append(scorer.base(pa, pb)[0])
        lg_nul.append(scorer.logo(pa, pb))
        if (i + 1) % 250 == 0:
            log(f"    replicate {i + 1} / {a.nperm}   [{log.elapsed()}]")
    nul = np.vstack(nul)
    mu, sd = np.nanmean(nul, 0), np.nanstd(nul, 0)
    with np.errstate(invalid="ignore", divide="ignore"):
        z = np.where(sd > 0, (obs - mu) / sd, np.nan)
    pemp = (np.sum(np.abs(nul - mu) >= np.abs(obs - mu), 0) + 1) / (len(nul) + 1)
    LG = np.vstack(lg_nul)
    with warnings.catch_warnings():       # one-gene sets have no leave-one-out value: all-empty columns
        warnings.simplefilter("ignore", RuntimeWarning)
        lmu, lsd = np.nanmean(LG, 0), np.nanstd(LG, 0)
    with np.errstate(invalid="ignore", divide="ignore"):
        lg_z = np.where(lsd > 0, (lg_obs - lmu) / lsd, np.nan)
    log(f"  pairs {len(aa):,}, {len(nul)} replicates   [{log.elapsed()}]")

    # ------------------------------------------------ comparison genes: control splitter genes
    # genes detected within 20 % (DET_TOL) of TXNRD1's detection rate, widened by 5 % until at least MINCTRL genes
    # (97.5th percentile of their effects is the upper threshold; the set's own tested genes are left out per set)
    X, keepg, det, genes, GI = S["X"], S["keepg"], S["det"], S["genes"], S["GI"]
    tdet = det[GI[L.SPLITTER]]
    tol = L.DET_TOL
    while True:
        pool = [g for g in keepg if abs(det[g] - tdet) <= tol * tdet and genes[g] != L.SPLITTER]
        if len(pool) >= MINCTRL or tol >= 1.0:
            break
        tol += 0.05
    log(f"  rule 6: {len(pool):,} control splitter genes within {tol:.0%} of TXNRD1 detection {tdet:.4f}")
    REF = np.full((len(pool), len(names)), np.nan)
    for i, g in enumerate(pool):
        s0, s1 = X.indptr[g], X.indptr[g + 1]
        ex = np.zeros(N, bool)
        ex[X.indices[s0:s1][X.data[s0:s1] > 0]] = True
        pa, pb = L.match_fast(depth, ex, strata, np.random.default_rng(SEED))
        if len(pa) >= L.MINPAIR:
            REF[i] = scorer.base(pa, pb)[0]
    pool_names = genes[np.array(pool, int)]
    log(f"  rule 6 done   [{log.elapsed()}]")

    # ------------------------------------------------ the table
    off = P["off"]
    npat = len(np.unique(pat))
    nstudies = len({L.study_of_batch(b) for b in studies})
    rows = []
    for s in SETS:
        nm = s["set_id"]
        p = PROG[nm]
        rec = dict(cell=cell, cell_type=spec["label"], set_id=nm, group=s["group"], kind=s["kind"],
                   label=s["label"], n_listed=len(p["listed"]), n_tested=len(p["tested"]),
                   tested=",".join(p["tested"]), not_testable=",".join(p["not_testable"]))
        testable = nm in names and (s["kind"] == "gene" or len(p["tested"]) >= MODULE_FLOOR)
        rec["testable"] = testable
        if nm not in names:
            rec["reading"] = NT_GENE if s["kind"] == "gene" else NT_READING
            rows.append(rec)
            continue
        k = names.index(nm)
        rec.update(effect=obs[k] if testable else np.nan, z=z[k], p_gauss=2 * stats.norm.sf(abs(z[k])),
                   p_emp=pemp[k], perm_sd=sd[k], n_pairs=int(cnt[k]), matched_pairs=len(aa), n_replicates=len(nul))
        if s["kind"] != "gene" and len(p["tested"]) >= 2:
            zs = [(g, lg_z[off[k] + j]) for j, g in enumerate(p["tested"])]
            if np.isfinite(z[k]) and abs(z[k]) >= L.Z95:
                signed = np.array([np.sign(obs[k]) * zz for _, zz in zs], float)
                j = int(np.nanargmin(signed))
                rec.update(logo_holds="yes" if np.all(signed >= L.Z95) else "no", logo_weakest_gene=zs[j][0],
                           logo_weakest_z=float(zs[j][1]), logo_genes_that_break_it=int(np.sum(~(signed >= L.Z95))))
        rec.update(n_controls_pool=len(pool), det_tolerance=tol)
        if not testable:
            rec["reading"] = NT_GENE if s["kind"] == "gene" else NT_READING
        else:
            keep = ~np.isin(pool_names, p["tested"])
            ref = REF[keep, k]
            ref = ref[np.isfinite(ref)]
            pct = 100 * float((ref < obs[k]).mean())
            rec.update(percentile=pct, n_controls=len(ref), ctrl_p2_5=np.percentile(ref, 2.5),
                       ctrl_median=np.median(ref), ctrl_p97_5=np.percentile(ref, 97.5),
                       reading=("unusual, higher" if pct >= 97.5 else "unusual, lower" if pct <= 2.5 else "ordinary"))
        rows.append(rec)
    T = pd.DataFrame(rows)
    for c, v in (("cells", N), ("patients", npat), ("batches", len(studies)), ("studies", nstudies)):
        T[c] = v
    o = L.out(f"trx_tests_{cell}{sfx}.csv")
    T.to_csv(o, index=False)
    L.sidecar(o, f"{spec['label']}: the TXNRD1 contrast for the sets of {SHEET} (cap {CAP} per batch x histology, "
              f"seed {SEED}), with leave one gene out and rule 6"
              + (f", without the cells of {' and '.join(L.leave_out_prefixes(a.leave_out_study))}." if sfx else "."),
              TEST_COLUMNS)
    with pd.option_context("display.width", 220):
        log(T[["set_id", "n_tested", "effect", "z", "percentile", "n_controls", "logo_holds", "reading"]].to_string(
            index=False, float_format=lambda v: f"{v:.4g}"))
    log(f"done   [{log.elapsed()}]")
    log.write(LOGS / f"06_{cell}{sfx}.log")


TEST_COLUMNS = [
    ("cell", "fib, mac or tnk"), ("cell_type", "cell set as printed"), ("set_id", "set"),
    ("group", f"row group of {SHEET}"), ("kind", "module (hand-made list), program (package gmt) or gene"),
    ("label", f"row label of {SHEET}"), ("n_listed", "genes in the list"),
    ("n_tested", "genes tested (package detection rule; TXNRD1 removed)"), ("tested", "the tested genes"),
    ("not_testable", "listed genes not tested"),
    ("testable", "False for a module with fewer than 3 tested genes or a gene that fails the detection rule (n.t.)"),
    ("effect", "mean paired difference in the set score, TXNRD1-detecting minus partner (empty when not testable)"),
    ("z", "(effect minus the permutation mean) / permutation SD; written for every scored set, not shown when not testable"),
    ("p_gauss", "two-sided normal p of z"), ("p_emp", "empirical permutation p (centered)"),
    ("perm_sd", "permutation SD"), ("n_pairs", "valid pairs of the set"), ("matched_pairs", "pairs of the observed pairing"),
    ("n_replicates", "permutations used"),
    ("logo_holds", "yes when |z| stays at or above 1.96 in the effect's direction with each gene left out (set when |z| >= 1.96)"),
    ("logo_weakest_gene", "the gene whose removal gives the weakest z"), ("logo_weakest_z", "that z"),
    ("logo_genes_that_break_it", "genes whose removal breaks it"),
    ("n_controls_pool", "control splitter genes (rule 6)"), ("det_tolerance", "their detection window around TXNRD1's rate"),
    ("percentile", "percentile of the effect among the control effects (tested members of the set left out)"),
    ("n_controls", "control genes used for the set"), ("ctrl_p2_5", "2.5th percentile of the control effects"),
    ("ctrl_median", "their median"), ("ctrl_p97_5", "their 97.5th percentile"),
    ("reading", "unusual, higher (>= 97.5th), unusual, lower (<= 2.5th), ordinary, or not testable"),
    ("cells", "cells after the cap"), ("patients", "patients"), ("batches", "sequencing batches"),
    ("studies", "studies: batch names without the platform part")]


# ====================================================================== s2
def state(r):
    """Drawing class of one entry: hi, hi_fragile, lo, lo_fragile, ord (ordinary) or nt (not testable). hi: percentile
    at or above 97.5 with |z| >= 1.96 and the effect holding with each gene left out; lo: the same at or below 2.5."""
    if not np.isfinite(r.percentile):
        return "nt"
    robust = abs(r.z) >= L.Z95 and str(r.logo_holds) != "no"
    if r.percentile >= 97.5:
        return "hi" if robust else "hi_fragile"
    if r.percentile <= 2.5:
        return "lo" if robust else "lo_fragile"
    return "ord"


def s2_combine(a):
    """Combine the cell sets: trx_results.csv (drawing class per entry) and trx_numbers.csv (counts of text and legends)."""
    log = L.Log()
    sfx = L.leave_out_suffix(a.leave_out_study)
    parts = [pd.read_csv(L.out(f"trx_tests_{c}{sfx}.csv"), float_precision="round_trip") for c in L.CELL_SETS]
    R = pd.concat(parts, ignore_index=True)
    if "percentile" not in R:
        R["percentile"] = np.nan
    R["state"] = [state(r) for r in R.itertuples()]
    genes = set_genes()
    R["listed"] = R.set_id.map(lambda k: listed_text(SET_BY_ID[k], genes[k]))
    o = L.out(f"trx_results{sfx}.csv")
    R.to_csv(o, index=False)
    L.sidecar(o, "The three cell sets of trx_tests_<cell>.csv together, with the drawing class of "
              f"Figure 2c (bold in {SHEET} when hi): the columns of trx_tests_<cell>.csv, then", TEST_COLUMNS + [
                  ("state", "hi: at or above the 97.5th percentile, |z| >= 1.96 and not resting on one gene (filled "
                            "blue, bold); hi_fragile: at or above the 97.5th but not robust (open blue); lo, lo_fragile: "
                            "the same at or below the 2.5th; ord: within the range; nt: not testable"),
                  ("listed", f"the genes column of {SHEET}")])
    # a cell set of ch2lib.EXTRA_CELL_SETS that was tested (06 --only s1 --cells can) is added in a
    # separate table, so that trx_results.csv and the numbers below stay those of the three cell sets
    extra = [c for c in L.EXTRA_CELL_SETS if L.out(f"trx_tests_{c}{sfx}.csv").exists()]
    if extra:
        X = pd.concat([pd.read_csv(L.out(f"trx_tests_{c}{sfx}.csv"), float_precision="round_trip") for c in extra],
                      ignore_index=True)
        if "percentile" not in X:
            X["percentile"] = np.nan
        X["state"] = [state(r) for r in X.itertuples()]
        X["listed"] = X.set_id.map(lambda k: listed_text(SET_BY_ID[k], genes[k]))
        ox = L.out(f"trx_results_with_extra{sfx}.csv")
        pd.concat([R, X], ignore_index=True).to_csv(ox, index=False)
        L.sidecar(ox, f"trx_results{sfx}.csv with the tested sets of ch2lib.EXTRA_CELL_SETS added ({', '.join(extra)}; "
                  "can: the cells the original authors labeled as cancer cells), same columns", [])
        log(f"  wrote {ox.name}: {len(R)} rows of the three cell sets and {len(X)} of {', '.join(extra)}")

    ix = R.set_index(["cell", "set_id"])
    items = []

    def add(item, value, meaning):
        """Append one counted item (key, value, what it counts)."""
        items.append((item, value, meaning))

    for c, spec in L.CELL_SETS.items():
        r = ix.loc[(c, M1)]
        add(f"matched_pairs_{c}", int(r.matched_pairs), f"matched pairs, {spec['label']}")
        add(f"controls_m1_{c}", int(r.n_controls), f"control genes behind the thioredoxin cycle percentile, {spec['label']}")
    studies = [int(ix.loc[(c, M1)].studies) for c in L.CELL_SETS]
    add("studies_min", min(studies), "studies per cell set, smallest")
    add("studies_max", max(studies), "studies per cell set, largest")
    add("n_replicates", int(R.n_replicates.min()), "permutations used, smallest over every set and cell set")
    unusual = R.reading.str.startswith("unusual")
    add("m1_z_high_cells", int((R[R.set_id == M1].z >= L.Z95).sum()),
        "cell sets whose thioredoxin cycle score is higher in TXNRD1-detecting cells at z >= 1.96 ('expressed more')")
    add("m1_unusual_higher_tnk", int(ix.loc[("tnk", M1)].reading == "unusual, higher"),
        "thioredoxin cycle at or above the 97.5th percentile, T and NK cells (1 yes; 'not specific to TXNRD1')")
    add("trx_other_unusual", int((unusual & R.set_id.isin(["M2_mito_trx", "M4_selenocysteine", "G_TXNIP"])).sum()),
        "mitochondrial thioredoxin, selenocysteine incorporation, TXNIP: entries beyond the 2.5th to 97.5th range")
    mac_targets = ["REF_NRF2_OTHER", "G_HMOX1", "G_SQSTM1", "M3_nadph_supply", "M5_glutathione"]
    sub = R[(R.cell == "mac") & R.set_id.isin(mac_targets)]
    add("mac_nrf2_min_percentile", float(sub.percentile.min()),
        "macrophages: smallest percentile of the other targets, HMOX1, SQSTM1, NADPH supply, glutathione "
        "('each difference larger than for 98 %')")
    add("mac_nrf2_hi", int((sub.state == "hi").sum()), "macrophages: of those five, entries drawn filled (robust)")
    fib_targets = ["REF_NRF2_OTHER", "G_HMOX1", "G_SQSTM1", "M5_glutathione"]
    add("fib_nrf2_unusual", int((unusual & (R.cell == "fib") & R.set_id.isin(fib_targets)).sum()),
        "fibroblasts: other targets, HMOX1, SQSTM1, glutathione beyond the range ('but not in fibroblasts')")
    prolif = R[(R.cell == "tnk") & R.set_id.isin(["C7a_cell_cycle", "C7b_e2f", "C7c_myc"])]
    add("tnk_prolif_z_high", int((prolif.z >= L.Z95).sum()), "T and NK cells: cell cycle, E2F, MYC with z >= 1.96")
    add("tnk_prolif_unusual", int(prolif.reading.str.startswith("unusual").sum()),
        "T and NK cells: of those three, beyond the range ('were within the range')")
    up = ["C5_hypoxia", "C6a_nfkb_tnf", "C6b_ifn_gamma", "C6c_tgfb", "C6d_il6_stat3", "C1_aldehyde"]
    add("upstream_unusual", int((unusual & R.set_id.isin(up)).sum()),
        "upstream programs other than the aryl hydrocarbon receptor, testable entries beyond the range ('none')")
    add("ahr_not_testable_cells", int((R.set_id == "C2_ahr").sum() - R[R.set_id == "C2_ahr"].testable.sum()),
        "cell sets in which the aryl hydrocarbon receptor targets were not testable")
    new = ["G_TXNDC17", "G_TXNL1", "M6_selenium", "G_GPX1", "G_GPX4", "G_PRDX6", "G_AIFM2"]
    add("added_unusual", int((unusual & R.set_id.isin(new)).sum()),
        "TXNDC17, TXNL1 and the group Selenium status and ferroptosis: entries beyond the 2.5th to "
        "97.5th range (printed 'The same was true for TRP14 (TXNDC17) and TRP32 (TXNL1) ... and for the ferroptosis "
        "regulators GPX4 and PRDX6' and 'FSP1 (AIFM2) ... was also higher in these macrophages'; 1, FSP1 in macrophages)")
    add("aifm2_mac_hi", int(ix.loc[("mac", "G_AIFM2")].state == "hi"),
        "AIFM2 (FSP1) in macrophages drawn filled: above the range, |z| >= 1.96 (1 yes)")
    add("gpx1_not_testable_cells", int((R.set_id == "G_GPX1").sum() - R[R.set_id == "G_GPX1"].testable.sum()),
        "cell sets in which GPX1 was not testable (the largest studies do not report it)")
    add("supptable2_bold", int((R.state == "hi").sum()), f"entries in bold in {SHEET} (program x cell set)")
    add("supptable2_nt", int((R.state == "nt").sum()), f"entries marked n.t. in {SHEET} (program x cell set)")
    N2 = pd.DataFrame(items, columns=["item", "value", "meaning"])
    o = L.out(f"trx_numbers{sfx}.csv")
    N2.to_csv(o, index=False)
    L.sidecar(o, "Counts quoted in the thioredoxin section and in the legends of Figure 2 and "
              f"{SHEET}.", [("item", "key"), ("value", "value"), ("meaning", "what it counts")])
    log(N2.to_string(index=False))
    log.write(LOGS / f"06_combine{sfx}.log")


# ====================================================================== s3
def s3_panels(a):
    """Figure 2c (file suppfig2c_programs_calibrated): every set as its effect (circle) against the range of the
    comparison genes (gray bar, tick at the median), one column per cell set."""
    RP.style()
    plt.rcParams["axes.edgecolor"] = AXIS          # the axis color of the panels
    sfx = L.leave_out_suffix(a.leave_out_study)
    R = pd.read_csv(L.out(f"trx_results{sfx}.csv"), float_precision="round_trip")
    # the cancer cells are a fourth column of Figure 2c, read from
    # trx_results_with_extra.csv when s2 wrote it (an ordinary run only; a leave-out run keeps three columns)
    if not sfx and L.out("trx_results_with_extra.csv").exists():
        R = pd.read_csv(L.out("trx_results_with_extra.csv"), float_precision="round_trip")
    IX = R.set_index(["cell", "set_id"])
    cells = list(L.CELL_SETS.items()) + [(c, s) for c, s in L.EXTRA_CELL_SETS.items() if c in set(R.cell)]
    four = len(cells) == 4

    # ---------------------------------------------------------------- c: every set against the controls
    rows = []
    for g in GROUPS:
        rows.append(("H", g))
        rows += [("R", s["set_id"], s["panel"]) for s in SETS if s["group"] == g]
    nrow = len(rows)
    # 4.75 in for the 23 rows (19 sets, 4 group headings) of the first layout; each further row adds the
    # same height, so the rows, margins, title and legend keep their size
    W, ROWH, TOP, BOT = 7.2, (0.865 - 0.17) * 4.75 / 23, 0.135 * 4.75, 0.17 * 4.75
    H = TOP + BOT + nrow * ROWH
    fig = plt.figure(figsize=(W, H))
    cy0, cy1 = BOT / H, 1 - TOP / H
    dy = (cy1 - cy0) / nrow
    cx0, cw, cgap = 0.355, 0.195, 0.02
    if four:                                       # four columns in the width of three
        cw = (0.98 - cx0 - 3 * cgap) / 4
    shown = R.dropna(subset=["effect", "ctrl_p2_5", "ctrl_p97_5"])
    big = float(shown[["effect", "ctrl_p2_5", "ctrl_p97_5"]].abs().max().max())
    XMAX = 0.25 if big < 0.25 else math.ceil(big * 20) / 20 + 0.05   # 0.25 as published; wider only if needed
    ticks = [-0.2, -0.1, 0, 0.1, 0.2]
    for i, (c, spec) in enumerate(cells):
        ax = fig.add_axes([cx0 + i * (cw + cgap), cy0, cw, cy1 - cy0])
        ax.set_xlim(-XMAX, XMAX)
        ax.set_ylim(nrow - 0.5, -0.5)
        ax.axvline(0, color=AXIS, lw=0.6, zorder=1)
        ax.set_yticks([])
        if four:
            ax.set_xticks([-0.2, 0, 0.2])
            ax.set_xticklabels(["−0.2", "0", "0.2"])
        else:
            ax.set_xticks(ticks)
            ax.set_xticklabels(["−0.2", "−0.1", "0", "0.1", "0.2"])
        ax.tick_params(labelsize=6, length=2)
        ax.spines[["top", "right", "left"]].set_visible(False)
        npair = int(IX.loc[(c, M1)].matched_pairs)
        ax.set_title(f"{spec['label']}\n{npair:,} pairs", fontsize=7, pad=4)
        for j, row in enumerate(rows):
            if row[0] == "H":
                ax.axhline(j + 0.45, color=GRID, lw=0.5, zorder=1)
                continue
            r = IX.loc[(c, row[1])]
            s = r.state
            if s == "nt":
                ax.text(0, j, "not testable", ha="center", va="center", fontsize=5.6, color=MUTED, style="italic", zorder=3)
                continue
            ax.barh(j, r.ctrl_p97_5 - r.ctrl_p2_5, left=r.ctrl_p2_5, height=0.5, color=BAND, edgecolor=AXIS, lw=0.5,
                    zorder=2)
            ax.plot([r.ctrl_median] * 2, [j - 0.25, j + 0.25], color=INK2, lw=1.0, zorder=3, solid_capstyle="butt")
            kw = dict(s=24, zorder=4, clip_on=False)
            if s == "hi":
                ax.scatter(r.effect, j, color=BLUE, edgecolor=BLUE, lw=0.8, **kw)
            elif s == "hi_fragile":
                ax.scatter(r.effect, j, facecolor="white", edgecolor=BLUE, lw=0.9, **kw)
            elif s in ("lo", "lo_fragile"):
                ax.scatter(r.effect, j, facecolor=RED if s == "lo" else "white", edgecolor=RED, lw=0.9, **kw)
            else:
                ax.scatter(r.effect, j, facecolor="white", edgecolor=MUTED, lw=0.8, **kw)
        if i == 1 and not four:
            ax.set_xlabel("difference in program score, TXNRD1-positive minus partner cells", fontsize=6.4, labelpad=2)
    if four:                                       # centered under the four columns
        fig.text(cx0 + (4 * cw + 3 * cgap) / 2, cy0 - 0.21 / H,
                 "difference in program score, TXNRD1-positive minus partner cells", fontsize=6.4, ha="center", va="top",
                 color=INK2)
    for j, row in enumerate(rows):
        yy = cy1 - (j + 0.5) * dy
        if row[0] == "H":
            fig.text(0.03, yy, row[1], fontsize=6.9, fontweight="bold", va="center")
        else:
            fig.text(cx0 - 0.008, yy, row[2], fontsize=6.6, va="center", ha="right")
    fig.text(0.03, 1 - 0.025 * 4.75 / H, "TXNRD1 difference against genes detected as often", fontsize=7.4, fontweight="bold", va="top")
    key = ((0.0, 0.72, "band", "range of 95 % of genes detected as often as TXNRD1, tick: their median"),
           (0.53, 0.72, "hi", "TXNRD1, higher than 97.5 % of them, robust"),
           (0.0, 0.2, "frag", "TXNRD1, higher than 97.5 % of them, not robust"),
           (0.53, 0.2, "ord", "TXNRD1, within their range (ordinary)"))
    lg_h, bh = 0.075 * 4.75, 0.17
    if four:                                       # the cancer column has an entry below the range: a third row
        key = ((0.0, 0.83, "band", key[0][3]), (0.0, 0.5, "hi", key[1][3]), (0.53, 0.5, "frag", key[2][3]),
               (0.0, 0.17, "lo", "TXNRD1, lower than 2.5 % of them, robust"), (0.53, 0.17, "ord", key[3][3]))
        lg_h, bh = 0.45, 0.17 * 0.075 * 4.75 / 0.45
    lg = fig.add_axes([0.03, 0.0, 0.94, lg_h / H])
    lg.set_xlim(0, 1)
    lg.set_ylim(0, 1)
    lg.axis("off")
    for (x, y, kind, txt) in key:
        if kind == "band":
            lg.add_patch(plt.Rectangle((x, y - bh), 0.026, 2 * bh, color=BAND, ec=AXIS, lw=0.5))
            lg.plot([x + 0.013] * 2, [y - bh, y + bh], color=INK2, lw=1.0)
        else:
            lg.scatter([x + 0.013], [y], s=22, color=BLUE if kind == "hi" else RED if kind == "lo" else "white",
                       edgecolor=BLUE if kind in ("hi", "frag") else RED if kind == "lo" else MUTED, lw=0.9)
        lg.text(x + 0.035, y, txt, va="center", fontsize=6.1, color=INK2)
    RP.save(fig, 2, f"suppfig2c_programs_calibrated{sfx}")
    print(f"wrote suppfig2c_programs_calibrated{sfx}.png and .pdf")


# ====================================================================== s4
KEY_TEXT = dict(band="range of 95 % of the comparison genes, tick: their median",
                hi="higher than 97.5 % of them, robust", frag="higher than 97.5 % of them, not robust",
                lo="lower than 2.5 % of them, robust", ord="within their range (ordinary)")


def strip_panel(R, columns, rows, title, name, xlabel, label_x=0.355):
    """One strip per column and one row per set, drawn as the panel of s3: the effect (circle) against the range of
    the control splits (gray bar, tick at their median). R: one row per (col, set_id) with effect, ctrl_p2_5,
    ctrl_median, ctrl_p97_5 and state; columns: (key, heading); rows: ("H", group) or ("R", set_id, label)."""
    IX = R.set_index(["col", "set_id"])
    nrow = len(rows)
    W, ROWH, TOP, BOT = 7.2, (0.865 - 0.17) * 4.75 / 23, 0.135 * 4.75, 0.17 * 4.75
    H = TOP + BOT + nrow * ROWH
    fig = plt.figure(figsize=(W, H))
    cy0, cy1 = BOT / H, 1 - TOP / H
    dy = (cy1 - cy0) / nrow
    cgap = 0.02
    cw = (0.98 - 0.355 - 3 * cgap) / 4             # the column width of the four-column panel of s3
    for i, (c, head) in enumerate(columns):
        ax = fig.add_axes([label_x + i * (cw + cgap), cy0, cw, cy1 - cy0])
        ax.set_xlim(-0.25, 0.25)
        ax.set_ylim(nrow - 0.5, -0.5)
        ax.axvline(0, color=AXIS, lw=0.6, zorder=1)
        ax.set_yticks([])
        ax.set_xticks([-0.2, 0, 0.2])
        ax.set_xticklabels(["−0.2", "0", "0.2"])
        ax.tick_params(labelsize=6, length=2)
        ax.spines[["top", "right", "left"]].set_visible(False)
        ax.set_title(head, fontsize=7, pad=4)
        for j, row in enumerate(rows):
            if row[0] == "H":
                ax.axhline(j + 0.45, color=GRID, lw=0.5, zorder=1)
                continue
            r = IX.loc[(c, row[1])]
            if r.state == "nt":
                ax.text(0, j, "not testable", ha="center", va="center", fontsize=5.6, color=MUTED, style="italic",
                        zorder=3)
                continue
            ax.barh(j, r.ctrl_p97_5 - r.ctrl_p2_5, left=r.ctrl_p2_5, height=0.5, color=BAND, edgecolor=AXIS, lw=0.5,
                    zorder=2)
            ax.plot([r.ctrl_median] * 2, [j - 0.25, j + 0.25], color=INK2, lw=1.0, zorder=3, solid_capstyle="butt")
            edge = BLUE if r.state.startswith("hi") else RED if r.state.startswith("lo") else MUTED
            face = edge if r.state in ("hi", "lo") else "white"
            ax.scatter(r.effect, j, facecolor=face, edgecolor=edge, lw=0.9, s=24, zorder=4, clip_on=False)
    fig.text(label_x + (len(columns) * cw + (len(columns) - 1) * cgap) / 2, cy0 - 0.21 / H, xlabel, fontsize=6.4,
             ha="center", va="top", color=INK2)
    for j, row in enumerate(rows):
        yy = cy1 - (j + 0.5) * dy
        if row[0] == "H":
            fig.text(0.03, yy, row[1], fontsize=6.9, fontweight="bold", va="center")
        else:
            fig.text(label_x - 0.008, yy, row[2], fontsize=6.6, va="center", ha="right")
    fig.text(0.03, 1 - 0.025 * 4.75 / H, title, fontsize=7.4, fontweight="bold", va="top")
    lg = fig.add_axes([0.03, 0.0, 0.94, 0.45 / H])
    lg.set_xlim(0, 1)
    lg.set_ylim(0, 1)
    lg.axis("off")
    bh = 0.17 * 0.075 * 4.75 / 0.45
    for (x, y, kind) in ((0.0, 0.83, "band"), (0.0, 0.5, "hi"), (0.53, 0.5, "frag"), (0.0, 0.17, "lo"),
                         (0.53, 0.17, "ord")):
        if kind == "band":
            lg.add_patch(plt.Rectangle((x, y - bh), 0.026, 2 * bh, color=BAND, ec=AXIS, lw=0.5))
            lg.plot([x + 0.013] * 2, [y - bh, y + bh], color=INK2, lw=1.0)
        else:
            lg.scatter([x + 0.013], [y], s=22, color=BLUE if kind == "hi" else RED if kind == "lo" else "white",
                       edgecolor=BLUE if kind in ("hi", "frag") else RED if kind == "lo" else MUTED, lw=0.9)
        lg.text(x + 0.035, y, KEY_TEXT[kind], va="center", fontsize=6.1, color=INK2)
    RP.save(fig, 2, name)
    print(f"wrote {name}.png and .pdf")


def s4_additional(a):
    """The two additional tests of Supplementary Figure 2 (2a: the cycle genes one by one; 2c: the cancer cells by
    histology), combined from the s1 tables of the variant runs and drawn as the panel of s3.
      trx_cycle_genes.csv       TXN, PRDX1, PRDX2 and SRXN1 one by one in every cell set (06 --split-cycle)
      trx_cancer_histology.csv  the 26 sets in the cancer cells of LUSC and of LUAD (06 --cells can --histology)
      trx_additional_numbers.csv  the counts the text quotes from the two
    The tables of the paper's main run (trx_results.csv, trx_numbers.csv, trx_sets.gmt) are not read or written."""
    RP.style()
    plt.rcParams["axes.edgecolor"] = AXIS
    log = L.Log()
    cells = list(L.CELL_SETS.items()) + [(c, s) for c, s in L.EXTRA_CELL_SETS.items()
                                         if L.out(f"trx_tests_{c}_split.csv").exists()]
    ids = [s["set_id"] for s in SPLIT_CYCLE]
    parts = []
    for c, spec in cells:
        T = pd.read_csv(L.out(f"trx_tests_{c}_split.csv"), float_precision="round_trip")
        main = pd.read_csv(L.out(f"trx_tests_{c}.csv"), float_precision="round_trip")
        same = T[T.set_id.isin(main.set_id)].reset_index(drop=True)
        if not same.equals(main):          # the single genes are extra sets: the 26 sets must be those of the main run
            raise SystemExit(f"trx_tests_{c}_split.csv: the 26 sets differ from trx_tests_{c}.csv; rerun --split-cycle")
        parts.append(T[T.set_id.isin(ids)])
    G = pd.concat(parts, ignore_index=True)
    G["state"] = [state(r) for r in G.itertuples()]
    o = L.out("trx_cycle_genes.csv")
    G.to_csv(o, index=False)
    L.sidecar(o, "The four genes of the thioredoxin cycle scored one by one under the TXNRD1 split, per cell set "
              "(06 --split-cycle; Supplementary Figure 2a): the columns of trx_tests_<cell>.csv, then",
              TEST_COLUMNS + [("state", "hi: at or above the 97.5th percentile and |z| >= 1.96; hi_fragile: at or "
                                        "above the 97.5th with |z| below 1.96; lo, lo_fragile: the same at or below "
                                        "the 2.5th; ord: within the range; nt: not testable")])
    hists = [h for h in L.HISTS if L.out(f"trx_tests_can_{h}.csv").exists()]
    parts = []
    for h in hists:
        T = pd.read_csv(L.out(f"trx_tests_can_{h}.csv"), float_precision="round_trip")
        T.insert(1, "histology", h)
        parts.append(T)
    Hh = pd.concat(parts, ignore_index=True)
    if "percentile" not in Hh:
        Hh["percentile"] = np.nan
    Hh["state"] = [state(r) for r in Hh.itertuples()]
    o = L.out("trx_cancer_histology.csv")
    Hh.to_csv(o, index=False)
    L.sidecar(o, "The sets of trx_tests_can.csv in the cancer cells of one histology (06 --cells can --histology; "
              "Supplementary Figure 2c). The cells are those of the ordinary run with that histology, so the two "
              "histologies together are its cells and patients. Columns: histology, the columns of "
              "trx_tests_<cell>.csv, then", [("state", "drawing class, as in trx_results.csv")])
    full = pd.read_csv(L.out("trx_tests_can.csv"), float_precision="round_trip").iloc[0]
    one = Hh.groupby("histology").first()
    if int(one.cells.sum()) != int(full.cells) or int(one.patients.sum()) != int(full.patients):
        raise SystemExit("the cancer cells of the two histologies are not the cells of the ordinary run")
    log(f"  cycle genes: {len(G)} rows ({', '.join(c for c, _ in cells)}); cancer cells by histology: "
        + ", ".join(f"{h} {int(one.loc[h].cells):,} cells, {int(one.loc[h].patients)} patients, "
                    f"{int(one.loc[h].matched_pairs):,} pairs" for h in hists))
    with pd.option_context("display.width", 220):
        log(G[["cell", "set_id", "effect", "z", "percentile", "state"]].to_string(
            index=False, float_format=lambda v: f"{v:.4g}"))
        log(Hh[["histology", "set_id", "effect", "z", "percentile", "state"]].to_string(
            index=False, float_format=lambda v: f"{v:.4g}"))
    # the counts that the text and the legend of Supplementary Figure 2 quote
    gx, hx = G.set_index(["cell", "set_id"]), Hh.set_index(["histology", "set_id"])
    core = ["REF_NRF2_OTHER", "M3_nadph_supply", "M5_glutathione", "C1_aldehyde"]
    items = [
        ("prdx1_fib_hi", int(gx.loc[("fib", "G_PRDX1")].state == "hi"),
         "PRDX1 alone in fibroblasts: at or above the 97.5th percentile with |z| >= 1.96 (1 yes)"),
        ("txn_mac_hi", int(gx.loc[("mac", "G_TXN")].state == "hi"),
         "TXN alone in macrophages: at or above the 97.5th percentile with |z| >= 1.96 (1 yes)"),
        ("txn_fib_hi", int(gx.loc[("fib", "G_TXN")].state == "hi"),
         "TXN alone in fibroblasts: the same (0 no; 'just below this threshold')"),
        ("prdx1_mac_hi", int(gx.loc[("mac", "G_PRDX1")].state == "hi"),
         "PRDX1 alone in macrophages: the same (0 no; 'just below this threshold')"),
        ("prdx2_fib_mac_beyond", int(sum(gx.loc[(c, "G_PRDX2")].state != "ord" for c in ("fib", "mac"))),
         "PRDX2 alone, fibroblasts and macrophages: entries beyond the 2.5th to 97.5th range (0)"),
        ("single_genes_tnk_beyond", int(sum(gx.loc[("tnk", k)].state not in ("ord", "nt") for k in ids)),
         "T and NK cells: single genes of the cycle beyond the range (0)"),
        ("srxn1_testable_cells", int(sum(gx.loc[(c, "G_SRXN1")].state != "nt" for c, _ in cells)),
         "cell sets in which SRXN1 was testable (1, the cancer cells)"),
        ("can_txn_prdx1_hi", int(sum(gx.loc[("can", k)].state == "hi" for k in ("G_TXN", "G_PRDX1"))) if "can" in dict(cells) else 0,
         "cancer cells: TXN and PRDX1 alone at or above the 97.5th percentile with |z| >= 1.96 (2; 'through TXN and PRDX1')"),
        ("can_core_hi_both", int(sum(all(hx.loc[(h, k)].state == "hi" for h in hists) for k in core)),
         "cancer cells by histology: of the other NRF2 targets, NADPH supply, glutathione synthesis and the "
         "aldehyde-metabolizing enzymes, sets drawn filled in both LUSC and LUAD (4; 'higher in both')"),
        ("can_txnip_lo_LUSC", int(hx.loc[("LUSC", "G_TXNIP")].state == "lo") if "LUSC" in hists else 0,
         "cancer cells of LUSC: TXNIP at or below the 2.5th percentile with |z| >= 1.96 (0 no)"),
        ("can_txnip_lo_LUAD", int(hx.loc[("LUAD", "G_TXNIP")].state == "lo") if "LUAD" in hists else 0,
         "cancer cells of LUAD: the same (1 yes; 'more than for 97.5 % of the comparison genes only in the larger LUAD "
         "group')"),
        # TXNIP in the cancer cells of each histology: lower in both, by the permutation z
        ("can_txnip_z_low_both", int(sum(hx.loc[(h, "G_TXNIP")].z <= -L.Z95 for h in hists)),
         "cancer cells: histologies in which TXNIP is lower in TXNRD1-positive cells at z <= -1.96 (2; 'TXNIP was lower "
         "in both')"),
    ]
    N4 = pd.DataFrame(items, columns=["item", "value", "meaning"])
    o = L.out("trx_additional_numbers.csv")
    N4.to_csv(o, index=False)
    L.sidecar(o, "Counts quoted with the additional tests of Supplementary Figure 2a and 2c (06 s4).",
              [("item", "key"), ("value", "value"), ("meaning", "what it counts")])
    log(N4.to_string(index=False))
    xlabel = "difference in program score, TXNRD1-positive minus partner cells"
    G["col"] = G.cell
    npair = {c: int(G[G.cell == c].matched_pairs.max()) for c, _ in cells}
    strip_panel(G, [(c, f"{spec['label']}\n{npair[c]:,} pairs") for c, spec in cells],
                [("H", "Thioredoxin cycle, gene by gene")] + [("R", s["set_id"], s["panel"]) for s in SPLIT_CYCLE],
                "TXNRD1 difference for each gene of the thioredoxin cycle", "panel_trx_cycle_genes", xlabel)
    Hh["col"] = Hh.histology
    rows = []
    for g in GROUPS:
        rows.append(("H", g))
        rows += [("R", s["set_id"], s["panel"]) for s in SETS if s["group"] == g]
    strip_panel(Hh, [(h, f"Cancer cells, {h}\n{int(one.loc[h].matched_pairs):,} pairs") for h in hists], rows,
                "TXNRD1 difference in the cancer cells of each histology", "panel_trx_cancer_histology", xlabel)
    log.write(LOGS / "06_additional.log")


def main():
    """Run the steps named with --only (default s1 to s3), over the cell sets named with --cells."""
    ap = argparse.ArgumentParser(description=__doc__.split("RUN:")[0].strip(),
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", nargs="*", default=["s1", "s2", "s3"], choices=["s1", "s2", "s3", "s4"])
    ap.add_argument("--cells", nargs="*", default=list(L.CELL_SETS),
                    choices=list(L.CELL_SETS) + list(L.EXTRA_CELL_SETS))
    ap.add_argument("--nperm", type=int, default=1000)
    ap.add_argument("--split-cycle", action="store_true",
                    help="also score TXN, PRDX1, PRDX2 and SRXN1 one by one (writes *_split files, s1 only)")
    ap.add_argument("--histology", choices=list(L.HISTS), default=None,
                    help="the cells of one histology only (writes *_<histology> files, s1 only)")
    ap.add_argument("--leave-out-study", nargs="+", default=None, metavar="PREFIX",
                    help="batch-name prefix of a study to leave out, or several (writes *_without_<prefix> files)")
    a = ap.parse_args()
    if a.split_cycle:
        SETS.extend(SPLIT_CYCLE)
        SET_BY_ID.update({s["set_id"]: s for s in SPLIT_CYCLE})
    if a.split_cycle or a.histology:     # these variants write their s1 tables only; s2 and s3 stay those of the paper
        a.only = [x for x in a.only if x == "s1"]
    LOGS.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    if "s1" in a.only:
        for cell in a.cells:
            s1_tests(cell, a)
    if "s2" in a.only:
        s2_combine(a)
    if "s3" in a.only:
        s3_panels(a)
    if "s4" in a.only:
        s4_additional(a)
    print(f"06_thioredoxin.py done in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
