#!/usr/bin/env python3
"""
Chapter 2, step 12: the four programs of T and NK cells in Figure 3, tested within the T- and NK-cell subtypes.

PAPER: Supplementary Data 1, sheet "T and NK subtypes", and the Results sentence that the tumor-infiltrating Treg
       signature difference was found in conventional CD4 T cells and in regulatory T cells, but not in CD8 T cells
       or NK cells (gene-program section). The Methods and Supplementary Methods describe these tests.

WHAT: the four programs drawn for T and NK cells in Figure 3 (MYC target genes, mTORC1 signaling genes, cell-cycle
      genes, tumor-infiltrating Treg signature), each tested
        - in all T and NK cells, the Figure 3 run (the first block; checked against fig2_rows.csv of 03), and
        - within conventional CD4 T cells, regulatory T cells, CD8 T cells and NK cells: the capped cells of the
          T and NK working set restricted to one subtype (the level-2 label of the atlas), with the pipeline of 03
          re-run on them (gene universe, within-batch z, pairs, permutations, the library-size control), so every
          pair joins two cells of the same subtype, batch and histology.
      Both blocks are repeated without the cells of Guo and colleagues (PMID 29942094; batch prefix Guo_Zhang_2018),
      with the cap redrawn as in the leave-out runs of 03.
      Cells without one of the four subtype labels are not in any subtype: the cells labeled as other T cells and
      all cells of the one sequencing batch that has no subtype labels (E-MTAB-6653).
      Statistics, per program and cell set, as in 03 s1:
        - effect: the mean paired difference in the program score (mean within-batch z of the count fraction of
          the tested genes, TXNRD1 removed), TXNRD1-detecting cell minus a non-detecting cell of the same batch and
          histology matched on genes detected (match_fast, seed 6);
        - permutation SD, z and p: 1,000 permutations of TXNRD1 detection among the cells of each patient, the
          matcher re-run in each, all drawn from one generator seeded 6; z = (effect minus the permutation mean) /
          permutation SD; p = two-sided normal p of z (as printed in Figure 3), not corrected for multiple tests;
          p_emp, the empirical p, is written as well;
        - batches in the same direction: the batches whose mean paired difference has the sign of the effect,
          among the batches with at least 20 valid pairs;
        - control A, the library-size split at TXNRD1's prevalence (control_a of ch2lib), with its own generator
          seeded 6, as in 03.
FROM: work_root/intermediate/chapter2/work_tnk/ (02_working_sets.py), work_root/output/chapter2/gene_programs.gmt
      (01_programs.py) and fig2_rows.csv (03_program_tests.py; the program labels and the check of the first block).
TO:   work_root/output/chapter2/tnk_subtypes.csv (one row per run, cell set and program) and
      tnk_subtypes_sets.gmt (the four programs, copied line for line from gene_programs.gmt), each with a sibling
      .md, and the log work_root/logs/chapter2/12_tnk_subtypes.log. 10_supplementary_data1.py writes the sheet.
RUN:  python3 chapter2/12_tnk_subtypes.py [--nperm 1000]
      Run after 03_program_tests.py. Seeds: SEED = 6 (cap, observed pairing, permutations, control A), as in 03.
      About 5 minutes on an Apple M2.
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.sparse import csc_matrix

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ch2lib as L  # noqa: E402
from redox import paths  # noqa: E402

SEED = L.SEED                       # 6, the seed of 03
CAP = 2000                          # cells per batch x histology, as in 03
LEAVE = "Guo_Zhang_2018"            # the study of PMID 29942094
PROGRAMS = ["myc_targets", "mtorc1", "cell_cycle", "treg"]     # the T and NK cell rows of Figure 3, in their order
SUBTYPES = [("CD4_T", "Conventional CD4 T cells"), ("Treg", "Regulatory T cells"), ("CD8_T", "CD8 T cells"),
            ("NK", "NK cells")]
ALL = ("all", "T and NK cells (Figure 3)")
LOGS = paths.WORK_ROOT / "logs" / "chapter2"
COLUMNS = [
    ("run", "all cells, or without_Guo_Zhang_2018 (the cells of PMID 29942094 left out)"),
    ("subtype", "all, CD4_T, Treg, CD8_T or NK (level-2 label of the atlas)"),
    ("subtype_label", "the cell set as printed"), ("program", "set id of gene_programs.gmt"),
    ("program_label", "the program as printed in Figure 3"), ("cells", "cells of the set after the cap"),
    ("txnrd1_positive", "cells with TXNRD1 detected"), ("patients", "patients"),
    ("batches", "sequencing batches"), ("pairs", "pairs of the observed pairing"),
    ("n_tested_genes", "tested genes of the program in this set"),
    ("effect", "mean paired difference in the program score, TXNRD1-positive minus partner"),
    ("perm_sd", "SD of the effect over the permutations"),
    ("z", "(effect minus the permutation mean) / permutation SD"),
    ("p_gauss", "two-sided normal p of z (as printed)"), ("p_emp", "empirical permutation p (centered)"),
    ("n_replicates", "permutations used"), ("valid_pairs", "pairs in which both cells have a score"),
    ("batches_same_direction", "batches whose mean paired difference has the sign of the effect"),
    ("batches_with_20_pairs", "batches with at least 20 valid pairs"),
    ("controlA_effect", "effect of the library-size split (control A)"), ("controlA_z", "its permutation z")]


def four_program_gmt():
    """The four programs, copied line for line and in the order of gene_programs.gmt. A program's score depends
    only on its own genes, so the scores equal those of the full gmt; the output rows follow PROGRAMS."""
    lines = [ln for ln in open(L.out("gene_programs.gmt")) if ln.split("\t")[0] in PROGRAMS]
    if sorted(ln.split("\t")[0] for ln in lines) != sorted(PROGRAMS):
        raise SystemExit("gene_programs.gmt does not hold each of the four programs exactly once")
    path = L.out("tnk_subtypes_sets.gmt")
    path.write_text("".join(lines))
    L.sidecar(path, "The four programs of T and NK cells in Figure 3 that 12_tnk_subtypes.py tests, copied from "
              "gene_programs.gmt (fields: set id, source, genes).", [])
    return path


def subtype_labels(S):
    """The level-2 label of every cell of a working set, in the row order of load_set ('' when the study gave none)."""
    mA, mB = L.load_meta(S["W"])
    l2 = np.concatenate([mA["l2"].astype(str), mB["l2"].astype(str)])
    return l2[S["sel"]]


def subset_set(S, keep):
    """The working set restricted to the cells in keep, with every derived quantity (gene universe, arm presence,
    detection rule, batches, patients) recomputed from those cells as load_set computes it for a whole set."""
    idx = np.where(keep)[0]
    X = csc_matrix(S["X"][idx])
    block, pat, arm = S["block"][idx], S["pat"][idx], S["arm"][idx]
    depth, hist, strat = S["depth"][idx], S["hist"][idx], S["strat"][idx]
    genes, GI = S["genes"], S["GI"]
    present = pd.read_pickle(S["W"] / "present.pkl")
    N = X.shape[0]
    tot = np.asarray(X.sum(1)).ravel().astype(np.float32)
    TX = np.asarray(X[:, GI[L.SPLITTER]].todense()).ravel() > 0
    arms_u = np.array(sorted(set(arm)))
    acode = np.searchsorted(arms_u, arm)
    PRES = np.zeros((len(arms_u), len(genes)), bool)
    for i, a in enumerate(arms_u):
        if a in present.columns:
            PRES[i] = present[a].values
    nmeas = np.bincount(acode, minlength=len(arms_u)).astype(float) @ PRES
    colnnz = np.diff(X.indptr)
    det = np.where(nmeas > 0, colnnz / np.maximum(nmeas, 1), 0.0)
    # a gene is tested when measured in at least 60 % of the cells and detected in at least 2 % of those measured
    keepg = np.where((nmeas >= L.PRESENCE * N) & (det >= L.DETFLOOR))[0]
    blks = np.unique(block)
    blk_of = np.searchsorted(blks, block)
    blk_arm = np.array([np.unique(acode[blk_of == b])[0] for b in range(len(blks))])
    pidx = [np.where(pat == p)[0] for p in np.unique(pat)]
    return dict(S, X=X, N=N, block=block, pat=pat, arm=arm, depth=depth, hist=hist, strat=strat, tot=tot,
                TX=TX, PRES=PRES, acode=acode, nmeas=nmeas, det=det, keepg=keepg, UNI=set(genes[keepg]),
                blks=blks, blk_of=blk_of, blk_arm=blk_arm, pidx=pidx, studies=blks, study_of=blk_of)


def run_arm(scorer, e, S, nperm):
    """The observed pairing and nperm permuted pairings of one split e, all from one generator seeded 6 (the
    stream of 03 s1). A permuted pairing with fewer than 50 pairs is skipped, as in 03."""
    strata = L.Strata(S["strat"])
    aa, bb = L.match_fast(S["depth"], e, strata, np.random.default_rng(SEED))
    obs, D, V, cnt = scorer.base(aa, bb)
    r = np.random.default_rng(SEED)
    nul = []
    for _ in range(nperm):
        pa, pb = L.match_fast(S["depth"], L.perm_expo(e, S["pidx"], r), strata, r)
        if len(pa) >= L.MINPAIR:
            nul.append(scorer.base(pa, pb)[0])
    return dict(aa=aa, obs=obs, D=D, V=V, cnt=cnt, nul=np.vstack(nul))


def summarize(res, names, block):
    """Per program: effect, permutation SD, z, normal and empirical p, valid pairs and the batch vote (batches
    with at least 20 valid pairs whose mean paired difference has the sign of the effect), as 03 s1."""
    obs, nul = res["obs"], res["nul"]
    mu, sd = np.nanmean(nul, 0), np.nanstd(nul, 0)
    with np.errstate(invalid="ignore", divide="ignore"):
        z = np.where(sd > 0, (obs - mu) / sd, np.nan)
    pemp = (np.sum(np.abs(nul - mu) >= np.abs(obs - mu), 0) + 1) / (len(nul) + 1)
    ba = block[res["aa"]]
    out = {}
    for k, nm in enumerate(names):
        v = res["V"][:, k]
        d, b_ = res["D"][v, k], ba[v]
        bm = np.array([d[b_ == u].mean() for u in np.unique(b_) if (b_ == u).sum() >= L.MINBLKPAIR])
        nsame = int((bm > 0).sum()) if obs[k] >= 0 else int((bm < 0).sum())
        out[nm] = dict(effect=obs[k], perm_sd=sd[k], z=z[k],
                       p_gauss=2 * stats.norm.sf(abs(z[k])) if np.isfinite(z[k]) else np.nan, p_emp=pemp[k],
                       n_replicates=len(nul), valid_pairs=int(res["cnt"][k]), batches_same_direction=nsame,
                       batches_with_20_pairs=len(bm))
    return out


def test_set(S, gmt, run, subtype, label, labels, nperm, log):
    """The pipeline of 03 s1 on one cell set: scores, observed and permuted pairings, control A. One row per
    program of PROGRAMS."""
    P = L.score_programs(S, gmt, log=lambda *a, **k: None)
    names = P["names"]
    L.check_matcher(S["depth"], S["TX"], S["strat"])
    scorer = L.Scorer.build(P, S)
    real = run_arm(scorer, S["TX"], S, nperm)
    R = summarize(real, names, S["block"])
    C = summarize(run_arm(scorer, L.control_a(S), S, nperm), names, S["block"])
    rows = []
    for nm in PROGRAMS:
        if nm not in R:
            raise SystemExit(f"{label}: {nm} has fewer than {L.MINSET} tested genes")
        rows.append(dict(run=run, subtype=subtype, subtype_label=label, program=nm, program_label=labels[nm],
                         cells=S["N"], txnrd1_positive=int(S["TX"].sum()), patients=len(np.unique(S["pat"])),
                         batches=len(np.unique(S["block"])), pairs=len(real["aa"]), n_tested_genes=len(P["PROG"][nm]["tested"]), **R[nm],
                         controlA_effect=C[nm]["effect"], controlA_z=C[nm]["z"]))
    log(f"  {run:26s} {label:27s} cells {S['N']:6,d}  TXNRD1+ {int(S['TX'].sum()):5,d}  "
        f"patients {rows[0]['patients']}  batches {rows[0]['batches']}   [{log.elapsed()}]")
    for rw in rows:
        log(f"      {rw['program']:12s} effect {rw['effect']:+.4f}  SD {rw['perm_sd']:.4f}  z {rw['z']:+.2f}  "
            f"p {rw['p_gauss']:.3g}  {rw['batches_same_direction']} of {rw['batches_with_20_pairs']}  "
            f"valid pairs {rw['valid_pairs']}")
    return rows


def check_first_block(rows):
    """The first block (all T and NK cells, all studies) must equal the T and NK rows of Figure 3 (fig2_rows.csv
    of 03) exactly: effect, z, p and the batch vote."""
    F = pd.read_csv(L.out("fig2_rows.csv"), float_precision="round_trip")
    F = F[F.cell_set == "tnk"].set_index("set_id")
    for rw in rows:
        f = F.loc[rw["program"]]
        same = (rw["effect"] == f.effect and rw["z"] == f.z and rw["p_gauss"] == f.p_gauss
                and rw["batches_same_direction"] == f.batches_same_direction
                and rw["batches_with_20_pairs"] == f.batches_with_20_pairs and rw["pairs"] == f.n_pairs)
        if not same:
            raise SystemExit(f"the first block does not reproduce Figure 3 for {rw['program']}")


def main():
    """Both runs (all cells, then without the study PMID 29942094), each with all T and NK cells and the four
    subtypes; writes tnk_subtypes.csv, its sidecar and the log."""
    ap = argparse.ArgumentParser(description=__doc__.split("RUN:")[0].strip(),
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nperm", type=int, default=1000, help="permutations (1,000, the count the paper prints)")
    a = ap.parse_args()
    LOGS.mkdir(parents=True, exist_ok=True)
    log = L.Log()
    t0 = time.time()
    gmt = four_program_gmt()
    F = pd.read_csv(L.out("fig2_rows.csv"))
    labels = dict(F[F.cell_set == "tnk"][["set_id", "label_as_printed"]].values)
    log(f"T AND NK SUBTYPES.  programs {', '.join(PROGRAMS)}; nperm {a.nperm}, cap {CAP}, seed {SEED}")
    rows = []
    for run, leave in (("all", None), (f"without_{LEAVE}", LEAVE)):
        S = L.load_set("tnk", CAP, log=lambda *x, **k: None, leave_out=leave)
        l2 = subtype_labels(S)
        outside = ~np.isin(l2, [s for s, _ in SUBTYPES])
        log(f"{run}: {S['N']:,} cells; outside the four subtypes {int(outside.sum()):,} (other T cells "
            f"{int((l2 == 'T_other').sum()):,}, no subtype label {int((l2 == '').sum()):,} in "
            f"{len(np.unique(S['block'][l2 == '']))} batch)")
        block = test_set(S, gmt, run, *ALL, labels, a.nperm, log)
        if leave is None and a.nperm == 1000:
            check_first_block(block)
            log("  the first block equals the T and NK cell rows of Figure 3")
        rows += block
        for st, label in SUBTYPES:
            sub = subset_set(S, l2 == st)
            rows += test_set(sub, gmt, run, st, label, labels, a.nperm, log)
    T = pd.DataFrame(rows)[[c for c, _ in COLUMNS]]
    o = L.out("tnk_subtypes.csv")
    T.to_csv(o, index=False)
    L.sidecar(o, "The four programs of T and NK cells in Figure 3, in all T and NK cells and within conventional CD4 "
              "T cells, regulatory T cells, CD8 T cells and NK cells, with all cells and without the cells of "
              "PMID 29942094 (12_tnk_subtypes.py; the sheet T and NK subtypes of Supplementary Data 1)", COLUMNS)
    log.write(LOGS / "12_tnk_subtypes.log")
    print(f"12_tnk_subtypes.py done in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
