#!/usr/bin/env python3
"""
Chapter 2, step 03: the TXNRD1 contrast per gene program and cell set, and the histology arm.

WHAT: s1, per cell set (fibroblasts, macrophages, T and NK cells), for every program of the cell set's panel
      (chapter2/data/figure_panels.json; each cell set scores only the programs of its panel) with at least 3
      tested genes (TXNRD1 removed from every list):
        - effect: mean paired difference in the program score (mean within-batch z of the count
          fraction of its tested genes), TXNRD1-detecting minus a non-detecting cell of the same batch
          and histology matched on genes detected (1:1, nearest depth, ties at random, SEED = 6);
        - z and p from 1,000 permutations of TXNRD1 detection within patient (the cells of each patient keep
          their patient, the detection labels are shuffled among them) with the matcher re-run in each
          (p_gauss: normal p of z, as printed in Figure 3; p_emp: empirical);
        - control A (cells split by residual library size at TXNRD1's prevalence) and the verdict rule;
        - batches in the same direction (batches with at least 20 valid pairs);
        - the histology arm: the same pairs split by histology, diff = effect in LUSC minus effect in
          LUAD, tested against the same permutations (z = diff / SD of the permuted diffs, normal p,
          empirical p), with Benjamini-Hochberg q over the Figure 3 programs of the cell set.
      Each cell set is capped at 2,000 cells per batch and histology (drawn with SEED = 6).
      s2: the results table of the three cell sets with the program definitions (drawn as Supplementary
          Figure 3), the Figure 3 rows, the histology summary and the chapter 2 counts.
ITEMS: Figure 3 (files fig2_programs, fig2_rows.csv; plotted by 05_plots.py), Supplementary Figures 3 and 4
      (program_results.csv; 05_plots.py), the numbers of the Figure 3 legend and the Results text
      (ch2_numbers.csv, histology_arm_summary.csv). The leave-out runs that back the paper are
      --leave-out-study Jimenez_Bruno_2026 (the ablation study GSE325414; Methods sentence on it) and, for T and
      NK cells, --leave-out-study Guo_Zhang_2018 (the regulatory T-cell sentence).
FROM: work_root/intermediate/chapter2/work_<cell set>/ (02_working_sets.py) and
      work_root/output/chapter2/gene_programs.gmt, gene_programs_definitions.csv,
      gene_programs_nested.csv (01_programs.py).
TO:   work_root/output/chapter2/, every table with a sibling .md:
        s1  program_tests_<cell>.csv            one row per program, every check
            histology_arm_<cell>.csv            one row per program, the histology arm
            cellset_<cell>.csv                  patients, matched pairs and batches of the set (read by s2)
            with --leave-out-study, the three tables with the suffix _without_<prefix> and
            leftout_<cell>_without_<prefix>.csv  cells, patients and batches left out, patients and pairs
                                                with and without them
        s2  program_results.csv, fig2_rows.csv, histology_arm_summary.csv, ch2_numbers.csv
      and the logs work_root/logs/chapter2/03_<cell>.log and 03_assemble.log.
WHY:  one pass computes the program tests and the histology arm from the same permuted pairings. The
      permuted pairings are drawn first, in order, from the one generator; the per-pairing statistics are
      then computed by --workers processes and accumulated in the original order, so the result does
      not depend on the number of workers. The matcher is match_fast (ch2lib), which is checked against
      the reference matcher match_random_ties at the start of every run; --slow-matcher uses the
      reference matcher throughout.
RUN:  python3 chapter2/03_program_tests.py [--only s1 s2] [--cells fib mac tnk] [--permutations 1000]
      [--workers N] [--slow-matcher] [--leave-out-study PREFIX [PREFIX ...]]
      --permutations (alias --nperm) defaults to 1,000, the count the manuscript prints; a smaller count
      (for example 250) only checks that the code runs: the effects, pairs, control A effects and batch
      votes do not depend on it, the z, p, bands and verdicts do.
      Seeds: SEED = 6 (cap, observed pairing, permutations).
"""
import argparse
import multiprocessing as mp
import os
import shutil
import sys
import tempfile
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ch2lib as L  # noqa: E402
from redox import paths  # noqa: E402

SEED = L.SEED
CAP = 2000                                    # cells per batch x histology in the tests
LOGS = paths.WORK_ROOT / "logs" / "chapter2"

VARIABLE = ("mean paired difference in the program score (mean within-batch z of the count fraction of its "
            "tested genes)")
DIRECTION = "positive = higher in the TXNRD1-detecting cell of a depth-matched pair"


def definitions():
    """The program definitions table written by 01_programs.py (gene_programs_definitions.csv)."""
    return pd.read_csv(L.out("gene_programs_definitions.csv"), dtype=str).fillna("")


def figure2_selection(T, cell):
    """Figure 3 rows of one cell set: primary or contrast programs of its panel whose verdict is
    "specific to TXNRD1", ordered by z, a nested program placed directly under its parent."""
    cfg = L.panels()
    d = T[T.role.isin(cfg["figure_roles"]) & T.verdict.str.startswith("specific")
          & T.set_id.map(L.stem).isin(cfg["panel"][cell])].sort_values("z", ascending=False)
    order = list(d.set_id)
    for child, spec in cfg["nested"].items():
        if child in order and spec["parent"] in order:
            order.remove(child)
            order.insert(order.index(spec["parent"]) + 1, child)
    return d.set_index("set_id").loc[order].reset_index()


# ====================================================================== s1
def panel_gmt(cell):
    """The lines of gene_programs.gmt whose program is on the panel of one cell set, written to a temporary
    file. A program's score depends only on its own genes, so the scores equal those of the full file."""
    stems = set(L.panels()["panel"][cell])
    lines = [ln for ln in open(L.out("gene_programs.gmt")) if L.stem(ln.split("\t")[0]) in stems]
    path = Path(tempfile.mkdtemp()) / f"panel_{cell}.gmt"
    path.write_text("".join(lines))
    return path


def s1_tests(cell, a):
    """Tests of one cell set: observed pairing, 1,000 within-patient permutations (SEED = 6) for the TXNRD1 split
    and for control A, the verdicts, the histology arm, and the counts of the cell set; writes the s1 tables."""
    spec = L.CELL_SETS[cell]
    log = L.Log()
    S = L.load_set(cell, CAP, leave_out=a.leave_out_study)
    sfx = L.leave_out_suffix(a.leave_out_study)
    # only the programs of this cell set's panel (chapter2/data/figure_panels.json) are scored and tested
    P = L.score_programs(S, panel_gmt(cell))
    PROG, names, R = P["PROG"], P["names"], P["R"]
    N, TX, depth, strat, pidx = S["N"], S["TX"], S["depth"], S["strat"], S["pidx"]
    block, pat, hist = S["block"], S["pat"], S["hist"]
    studies = S["studies"]
    strata = L.Strata(strat)
    log("=" * 110)
    log(f"GENE PROGRAMS, TXNRD1 CONTRAST.   CELL SET: {spec['label']}")
    log(f"NPERM={a.nperm}, cap {CAP} per batch x histology, "
        f"matcher {'verbatim' if a.slow_matcher else 'match_fast'}, workers {a.workers}")
    log(f"  {N:,} cells ({S['n_before']:,} before the cap), {len(np.unique(pat))} patients, {len(studies)} batches, "
        f"TXNRD1 detected in {TX.mean():.4f}")
    log(f"  {len(S['keepg']):,} genes pass arm presence {L.PRESENCE} and detection {L.DETFLOOR}")
    log(f"  {len(PROG)} programs in the file; {R} have at least {L.MINSET} tested genes and are scored")
    for n, p in PROG.items():
        if n not in names:
            log(f"    not scored, {len(p['tested'])} tested genes: {n}")
    log(f"  {P['G']:,} member genes z-scored within batch   [{log.elapsed()}]")

    if a.slow_matcher:
        def match(e, r):
            """Pairs of exposed cells e, by the reference matcher."""
            return L.match_random_ties(depth, e, strat, r)
    else:
        L.check_matcher(depth, TX, strat)
        log("  match_fast returns the pairs of the verbatim matcher on 3 successive draws (checked)")

        def match(e, r):
            """Pairs of exposed cells e, by match_fast."""
            return L.match_fast(depth, e, strata, r)

    both = sorted(b for b in np.unique(block) if {"LUSC", "LUAD"} <= set(hist[block == b]))
    scorer = L.Scorer.build(P, S)
    CTRLA = L.control_a(S)

    log(f"\n  matching, {a.nperm} within-patient permutation replicates for each of two arms ...")
    RES = {}
    HKEYS = [f"all_{h}" for h in ("ALL", "LUSC", "LUAD", "DIFF")]
    for nm, e in (("REAL", TX), ("CONTROL_A", CTRLA)):
        aa, bb = match(e, np.random.default_rng(SEED))
        obs, D, V, cnt = scorer.base(aa, bb)
        # the permuted pairings, drawn in order from one generator
        r = np.random.default_rng(SEED)
        pairs = []
        for i in range(a.nperm):
            pa, pb = match(L.perm_expo(e, pidx, r), r)
            if len(pa) >= L.MINPAIR:
                pairs.append((pa, pb))
        log(f"    {nm:10s} {len(pairs):,} permuted pairings drawn   [{log.elapsed()}]")
        nul = []
        if nm == "REAL":
            OBS_H = scorer.hist(aa, D, V)
            PERM_H = {k: [] for k in HKEYS}
            for k, res in enumerate(real_stats(scorer, pairs, a.workers, cell)):
                nul.append(res["e"])
                for key in HKEYS:
                    PERM_H[key].append(res["hist"][key])
                if (k + 1) % 100 == 0:
                    log(f"    {nm:10s} replicate {k + 1:,} / {len(pairs):,}   [{log.elapsed()}]")
            REAL_AB = (aa, bb, D, V)
        else:
            for pa, pb in pairs:
                nul.append(scorer.base(pa, pb)[0])
        nul = np.vstack(nul)
        mu, sd = np.nanmean(nul, 0), np.nanstd(nul, 0)
        with np.errstate(invalid="ignore", divide="ignore"):
            z = np.where(sd > 0, (obs - mu) / sd, np.nan)
        pemp = (np.sum(np.abs(nul - mu) >= np.abs(obs - mu), 0) + 1) / (len(nul) + 1)
        RES[nm] = dict(eff=obs, z=z, sd=sd, p=2 * stats.norm.sf(np.abs(z)), p_emp=pemp, n=cnt, npairs=len(aa),
                       nperm=len(nul))
        log(f"    {nm:10s} pairs {len(aa):,}, {len(nul):,} replicates used   [{log.elapsed()}]")

    # ---------------------------------------------------------------- the program table
    def verdict(e, z, ce, cz):
        """Verdict of one program from its effect e and z, and the effect ce and z cz of control A (|z| >= 1.96 counts)."""
        if not np.isfinite(z) or abs(z) < L.Z95:
            return "no effect, |z| below 1.96"
        same = np.isfinite(cz) and abs(cz) >= L.Z95 and np.sign(ce) == np.sign(e)
        if same and abs(ce) >= abs(e):
            return "withdrawn, RNA amount reproduces it at least as strongly"
        if same:
            return "not specific, RNA amount moves it the same way"
        return "specific to TXNRD1"

    Rr = RES["REAL"]
    aa, bb, Dv, Vv = REAL_AB
    ba = block[aa]
    rows = []
    for k, nm in enumerate(names):
        e, z = Rr["eff"][k], Rr["z"][k]
        sgn = np.sign(e)
        v = Vv[:, k]
        d = Dv[v, k]
        bb_ = ba[v]
        bm = np.array([d[bb_ == u].mean() for u in np.unique(bb_) if (bb_ == u).sum() >= L.MINBLKPAIR])
        p = PROG[nm]
        nsame = int((bm > 0).sum()) if sgn >= 0 else int((bm < 0).sum())
        rows.append(dict(
            cell_type=spec["label"], set_id=nm, n_listed=len(p["listed"]), n_tested=len(p["tested"]),
            txnrd1_removed=p["splitter_removed"], not_testable=";".join(p["not_testable"]),
            effect=e, band95_lo=e - L.Z95 * Rr["sd"][k], band95_hi=e + L.Z95 * Rr["sd"][k], z=z,
            p_gauss=Rr["p"][k], p_emp=Rr["p_emp"][k], n_pairs=int(Rr["n"][k]), n_replicates=Rr["nperm"],
            controlA_effect=RES["CONTROL_A"]["eff"][k], controlA_z=RES["CONTROL_A"]["z"][k],
            verdict=verdict(e, z, RES["CONTROL_A"]["eff"][k], RES["CONTROL_A"]["z"][k]),
            study_vote=f"{nsame} of {len(bm)}", batches_same_direction=nsame, batches_with_20_pairs=len(bm),
            study_vote_floor=float(2 / 2 ** len(bm)) if 2 <= len(bm) <= 18 else np.nan))
    T = pd.DataFrame(rows)
    npat = len(np.unique(pat))
    T.insert(2, "variable", VARIABLE)
    T.insert(3, "unit", "within-batch z")
    T.insert(4, "direction", DIRECTION)
    T.insert(5, "denominator", f"{Rr['npairs']:,} matched pairs, {N:,} cells, {npat} patients, {len(studies)} batches")
    write_program_tables(cell + sfx, T)

    # ---------------------------------------------------------------- histology arm
    # the permuted effects per histology come from the same pairings; z = (observed - permutation mean) / permutation SD
    Dd = definitions().set_index("set_id")
    T["role"] = T.set_id.map(Dd.role)
    f2 = set(figure2_selection(T, cell).set_id)
    PERM_H = {k: np.vstack(v) for k, v in PERM_H.items()}
    H = pd.DataFrame({"cell_set": spec["label"], "set_id": names,
                      "n_tested": [len(PROG[n]["tested"]) for n in names], "in_fig2": [n in f2 for n in names]})
    with np.errstate(invalid="ignore", divide="ignore"):
        for h in ("ALL", "LUSC", "LUAD"):
            k = f"all_{h}"
            Pm = PERM_H[k]
            mu, sd = np.nanmean(Pm, 0), np.nanstd(Pm, 0)
            H[f"all_effect_{h}"] = OBS_H[k]
            H[f"all_sd_{h}"] = sd
            H[f"all_z_{h}"] = np.where(sd > 0, (OBS_H[k] - mu) / sd, np.nan)
            H[f"all_n_pairs_{h}"] = OBS_H[f"all_n_{h}"]
            H[f"all_npairs_{h}_total"] = OBS_H[f"all_npairs_{h}"]
        Pm = PERM_H["all_DIFF"]
        d = OBS_H["all_DIFF"]
        sd = np.nanstd(Pm, 0)
        nn = np.isfinite(Pm).sum(0)
        H["all_diff"] = d
        H["all_perm_sd_diff"] = sd
        H["all_n_replicates"] = nn
        z = np.where(sd > 0, d / sd, np.nan)
        H["all_z"] = z
        H["all_p_gauss"] = 2 * stats.norm.sf(np.abs(z))
        H["all_p_emp"] = np.where(np.isfinite(d), ((np.abs(np.nan_to_num(Pm, nan=0)) >= np.abs(d))
                                                   & np.isfinite(Pm)).sum(0) + 1, np.nan) / (nn + 1)
        qf = np.full(R, np.nan)
        m = H.in_fig2.values
        qf[m] = L.bh(H.loc[m, "all_p_gauss"])
        H["all_q_fig2"] = qf
        H["all_n_batches"] = len(np.unique(block))
        H["all_n_patients_LUSC"] = len(np.unique(pat[hist == "LUSC"]))
        H["all_n_patients_LUAD"] = len(np.unique(pat[hist == "LUAD"]))
    o = L.out(f"histology_arm_{cell}{sfx}.csv")
    H.to_csv(o, index=False)
    histology_sidecar(o, spec["label"], a.nperm)
    log(f"\n  histology arm: {len(both)} of {len(studies)} batches carry both histologies; Figure 3 programs with "
        f"all-pairs q < 0.05: {int((H.all_q_fig2 < 0.05).sum())}")

    # ---------------------------------------------------------------- the cell set in numbers
    # the counts s2 turns into the Chapter 2 numbers (the log above has the rest)
    C = pd.DataFrame([dict(cell_set=cell, label=spec["label"], patients=npat, matched_pairs=Rr["npairs"],
                           batch_names=";".join(studies))])
    o = L.out(f"cellset_{cell}{sfx}.csv")
    C.to_csv(o, index=False)
    L.sidecar(o, f"The {spec['short']} working set of the tests (cap {CAP} per batch x histology, seed {SEED}): the "
              "counts ch2_numbers.csv is built from (s2).", [
        ("cell_set", "fib, mac or tnk"), ("label", "printed name"), ("patients", "patients after the cap"),
        ("matched_pairs", "matched pairs of the observed pairing"),
        ("batch_names", "the sequencing batches (block_id: a study, or one platform within a study), semicolon "
                        "separated; s2 counts the studies from them")])
    if sfx:
        write_left_out(cell, spec, sfx, a.leave_out_study, C)
    log(f"wrote program_tests_{cell}.csv and histology_arm_{cell}.csv   [{log.elapsed()}]")
    log.write(LOGS / f"03_{cell}{sfx}.log")


def write_left_out(cell, spec, sfx, leave_out, C):
    """leftout_<cell>_without_<prefix>.csv (leave-out runs only): what the left-out studies held in the
    working set before the cap, and the cell set with and without them (the ordinary cellset_<cell>.csv)."""
    mA, mB = L.load_meta(L.work_dir(cell))
    block = np.concatenate([mA["block"].astype(str), mB["block"].astype(str)])
    pat = np.concatenate([mA["patient"].astype(str), mB["patient"].astype(str)])
    gone = L.left_out(block, leave_out)
    full = L.out(f"cellset_{cell}.csv")
    F = pd.read_csv(full) if full.exists() else None
    row = dict(cell_set=cell, label=spec["label"], left_out=" and ".join(L.leave_out_prefixes(leave_out)),
               cells_left_out=int(gone.sum()), cells_before_cap=len(block), patients_left_out=len(np.unique(pat[gone])),
               batches_left_out=";".join(sorted(set(block[gone]))),
               patients_all=int(F.patients.iloc[0]) if F is not None else np.nan, patients_without=int(C.patients.iloc[0]),
               matched_pairs_all=int(F.matched_pairs.iloc[0]) if F is not None else np.nan,
               matched_pairs_without=int(C.matched_pairs.iloc[0]))
    o = L.out(f"leftout_{cell}{sfx}.csv")
    pd.DataFrame([row]).to_csv(o, index=False)
    L.sidecar(o, f"The {spec['short']} working set with and without the left-out studies (03 --leave-out-study).", [
        ("cell_set", "fib, mac or tnk"), ("label", "printed name"), ("left_out", "the left-out batch-name prefixes"),
        ("cells_left_out", "cells of the left-out studies in the working set, before the cap"),
        ("cells_before_cap", "all cells of the working set, before the cap"),
        ("patients_left_out", "patients of the left-out studies in the working set"),
        ("batches_left_out", "their batches, semicolon separated"),
        ("patients_all", "patients after the cap, all studies (cellset_<cell>.csv)"),
        ("patients_without", "patients after the cap, without the left-out studies"),
        ("matched_pairs_all", "matched pairs, all studies (cellset_<cell>.csv)"),
        ("matched_pairs_without", "matched pairs without the left-out studies (the cap and the pairing are redrawn "
                                  "from seed 6, so the difference is not exactly the pairs of the left-out studies)")])


def real_stats(scorer, pairs, workers, cell):
    """The REAL-arm statistics of every permuted pairing, in the order of the pairings."""
    if workers <= 1:
        for pa, pb in pairs:
            yield scorer.real(pa, pb)
        return
    folder = paths.intermediate("chapter2") / f"_scorer_{cell}"
    scorer.save(folder)
    ctx = mp.get_context("spawn")
    try:
        with ctx.Pool(workers, initializer=L.worker_init, initargs=(str(folder),)) as pool:
            for res in pool.imap(L.worker_real, pairs, chunksize=4):
                yield res
    finally:
        shutil.rmtree(folder, ignore_errors=True)


TEST_COLUMNS = [
    ("cell_type", "cell set"), ("set_id", "program id"), ("variable", "what the effect is"), ("unit", "unit"),
    ("direction", "sign convention"), ("denominator", "matched pairs, cells, patients and batches"),
    ("n_listed", "genes in the list"), ("n_tested", "genes tested (TXNRD1 removed)"),
    ("txnrd1_removed", "True when TXNRD1 was in the list"), ("not_testable", "listed genes not tested"),
    ("effect", "mean paired difference, TXNRD1-detecting minus not"),
    ("band95_lo", "effect minus 1.96 permutation SD"), ("band95_hi", "effect plus 1.96 permutation SD"),
    ("z", "(effect minus the permutation mean) / permutation SD"), ("p_gauss", "two-sided normal p of z (printed in Figure 3)"),
    ("p_emp", "empirical permutation p"), ("n_pairs", "valid pairs"), ("n_replicates", "permutations used"),
    ("controlA_effect", "effect of the library-size split (control A)"), ("controlA_z", "its z"),
    ("verdict", "no effect / withdrawn / not specific / specific to TXNRD1 (rule in README)"),
    ("study_vote", "k of n batches (with at least 20 valid pairs) with the sign of the effect"),
    ("batches_same_direction", "k of study_vote"), ("batches_with_20_pairs", "n of study_vote"),
    ("study_vote_floor", "exact floor 2 / 2^n")]


def write_program_tables(cell, T):
    """program_tests_<cell>.csv with its sidecar; the cell name may carry the _without_ suffix of a leave-out run."""
    label = L.CELL_SETS[cell.split("_without_")[0]]["label"]
    o = L.out(f"program_tests_{cell}.csv")
    T.to_csv(o, index=False)
    L.sidecar(o, f"The TXNRD1 contrast per program in {label}, one row per scored program of its panel.",
              TEST_COLUMNS)


def histology_sidecar(o, label, nperm):
    """The .md sidecar of histology_arm_<cell>.csv."""
    cols = [("cell_set", "cell set"), ("set_id", "program id"), ("n_tested", "genes tested"),
            ("in_fig2", "True for the programs drawn in Figure 3 for this cell set")]
    for h in ("ALL", "LUSC", "LUAD"):
        cols += [(f"all_effect_{h}", f"effect, all pairs, {h if h != 'ALL' else 'both histologies'}"),
                 (f"all_sd_{h}", "its permutation SD"), (f"all_z_{h}", "its permutation z"),
                 (f"all_n_pairs_{h}", "valid pairs"), (f"all_npairs_{h}_total", "pairs")]
    cols += [("all_diff", "effect in LUSC minus effect in LUAD, all pairs"),
             ("all_perm_sd_diff", "SD of the permuted diffs"), ("all_n_replicates", "permutations with a diff"),
             ("all_z", "diff / SD of the permuted diffs (primary statistic)"),
             ("all_p_gauss", "two-sided normal p of that z"), ("all_p_emp", "empirical p"),
             ("all_q_fig2", "BH q of all_p_gauss over the Figure 3 programs of the cell set"),
             ("all_n_batches", "batches"), ("all_n_patients_LUSC", "LUSC patients"),
             ("all_n_patients_LUAD", "LUAD patients")]
    L.sidecar(o, f"{label}: does the TXNRD1 effect on a program differ between LUSC and LUAD? "
              f"The pairs of the Figure 3 test split by histology, against the same {nperm:,} within-patient permutations. "
              "The prefix all_ means all pairs.", cols)


# ====================================================================== s2
RES_COLS = ["cell_type", "set_id", "program", "role", "group", "variable", "unit", "direction", "denominator", "n_listed",
            "n_tested", "effect", "band95_lo", "band95_hi", "z", "p_gauss", "p_emp", "n_pairs", "n_replicates",
            "controlA_effect", "controlA_z", "verdict", "study_vote", "batches_same_direction", "batches_with_20_pairs",
            "study_vote_floor"]


def s2_assemble(a):
    """Assemble the three cell sets: program_results.csv, fig2_rows.csv (the Figure 3 rows), the histology
    summary and the counts of the legend and text (ch2_numbers.csv)."""
    log = L.Log()
    D = definitions()
    cfg = L.panels()
    RS = pd.concat([pd.read_csv(L.out(f"program_tests_{cell}.csv"), float_precision="round_trip")
                    for cell in L.CELL_SETS], ignore_index=True)
    RS = RS.merge(D[["set_id", "program", "role", "group"]], on="set_id", how="left")
    RS = RS[RES_COLS]
    o = L.out("program_results.csv")
    RS.to_csv(o, index=False)
    meaning = dict(TEST_COLUMNS, program="label (catalog, US spelling)", role="role of the list", group="group")
    L.sidecar(o, "The TXNRD1 contrast for every scored program in the three cell sets, with the program's role and "
              "group (the published results table of all cell types without the pooled rows; drawn as Supplementary "
              "Figure 3). The columns of program_tests_<cell>.csv and the program's label, role and group:",
              [(c, meaning[c]) for c in RES_COLS])

    # the rows of Figure 3
    rows = []
    for cell, spec in L.CELL_SETS.items():
        d = figure2_selection(RS[RS.cell_type == spec["label"]], cell)
        d = d.assign(cell_set=cell, cell_block=spec["label"],
                     label_as_printed=[cfg["printed_labels"].get(s, L.us(p)) for s, p in zip(d.set_id, d.program)],
                     nested_under=[cfg["nested"].get(s, {}).get("parent", "") for s in d.set_id],
                     row=np.arange(1, len(d) + 1))
        rows.append(d)
    F = pd.concat(rows, ignore_index=True)
    F = F[["cell_set", "cell_block", "row", "set_id", "label_as_printed", "nested_under", "effect", "band95_lo",
           "band95_hi", "z", "p_gauss", "study_vote", "batches_same_direction", "batches_with_20_pairs",
           "n_pairs", "denominator"]]
    o = L.out("fig2_rows.csv")
    F.to_csv(o, index=False)
    L.sidecar(o, "The rows of Figure 3: primary or contrast programs of each cell set's panel "
              "(chapter2/data/figure_panels.json) whose verdict is 'specific to TXNRD1', ordered by z, the nested "
              "row under its parent.", [
                  ("cell_set", "fib, mac or tnk"), ("cell_block", "block title"), ("row", "row within the block"),
                  ("set_id", "program id"), ("label_as_printed", "row label"), ("nested_under", "parent row of a nested row"),
                  ("effect", "bar: mean paired difference"), ("band95_lo", "line, left end"), ("band95_hi", "line, right end"),
                  ("z", "permutation z"), ("p_gauss", "p as printed"), ("study_vote", "batches in the same direction, as printed"),
                  ("batches_same_direction", "k"), ("batches_with_20_pairs", "n"), ("n_pairs", "valid pairs"),
                  ("denominator", "pairs, cells, patients, batches of the cell set")])

    # histology summary (all pairs, the Figure 3 programs of each cell set)
    srows = []
    for cell, spec in L.CELL_SETS.items():
        T = pd.read_csv(L.out(f"histology_arm_{cell}.csv"), float_precision="round_trip")
        sub = T[T.in_fig2]
        p, q = sub["all_p_gauss"], sub["all_q_fig2"]
        ok = p.notna()
        srows.append(dict(
            cell_set=spec["label"], analysis="all pairs", family="Figure 3 programs", n_programs=int(ok.sum()),
            pairs_LUSC=int(T["all_npairs_LUSC_total"].iloc[0]), pairs_LUAD=int(T["all_npairs_LUAD_total"].iloc[0]),
            batches=int(T["all_n_batches"].iloc[0]),
            n_p_below_005=int((p < 0.05).sum()), expected_by_chance=round(0.05 * ok.sum(), 2),
            n_p_emp_below_005=int((sub["all_p_emp"] < 0.05).sum()),
            n_q_below_005=int((q < 0.05).sum()), n_q_below_020=int((q < 0.20).sum()),
            min_p=float(p.min()), min_q=float(q.min()),
            smallest_p_program=str(sub.loc[p.idxmin(), "set_id"]) if ok.any() else "",
            programs_q_below_005=";".join(sub.loc[q < 0.05, "set_id"])))
    Sm = pd.DataFrame(srows)
    o = L.out("histology_arm_summary.csv")
    Sm.to_csv(o, index=False)
    L.sidecar(o, "Histology arm, per cell set: the Figure 3 programs, all pairs.", [
        ("cell_set", "cell set"), ("analysis", "all pairs"),
        ("family", "Figure 3 programs of the cell set (the Benjamini-Hochberg family)"),
        ("n_programs", "programs with a p"),
        ("pairs_LUSC", "LUSC pairs"), ("pairs_LUAD", "LUAD pairs"), ("batches", "batches"),
        ("n_p_below_005", "programs with p_gauss < 0.05"), ("expected_by_chance", "0.05 x programs"),
        ("n_p_emp_below_005", "programs with empirical p < 0.05"),
        ("n_q_below_005", "programs with BH q < 0.05 in the family"), ("n_q_below_020", "q < 0.20"),
        ("min_p", "smallest p_gauss"), ("min_q", "smallest q"), ("smallest_p_program", "program with the smallest p"),
        ("programs_q_below_005", "programs with q < 0.05")])

    # chapter 2 numbers
    C = pd.concat([pd.read_csv(L.out(f"cellset_{c}.csv")) for c in L.CELL_SETS], ignore_index=True)
    batches = sorted({b for s in C.batch_names for b in s.split(";")})
    studies = sorted({L.study_of_batch(b) for b in batches})
    nest = pd.read_csv(L.out("gene_programs_nested.csv"))
    items = [("studies_all", len(studies), "studies behind the three cell sets (batch names without the platform part)"),
             ("patients_min", int(C.patients.min()), "patients per cell set, smallest"),
             ("patients_max", int(C.patients.max()), "patients per cell set, largest")]
    for r in C.itertuples():
        items += [(f"pairs_{r.cell_set}", r.matched_pairs, f"matched pairs, {r.label}")]
    for r in nest.itertuples():
        items += [(f"genes_{r.child}", r.n_child, f"genes of {r.child}"),
                  (f"genes_{r.child}_in_{r.parent}", r.n_child_in_parent, f"of them in {r.parent}")]
    items.append(("hist_fig2_q005_all_cells", int(Sm.n_q_below_005.sum()),
                  "Figure 3 programs differing between LUSC and LUAD at BH q < 0.05 (all pairs, three cell sets)"))
    N2 = pd.DataFrame(items, columns=["item", "value", "meaning"])
    o = L.out("ch2_numbers.csv")
    N2.to_csv(o, index=False)
    L.sidecar(o, "Counts quoted in the Chapter 2 text and legends.", [("item", "key"), ("value", "value"),
                                                                        ("meaning", "what it counts")])
    log(f"s2: program_results.csv {len(RS)} rows, fig2_rows.csv {len(F)} rows, {len(studies)} studies in "
        f"{len(batches)} batches, patients {int(C.patients.min())} to {int(C.patients.max())}")
    log.write(LOGS / "03_assemble.log")


def main():
    """Run s1 for each cell set and/or s2; a leave-out run is s1 only."""
    ap = argparse.ArgumentParser(description=__doc__.split("RUN:")[0].strip(),
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", nargs="*", default=["s1", "s2"], choices=["s1", "s2"])
    ap.add_argument("--cells", nargs="*", default=list(L.CELL_SETS), choices=list(L.CELL_SETS))
    ap.add_argument("--permutations", "--nperm", dest="nperm", type=int, default=1000,
                    help="within-patient permutations per arm (default 1,000, the published count)")
    ap.add_argument("--workers", type=int, default=max(1, min(6, (os.cpu_count() or 2) - 1)))
    ap.add_argument("--slow-matcher", action="store_true", help="use the reference matcher throughout (match_fast is checked against it at the start of every run)")
    ap.add_argument("--leave-out-study", nargs="+", default=None, metavar="PREFIX",
                    help="batch-name prefix of a study to leave out, or several (writes *_without_<prefix> tables, "
                         "*_without_<prefix 1>_and_<prefix 2> for several; s1 only). The manuscript uses Guo_Zhang_2018 "
                         "for T and NK cells: the regulatory T-cell signature difference disappears without it "
                         "(Results); Jimenez_Bruno_2026 (GSE325414) for the sensitivity analysis of the Methods")
    a = ap.parse_args()
    if a.leave_out_study and "s2" in a.only:
        raise SystemExit("--leave-out-study runs s1 only: add --only s1")
    LOGS.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    if "s1" in a.only:
        for cell in a.cells:
            s1_tests(cell, a)
    if "s2" in a.only:
        s2_assemble(a)
    print(f"03_program_tests.py done in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
