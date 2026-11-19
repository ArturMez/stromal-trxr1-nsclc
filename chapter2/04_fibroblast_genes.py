#!/usr/bin/env python3
"""
Chapter 2, step 04: the per-gene fibroblast scan and the per-batch gene effects of Supplementary Figure 4.

WHAT: s1, the transcriptome-wide TXNRD1 scan in fibroblasts: for every gene of the tested universe, the
          mean paired difference in the within-batch z of its count fraction, TXNRD1-detecting minus a
          depth-matched non-detecting cell of the same batch and histology, with 100 within-patient
          permutations of TXNRD1 detection (matcher re-run in each, SEED = 6), control A (the library-size
          split) and Benjamini-Hochberg q over the genes.
      s2, the genes of Supplementary Figure 4: a gene is listed when its q (normal p) is 0.20 or below and
          the library-size split does not move it the same way at |z| 1.96 or more; the lower group is
          every tested member of the complement cascade list. Then the per-batch effect of each drawn
          gene over the observed pairs, and the rows, dots and counts of the panel.
ITEMS: the gene tables behind Supplementary Figure 4 (the figure itself is drawn by 05_plots.py), and the
      numbers of its legend and footnote; the values of the sheet Aldehyde enzyme genes of Supplementary Data 1
      (fibroblast_genes_scan.csv, read by 10_supplementary_data1.py). The leave-out run that backs the paper is
      --leave-out-study Jimenez_Bruno_2026 (the ablation study GSE325414; Methods sentence on it).
FROM: work_root/intermediate/chapter2/work_fib/ (02_working_sets.py); work_root/output/chapter2/
      gene_programs_definitions.csv (01) and program_tests_fib.csv (03, the not-testable members).
TO:   work_root/output/chapter2/, every table with a sibling .md:
        s1  fibroblast_genes_scan.csv       every tested gene: effect, permutation mean and SD, z, p, q,
                                            control A, verdict
        s2  suppfig4_genes.csv            every gene considered for the panel and why it is in or out
            suppfig4_rows.csv             the panel rows, in order, with the pooled effect and band
            suppfig4_dots.csv             gene x batch effects (the dots), with the batch label and pairs
            suppfig4_numbers.csv          the counts of the legend and the footnote
      and the log work_root/logs/chapter2/04_fibroblast_genes.log.
WHY:  Supplementary Figure 4 shows which single genes, beyond the programs of Figure 3, move with TXNRD1
      in fibroblasts, and whether the effect is seen batch by batch. The gene choice uses this script's own
      scan alone.
RUN:  python3 chapter2/04_fibroblast_genes.py [--only s1 s2] [--nperm 100] [--leave-out-study PREFIX [PREFIX ...]]
      --leave-out-study removes the cells of one study or several (batch-name prefixes) before the
      cap, as 03_program_tests.py does, reads program_tests_fib_without_<prefix>.csv of 03 and writes every table
      and the log with the suffix _without_<prefix> beside the ordinary ones, which it leaves untouched.
      Seeds: SEED = 6 (cap, pairings, permutations); 100 permutations; cap 2,000 cells per batch and histology.
      About two minutes.
"""
import argparse
import os
import sys
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402
from scipy.sparse import vstack as spvstack  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ch2lib as L  # noqa: E402
from redox import paths  # noqa: E402

SEED = L.SEED
CAP = 2000
Q_RULE = 0.20
CELL = "fib"
LOGS = paths.WORK_ROOT / "logs" / "chapter2"
LOG = L.Log()
UNIT = "within-batch z units"


def s1_scan(a):
    """The per-gene scan: effect, permutation z, p, q and control A for every tested gene (fibroblast_genes_scan.csv)."""
    sfx = L.leave_out_suffix(a.leave_out_study)
    S = L.load_set(CELL, CAP, leave_out=a.leave_out_study)
    X, N, genes = S["X"], S["N"], S["genes"]
    pat, depth, strat = S["pat"], S["depth"], S["strat"]
    TX, PRES, acode, blks, blk_of, pidx = S["TX"], S["PRES"], S["acode"], S["blks"], S["blk_of"], S["pidx"]
    tot = S["tot"]
    blk_arms = [np.unique(acode[blk_of == b]) for b in range(len(blks))]
    if not all(len(u) == 1 for u in blk_arms):
        raise SystemExit("a batch spans two arms; the scan assumes one arm per batch")
    blk_arm = S["blk_arm"]
    keepg = S["keepg"]
    det = S["det"]
    strata = L.Strata(strat)
    LOG(f"PER-GENE TXNRD1 SCAN, fibroblasts. NPERM={a.nperm}, cap {CAP} per batch x histology")
    LOG(f"  cells {S['n_before']:,} before the cap, {N:,} after; patients {len(np.unique(pat))}, batches {len(blks)}")
    LOG(f"  TXNRD1 detected in {TX.mean():.4f} of cells ({int(TX.sum()):,} of {N:,})")
    LOG(f"  {len(keepg):,} genes pass arm presence {L.PRESENCE} and detection {L.DETFLOOR}")
    L.check_matcher(depth, TX, strat)

    CTRLA = L.control_a(S)

    def pairsets(expo, n):
        """The observed pairing and n permuted pairings (detection permuted within patient) for an exposure vector."""
        r = np.random.default_rng(SEED)
        obs = L.match_fast(depth, expo, strata, np.random.default_rng(SEED))
        P = []
        for _ in range(n):
            pa, pb = L.match_fast(depth, L.perm_expo(expo, pidx, r), strata, r)
            if len(pa) >= L.MINPAIR:
                P.append((pa, pb))
        return obs, P

    ARMS = {}
    for nm, e in (("REAL", TX), ("CONTROL_A", CTRLA)):
        ARMS[nm] = pairsets(e, a.nperm)
        LOG(f"    {nm:10s} observed pairs {len(ARMS[nm][0][0]):,}, replicates kept {len(ARMS[nm][1])}   [{LOG.elapsed()}]")

    res = {k: [] for k in ARMS}
    CH = 1000
    for c0 in range(0, len(keepg), CH):
        idx = keepg[c0:c0 + CH]
        G = len(idx)
        A = X[:, idx].toarray().astype(np.float32)
        A /= tot[:, None]
        Z = np.full((N, G), np.nan, np.float32)
        for b in range(len(blks)):
            mb = blk_of == b
            if mb.sum() < L.MINBLK:
                continue
            has = PRES[blk_arm[b]][idx]
            sub = A[mb]
            sd = sub.std(0)
            mu = sub.mean(0)
            ok = has & (sd > 0)
            Zb = np.full(sub.shape, np.nan, np.float32)
            Zb[:, ok] = (sub[:, ok] - mu[ok]) / sd[ok]
            Z[mb] = Zb
        del A
        fin = np.isfinite(Z)

        def mp(pa, pb):
            """Per-gene mean paired difference over the pairs (pa, pb), and the number of valid pairs."""
            V = fin[pa] & fin[pb]
            D = np.where(V, np.nan_to_num(Z[pa] - Z[pb]), 0.0)
            c = V.sum(0)
            return np.where(c >= L.MINPAIR, D.sum(0) / np.maximum(c, 1), np.nan), c

        for nm, ((aa, bb), P) in ARMS.items():
            obs, cnt = mp(aa, bb)
            nul = np.vstack([mp(pa, pb)[0] for pa, pb in P])
            mu, sd = np.nanmean(nul, 0), np.nanstd(nul, 0)
            with np.errstate(invalid="ignore", divide="ignore"):
                z = np.where(sd > 0, (obs - mu) / sd, np.nan)
            pgz = 2 * stats.norm.sf(np.abs(z))
            pe = (np.sum(np.abs(nul - mu) >= np.abs(obs - mu), 0) + 1) / (len(nul) + 1)
            for j in range(G):
                if np.isfinite(z[j]):
                    res[nm].append((genes[idx[j]], obs[j], mu[j], sd[j], z[j], pgz[j], pe[j], int(cnt[j]), det[idx[j]]))

        LOG(f"    genes {min(c0 + CH, len(keepg)):,} / {len(keepg):,}   [{LOG.elapsed()}]")

    COLS = ["gene", "eff", "perm_mean", "perm_sd", "z", "p_gauss", "p_emp", "n_pairs", "detection"]
    tab = {}
    for nm in ARMS:
        d = pd.DataFrame(res[nm], columns=COLS)
        d["bh_gauss"] = L.bh_all(d.p_gauss)
        tab[nm] = d
    M = tab["REAL"].merge(tab["CONTROL_A"][["gene", "eff", "z", "p_gauss"]]
                          .rename(columns={"eff": "ctrlA_eff", "z": "ctrlA_z", "p_gauss": "ctrlA_p_gauss"}),
                          on="gene", how="left")
    M["ctrlA_same_sign"] = np.sign(M.eff) == np.sign(M.ctrlA_eff)
    M["ctrlA_counts"] = M.ctrlA_same_sign & (M.ctrlA_z.abs() >= L.Z95)
    M["ctrlA_stronger"] = M.ctrlA_counts & (M.ctrlA_eff.abs() >= M.eff.abs())
    M["verdict"] = np.where(M.bh_gauss > Q_RULE, "not called",
                            np.where(M.ctrlA_stronger, "withdrawn, RNA amount reproduces it at least as strongly",
                                     np.where(M.ctrlA_counts, "not specific, RNA amount moves it the same way",
                                              "specific to TXNRD1")))
    M.insert(1, "compartment", "Fibroblasts")
    M.insert(2, "variable", "mean paired difference in the within-batch z of the gene's count fraction")
    M.insert(3, "unit", UNIT)
    M.insert(4, "direction", "positive = higher in the TXNRD1-detecting cell of the depth-matched pair")
    M = M.sort_values("p_gauss")
    o = L.out(f"fibroblast_genes_scan{sfx}.csv")
    M.to_csv(o, index=False)
    L.sidecar(o, f"Per-gene TXNRD1 scan in fibroblasts ({N:,} cells after a cap of {CAP} per batch x histology, "
              f"{len(ARMS['REAL'][0][0]):,} observed pairs, {a.nperm} within-patient permutations).", [
                  ("gene", "gene"), ("compartment", "cell set"), ("variable", "what eff is"), ("unit", "unit"),
                  ("direction", "sign convention"), ("eff", "mean paired difference"), ("perm_mean", "permutation mean"),
                  ("perm_sd", "permutation SD (the band of Supplementary Figure 4 is eff +- 1.96 perm_sd)"),
                  ("z", "(eff - perm_mean) / perm_sd"), ("p_gauss", "normal p of z"), ("p_emp", "empirical p"),
                  ("n_pairs", "valid pairs"), ("detection", "detection rate"), ("bh_gauss", "BH q of p_gauss over the genes"),
                  ("ctrlA_eff", "effect of the library-size split"), ("ctrlA_z", "its z"), ("ctrlA_p_gauss", "its p"),
                  ("ctrlA_same_sign", "control A has the sign of eff"),
                  ("ctrlA_counts", "and |ctrlA_z| >= 1.96 (the rule that removes a gene)"),
                  ("ctrlA_stronger", "and at least as large"), ("verdict", "not called (q > 0.20), withdrawn, not specific, "
                                                                           "or specific to TXNRD1")])
    LOG(f"  REAL arm genes at BH q <= {Q_RULE}: {int((tab['REAL'].bh_gauss <= Q_RULE).sum())}; "
        f"CONTROL_A {int((tab['CONTROL_A'].bh_gauss <= Q_RULE).sum())}")
    LOG(f"wrote fibroblast_genes_scan{sfx}.csv   [{LOG.elapsed()}]")


def perbatch_effects(gene_list, leave_out=None):
    """Per-batch mean paired difference of each gene over the observed pairs (cap 2,000 per batch x
    histology, the pairing of SEED = 6); no permutations. Returns the table, the number of pairs and of cells."""
    S = L.load_set(CELL, CAP, leave_out=leave_out)
    W, sel, nA = L.work_dir(CELL), S["sel"], S["nA"]
    X = spvstack([L.load_rows(W, "LUSC", sel[sel < nA]), L.load_rows(W, "LUAD", sel[sel >= nA] - nA)]).tocsr()
    block, arm, depth, strat, GI = S["block"], S["arm"], S["depth"], S["strat"], S["GI"]
    present = pd.read_pickle(L.work_dir(CELL) / "present.pkl")
    tot = np.asarray(X.sum(1)).ravel().astype(np.float64)
    TX = S["TX"]
    a, b = L.match_fast(depth, TX, L.Strata(strat), np.random.default_rng(SEED))
    use = [g for g in gene_list if g in GI]
    A = X[:, [GI[g] for g in use]].toarray().astype(np.float64) / np.maximum(tot, 1)[:, None]
    Z = np.full(A.shape, np.nan)
    for bl in np.unique(block):
        m = block == bl
        if m.sum() < L.MINBLK:
            continue
        ar = np.unique(arm[m])[0]
        for k, g in enumerate(use):
            if ar in present.columns and not bool(present[ar].values[GI[g]]):
                continue
            s = A[m, k]
            sd = s.std()
            if sd > 0:
                Z[m, k] = (s - s.mean()) / sd
    rows = []
    for k, g in enumerate(use):
        d = Z[a, k] - Z[b, k]
        ok = np.isfinite(d)
        for bl in np.unique(block[a]):
            inb = block[a] == bl
            v = ok & inb
            rows.append(dict(gene=g, block=bl, n_pairs=int(inb.sum()), n_valid=int(v.sum()),
                             effect=float(d[v].mean()) if v.any() else np.nan))
        rows.append(dict(gene=g, block="(all batches)", n_pairs=len(a), n_valid=int(ok.sum()), effect=float(d[ok].mean())))
    return pd.DataFrame(rows), len(a), S["N"]


def s2_panel(a):
    """The genes of Supplementary Figure 4: selection rule, per-batch dots, rows and legend counts."""
    sfx = L.leave_out_suffix(a.leave_out_study)
    PG = pd.read_csv(L.out(f"fibroblast_genes_scan{sfx}.csv"), float_precision="round_trip").set_index("gene")
    D = pd.read_csv(L.out("gene_programs_definitions.csv"), dtype=str).fillna("").set_index("set_id")
    T = pd.read_csv(L.out(f"program_tests_fib{sfx}.csv"), dtype={"not_testable": str}).fillna({"not_testable": ""}).set_index("set_id")
    cset = L.panels()["complement_set"]

    def rna(g):
        """True when control A (library-size split) moves gene g the same way at |z| >= 1.96, which removes it."""
        return bool(np.sign(PG.loc[g, "ctrlA_eff"]) == np.sign(PG.loc[g, "eff"]) and abs(PG.loc[g, "ctrlA_z"]) >= L.Z95)

    cand = sorted(set(PG.index[PG.bh_gauss <= Q_RULE]) - {L.SPLITTER})
    called = sorted([g for g in cand if not rna(g)], key=lambda g: -PG.loc[g, "eff"])
    removed = sorted(g for g in cand if rna(g))
    nt = set(str(T.loc[cset, "not_testable"]).split(";"))
    members = sorted([g for g in D.loc[cset, "genes_used"].split() if g in PG.index and g not in nt and g != L.SPLITTER],
                     key=lambda g: -PG.loc[g, "eff"])
    G = pd.DataFrame([dict(gene=g, q=PG.loc[g, "bh_gauss"], eff=PG.loc[g, "eff"], ctrlA_eff=PG.loc[g, "ctrlA_eff"],
                           ctrlA_z=PG.loc[g, "ctrlA_z"], removed_by_rna_amount=rna(g), drawn=g in called)
                      for g in cand])
    o = L.out(f"suppfig4_genes{sfx}.csv")
    G.to_csv(o, index=False)
    L.sidecar(o, f"Every gene with BH q <= {Q_RULE} in the fibroblast scan (TXNRD1 excluded), and whether the "
              "library-size rule removes it. The rule uses the genome-wide scan alone.", [
                  ("gene", "gene"), ("q", "BH q of the normal p"), ("eff", "pooled effect"),
                  ("ctrlA_eff", "library-size split effect"), ("ctrlA_z", "its z"),
                  ("removed_by_rna_amount", "True when control A has the sign of eff with |z| >= 1.96"),
                  ("drawn", "True when the gene is a row of the upper group")])
    need = list(dict.fromkeys(called + members))
    MEM, npairs_all, ncells = perbatch_effects(need, a.leave_out_study)
    MEMb = MEM[MEM.block != "(all batches)"]
    npairs = MEMb.drop_duplicates("block").set_index("block").n_pairs
    cols = [b for b in npairs.sort_values(ascending=False).index if npairs[b] >= L.MINBLKPAIR]
    small = [b for b in npairs.index if npairs[b] < L.MINBLKPAIR]
    LAB = L.batch_labels(list(npairs.index))
    dots = MEMb.assign(label=MEMb.block.map(LAB), drawn=MEMb.block.isin(cols),
                       column=MEMb.block.map({b: i + 1 for i, b in enumerate(cols)}))
    o = L.out(f"suppfig4_dots{sfx}.csv")
    dots.to_csv(o, index=False)
    L.sidecar(o, "Per-batch effect of each drawn gene over the observed pairs (the dots).", [
        ("gene", "gene"), ("block", "batch"), ("n_pairs", "observed pairs in the batch (dot size)"),
        ("n_valid", "pairs with a value for the gene"), ("effect", "mean paired difference in the batch (dot color)"),
        ("label", "column label"), ("drawn", "True for batches with at least 20 pairs"),
        ("column", "column number, batches ordered by pairs")])
    rows = []
    for section, genes_ in (("Genes passing the transcriptome-wide test", called),
                            ("Complement cascade genes (Reactome), every tested member", members)):
        rows.append(dict(kind="header", section=section, gene="", eff=np.nan, lo=np.nan, hi=np.nan, q=np.nan,
                         batches_same_way=np.nan, batches_drawn=np.nan))
        for g in genes_:
            e, sd = PG.loc[g, "eff"], PG.loc[g, "perm_sd"]
            v = dots[(dots.gene == g) & dots.drawn].set_index("block").effect.reindex(cols)
            fin = v[np.isfinite(v)]
            rows.append(dict(kind="gene", section=section, gene=g, eff=e, lo=e - L.Z95 * sd, hi=e + L.Z95 * sd,
                             q=PG.loc[g, "bh_gauss"], batches_same_way=int((np.sign(fin) == np.sign(e)).sum()),
                             batches_drawn=len(fin)))
    Rw = pd.DataFrame(rows)
    Rw.insert(0, "row", np.arange(len(Rw)))
    o = L.out(f"suppfig4_rows{sfx}.csv")
    Rw.to_csv(o, index=False)
    L.sidecar(o, "Rows of Supplementary Figure 4 in order.", [
        ("row", "row from the top"), ("kind", "header or gene"), ("section", "group"), ("gene", "gene"),
        ("eff", "pooled effect (bar)"), ("lo", "eff - 1.96 permutation SD"), ("hi", "eff + 1.96 permutation SD"),
        ("q", "BH q in the scan"), ("batches_same_way", "drawn batches with the sign of eff (last column, k)"),
        ("batches_drawn", "drawn batches with a value (last column, n)")])
    nums = [("pairs_all", npairs_all, "matched pairs of the observed pairing"),
            ("cells", ncells, "fibroblasts after the cap"),
            ("batches_with_pairs", len(npairs), "batches contributing pairs"),
            ("batches_drawn", len(cols), "batches with at least 20 pairs (drawn columns)"),
            ("pairs_in_drawn_batches", int(npairs[cols].sum()), "pairs in the drawn batches"),
            ("batches_not_drawn", len(small), "batches with fewer than 20 pairs"),
            ("genes_removed_rna", len(removed), "genes removed by the library-size rule"),
            ("genes_scanned", len(PG), "genes in the scan")]
    Nm = pd.DataFrame(nums, columns=["item", "value", "meaning"])
    Nm["genes"] = ""
    Nm.loc[Nm.item == "genes_removed_rna", "genes"] = " ".join(removed)
    Nm.loc[Nm.item == "batches_not_drawn", "genes"] = "; ".join(f"{LAB[b]} ({npairs[b]})" for b in small)
    o = L.out(f"suppfig4_numbers{sfx}.csv")
    Nm.to_csv(o, index=False)
    L.sidecar(o, "Counts of the Supplementary Figure 4 legend and of its in-panel footnote.", [
        ("item", "key"), ("value", "value"), ("meaning", "what it counts"),
        ("genes", "the genes or batches behind the count")])
    LOG(f"s2: upper group {called}; removed by the library-size rule {removed}; {len(members)} complement members; "
        f"{len(cols)} batches drawn ({int(npairs[cols].sum()):,} of {npairs_all:,} pairs)")


def main():
    """Run s1 and/or s2 and write the log."""
    ap = argparse.ArgumentParser(description=__doc__.split("RUN:")[0].strip(),
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", nargs="*", default=["s1", "s2"], choices=["s1", "s2"])
    ap.add_argument("--nperm", type=int, default=100)
    ap.add_argument("--leave-out-study", nargs="+", default=None, metavar="PREFIX",
                    help="batch-name prefix of a study to leave out, or several (writes *_without_<prefix> tables)")
    a = ap.parse_args()
    LOGS.mkdir(parents=True, exist_ok=True)
    if "s1" in a.only:
        s1_scan(a)
    if "s2" in a.only:
        s2_panel(a)
    LOG.write(LOGS / f"04_fibroblast_genes{L.leave_out_suffix(a.leave_out_study)}.log")


if __name__ == "__main__":
    main()
