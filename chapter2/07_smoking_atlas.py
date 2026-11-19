#!/usr/bin/env python3
"""
Chapter 2, step 07: TXNRD1 detection and smoking in LUAD donors of the atlas (Supplementary Figure 2b).

WHAT: s1, per cell set (fibroblasts, macrophages, T and NK cells):
        - TXNRD1 detection of every cell of the working set (count > 0), and the detection expected
          from the cell's number of detected genes: per study, a binomial GLM of detection on a cubic
          polynomial of the standardized log genes detected, fitted on all cells of the cell set and
          study (LUAD and LUSC together); a linear term when the cubic fit fails or the study has fewer
          than 30 detecting cells, a constant when that fails too;
        - per LUAD donor with at least 20 cells: observed over expected detections (O/E) and the raw
          rate; smoking from the atlas obs (yes, ex, current, former, active, "ex; exposed to asbestos"
          -> ever; no -> never; unknown or empty -> left out); a donor whose cells disagree on study or
          smoking is left out;
        - ever against never smokers, van Elteren test stratified by study over the studies that hold
          both groups, with a permutation p (10,000 label shuffles within study) and the asymptotic p;
          the effect is the n1*n0/N-weighted mean of the within-study differences in median O/E with a
          2,000-resample bootstrap within study and group. The raw detection rate is tested the same way;
          that test is not printed and stays because it draws from the same random generators as the
          printed tests (removing it would change their permutation and bootstrap draws);
        - the smoking counts behind "nearly all donors with LUSC had smoked": donors per cell set,
          histology and smoking group.
      s2, the panel of Supplementary Figure 2b (one dot per donor, never and ever smokers, median bars,
          the permutation p), drawn from the tables of s1 re-read from disk.
ITEMS: Supplementary Figure 2b (file suppfig2d_smoking_atlas); the smoking sentence of the thioredoxin section
      and its counts (smoking_atlas_numbers.csv, smoking_atlas_counts.csv).
FROM: work_root/intermediate/chapter2/work_<cell set>/ (02_working_sets.py: meta_<histology>.npz, and
      raw_<histology>.npz, whose TXNRD1 entries are found by scanning the stored arrays memory-mapped,
      without loading them), and three obs columns (smoking, study, patient_id) of the atlas object
      (configuration key "atlas_object"), read with h5py; X is not touched.
TO:   work_root/output/chapter2/, every table with a sibling .md:
        s1  smoking_atlas_donors.csv    one row per LUAD donor and cell set (donor identifiers of the
                                        public atlas; stays under work_root)
            smoking_atlas_tests.csv     the six tests (three printed O/E tests, three raw-rate tests that
                                        are not printed)
            smoking_atlas_counts.csv    donors by cell set, histology and smoking group
            smoking_atlas_numbers.csv   the count behind "did not differ significantly"
        s2  suppfig2d_smoking_atlas.png/.pdf
      and the log work_root/logs/chapter2/07_smoking_atlas.log.
WHY:  the thioredoxin section asks whether TXNRD1 detection in stromal cells follows smoking. Detection
      depends strongly on the number of genes detected per cell, so each cell's expected detection is
      modelled study by study and the donors are compared on observed over expected. Statistics:
      van_elteren and strat_med_diff (redox.ranktests); panel layout 7.2 by 2.2 inches.
RUN:  python3 chapter2/07_smoking_atlas.py [--only s1 s2]
      The three cell sets always run together and in this order: the two generators are shared across
      them. Seeds: SEED_BOOT = 20260925 (bootstrap, 2,000 resamples), SEED_PERM = 20260926 (van Elteren
      permutations, 10,000), SEED_JITTER = 1 (horizontal jitter of the dots, drawing only). Needs
      statsmodels (the GLM).
"""
import argparse
import struct
import sys
import time
import zipfile
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import statsmodels.api as sm

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ch2lib as L  # noqa: E402
from redox import atlas_fixes, paths  # noqa: E402
from redox import plotting as RP  # noqa: E402  (chooses the Agg backend)
from redox.ranktests import strat_med_diff, van_elteren  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

SEED_BOOT, SEED_PERM, SEED_JITTER = 20260925, 20260926, 1
N_BOOT, N_PERM = 2000, 10000
MIN_CELLS = 20            # cells of the type per donor
MIN_DETECTING_CUBIC = 30  # detecting cells a study needs for the cubic fit
SMOKE = {"yes": "ever", "ex": "ever", "current": "ever", "former": "ever", "active": "ever",
         "ex; exposed to asbestos": "ever", "no": "never"}
HIST_ORDER = ("LUAD", "LUSC")   # the order in which the two histologies are stacked (fixes the order of the cells)
LOGS = paths.WORK_ROOT / "logs" / "chapter2"
AXIS = "#c3c2b7"
LOG = L.Log()


# ---------------------------------------------------------------- reading
def atlas_object():
    """The atlas object from the configuration (key "atlas_object")."""
    named = paths.CFG.get("atlas_object")
    if not named:
        raise SystemExit(f'{paths.CFG["_config_path"]} has no "atlas_object" entry (see README, Chapter 2)')
    p = Path(named).expanduser()
    if not p.exists():
        raise SystemExit(f"atlas object not found: {p}")
    return p


def npz_member(path, key):
    """One array of an .npz written by np.savez (members stored, not compressed), memory-mapped, so that
    the 1 to 2 GB count arrays of a working set are scanned from disk instead of read into memory."""
    with zipfile.ZipFile(path) as z:
        info = z.getinfo(f"{key}.npy")
    if info.compress_type != zipfile.ZIP_STORED:
        return np.load(path)[key]
    with open(path, "rb") as fh:
        fh.seek(info.header_offset)
        n_name, n_extra = struct.unpack("<HH", fh.read(30)[26:30])
        fh.seek(info.header_offset + 30 + n_name + n_extra)
        major, _ = np.lib.format.read_magic(fh)
        read = np.lib.format.read_array_header_1_0 if major == 1 else np.lib.format.read_array_header_2_0
        shape, fortran, dtype = read(fh)
        offset = fh.tell()
    return np.memmap(path, dtype=dtype, mode="r", shape=shape, offset=offset, order="F" if fortran else "C")


def txnrd1_detected(raw, gidx, chunk=50_000_000):
    """Per row of a working-set file: does the row hold a positive TXNRD1 count."""
    indptr = np.asarray(npz_member(raw, "indptr"))
    ind, dat = npz_member(raw, "indices"), npz_member(raw, "data")
    hits = []
    for s in range(0, len(ind), chunk):
        hits.append(np.nonzero(np.asarray(ind[s:s + chunk]) == gidx)[0] + s)
    pos = np.concatenate(hits)
    pos = pos[np.asarray(dat[pos]) > 0]
    row = np.searchsorted(indptr, pos, side="right") - 1
    det = np.zeros(len(indptr) - 1, bool)
    det[row] = True
    return det


def load_cells(cell, hist, gidx):
    """The cells of one cell set and histology: row in the object, patient, genes detected, TXNRD1 detected (0 or 1)."""
    W = L.work_dir(cell)
    m = np.load(W / f"meta_{hist}.npz", allow_pickle=True)
    raw = W / f"raw_{hist}.npz"
    return pd.DataFrame(dict(sub=np.asarray(npz_member(raw, "sub")), patient=m["patient"], ngenes=m["ngenes_obs"],
                             det=txnrd1_detected(raw, gidx).astype(int), hist=hist))


def obs_columns(names):
    """Whole obs columns of the atlas object as str arrays (h5py; categorical or plain strings)."""
    res = {}
    with h5py.File(atlas_object(), "r") as f:
        for k in names:
            x = f["obs"][k]
            if isinstance(x, h5py.Group):
                cats = x["categories"][:].astype(str)
                codes = x["codes"][:]
                res[k] = np.where(codes >= 0, cats[np.clip(codes, 0, None)], "")
            else:
                res[k] = x[:].astype(str)
    return res


def fit_expected(df):
    """Per study: binomial GLM of TXNRD1 detection on a cubic in standardized log genes detected; returns the fitted
    (expected) detection per cell and a note per study."""
    exp = np.full(len(df), np.nan)
    notes = []
    for st, idx in df.groupby("study").groups.items():
        s = df.loc[idx]
        z = np.log(s.ngenes.clip(lower=1))
        z = (z - z.mean()) / (z.std() + 1e-12)
        y = s.det.values
        kind = "cubic"
        X = np.column_stack([np.ones(len(z)), z, z ** 2, z ** 3])
        fitted = None
        if y.sum() >= MIN_DETECTING_CUBIC and y.sum() < len(y):
            try:
                fr = sm.GLM(y, X, family=sm.families.Binomial()).fit()
                if np.all(np.isfinite(fr.params)) and fr.converged:
                    fitted = fr.predict(X)
            except Exception:
                fitted = None
        if fitted is None:
            kind = "linear"
            try:
                fr = sm.GLM(y, X[:, :2], family=sm.families.Binomial()).fit()
                fitted = fr.predict(X[:, :2])
            except Exception:
                kind = "constant"
                fitted = np.full(len(y), y.mean())
        exp[df.index.get_indexer(idx)] = fitted
        notes.append(f"{st}: {len(s)} cells, {int(y.sum())} detecting, fit {kind}")
    return exp, notes


def donors(df, hist):
    """Per donor of one histology: cells, observed and expected detections, study, smoking."""
    g = df[df["hist"] == hist].groupby("patient")
    dd = pd.DataFrame(dict(n_cells=g.size(), O=g.det.sum(), E=g.expected.sum(), n_study=g.study.nunique(),
                           study=g.study.first(), n_smk=g.smk_raw.nunique(), smk_raw=g.smk_raw.first()))
    bad = dd[(dd.n_study > 1) | (dd.n_smk > 1)]
    dd = dd.drop(bad.index)
    dd["smoking"] = dd.smk_raw.map(SMOKE)
    dd["OE"] = dd.O / dd.E
    dd["raw_rate"] = dd.O / dd.n_cells
    return dd, len(bad)


# ====================================================================== s1
def s1_tests(a):
    """Observed over expected TXNRD1 detection per LUAD donor, and ever against never smokers (van Elteren test
    stratified by study; 10,000 label permutations within study, 2,000 bootstrap resamples)."""
    rng_b, rng_p = np.random.default_rng(SEED_BOOT), np.random.default_rng(SEED_PERM)
    genes = np.load(L.work_dir("fib") / "genes.npy", allow_pickle=True)
    gidx = int(np.where(genes == L.SPLITTER)[0][0])
    obs = obs_columns(["smoking", "study", "patient_id"])
    obs["study"], obs["patient_id"] = atlas_fixes.relabel(obs["study"]), atlas_fixes.relabel(obs["patient_id"])
    LOG(f"TXNRD1 column {gidx}; obs columns smoking, study, patient_id read   [{LOG.elapsed()}]")
    donor_rows, test_rows, count_rows = [], [], []
    for cell, spec in L.CELL_SETS.items():
        lab = spec["label"]
        g2 = np.load(L.work_dir(cell) / "genes.npy", allow_pickle=True)
        if str(g2[gidx]) != L.SPLITTER:
            raise SystemExit(f"TXNRD1 is not column {gidx} of work_{cell}/genes.npy")
        df = pd.concat([load_cells(cell, h, gidx) for h in HIST_ORDER], ignore_index=True)
        if not (obs["patient_id"][df["sub"].values] == df.patient.astype(str).values).all():
            raise SystemExit(f"work_{cell}: the patient of a cell differs from the atlas obs")
        df["study"] = obs["study"][df["sub"].values]
        df["smk_raw"] = obs["smoking"][df["sub"].values]
        df["expected"], notes = fit_expected(df)
        LOG(f"\n{lab}: {len(df):,} cells ({(df['hist'] == 'LUAD').sum():,} LUAD); detection {df.det.mean():.3f}"
            f"   [{LOG.elapsed()}]")
        for n in notes:
            LOG("  " + n)
        for hist in ("LUAD", "LUSC"):
            dh, nbad = donors(df, hist)
            for floor in (0, MIN_CELLS):
                d = dh[dh.n_cells >= floor]
                count_rows.append(dict(cell_set=cell, cell_type=lab, histology=hist, min_cells=floor,
                                       donors=len(d), ever=int((d.smoking == "ever").sum()),
                                       never=int((d.smoking == "never").sum()), unknown=int(d.smoking.isna().sum()),
                                       left_out_inconsistent=nbad))
        dd, nbad = donors(df, "LUAD")
        if nbad:
            LOG(f"  excluded for inconsistent study or smoking within donor: {nbad}")
        dd["cell_set"], dd["celltype"] = cell, lab
        LOG(f"  LUAD donors: {len(dd)}; with >= {MIN_CELLS} cells: {(dd.n_cells >= MIN_CELLS).sum()}; "
            f"smoking {dd.smk_raw.value_counts().to_dict()}")
        donor_rows.append(dd.reset_index())
        D = dd[(dd.n_cells >= MIN_CELLS) & dd.smoking.notna()]
        both = [s for s, x in D.groupby("study") if (x.smoking == "ever").any() and (x.smoking == "never").any()]
        Dt = D[D.study.isin(both)]
        LOG(f"  with known smoking: {len(D)} donors (ever {int((D.smoking == 'ever').sum())}, "
            f"never {int((D.smoking == 'never').sum())}); {len(both)} studies with both; "
            f"in the test ever {int((Dt.smoking == 'ever').sum())}, never {int((Dt.smoking == 'never').sum())}")
        # raw_rate: tests the raw TXNRD1 detection rate; it is not printed and stays because it draws from the same
        # random generators as the printed tests (removing it would change their permutation and bootstrap draws)
        for var, fam in (("OE", "B"), ("raw_rate", "sensitivity")):
            strata = [(Dt.loc[(Dt.study == s) & (Dt.smoking == "ever"), var].values,
                       Dt.loc[(Dt.study == s) & (Dt.smoking == "never"), var].values) for s in both]
            ve = van_elteren(strata, N_PERM, rng_p)
            est, lo, hi = strat_med_diff(strata, N_BOOT, rng_b)
            what = "O/E (depth-adjusted)" if var == "OE" else "raw rate"
            test_rows.append(dict(block="B2", family=fam, kind="main" if fam == "B" else "sensitivity",
                                  cohort="atlas", histology="LUAD", variable="Smoking ever vs never",
                                  outcome=f"TXNRD1 detection {what}: {lab}",
                                  test="van Elteren by study (permutation p)",
                                  n1=sum(len(x) for x, _ in strata), n0=sum(len(y) for _, y in strata),
                                  group1="ever", group0="never",
                                  median1=np.median(Dt.loc[Dt.smoking == "ever", var]),
                                  median0=np.median(Dt.loc[Dt.smoking == "never", var]),
                                  estimate=est, ci_lo=lo, ci_hi=hi, stat=ve["Z"], p=ve["p_perm"],
                                  p_asym=ve["p_asym"], note=f"{ve['n_strata']} studies",
                                  cell_set=cell, n_studies=ve["n_strata"]))
            LOG(f"  {var:8s} ever {test_rows[-1]['n1']} never {test_rows[-1]['n0']}: Z {ve['Z']:+.3f}, "
                f"permutation p {ve['p_perm']:.4f}, estimate {est:+.3f} ({lo:+.3f} to {hi:+.3f})")

    D = pd.concat(donor_rows, ignore_index=True)
    o = L.out("smoking_atlas_donors.csv")
    D.to_csv(o, index=False)
    L.sidecar(o, "LUAD donors of the atlas working sets, per cell set: TXNRD1 detection observed and expected "
              "from the genes detected per cell (study-wise GLM), and smoking.", [
                  ("patient", "donor (patient_id of the public atlas)"), ("n_cells", "cells of the type"),
                  ("O", "cells that detected TXNRD1"), ("E", "detections expected from each cell's genes detected"),
                  ("n_study", "studies of the donor's cells (1 after the exclusion)"), ("study", "study"),
                  ("n_smk", "smoking values of the donor's cells (1 after the exclusion)"),
                  ("smk_raw", "smoking as recorded in the atlas"), ("smoking", "ever, never or empty (unknown)"),
                  ("OE", "O / E"), ("raw_rate", "O / n_cells"), ("cell_set", "fib, mac or tnk"),
                  ("celltype", "cell set as printed")])
    T = pd.DataFrame(test_rows)
    o = L.out("smoking_atlas_tests.csv")
    T.to_csv(o, index=False)
    L.sidecar(o, "Ever against never smokers, LUAD donors with at least 20 cells, studies with both groups "
              "(the atlas part of the driver analysis, block B2).", [
                  ("block", "B2"), ("family", "B (printed tests) or sensitivity (raw rate, not printed)"),
                  ("kind", "main or sensitivity"), ("cohort", "atlas"), ("histology", "LUAD"),
                  ("variable", "the comparison"), ("outcome", "what is compared"), ("test", "the test"),
                  ("n1", "ever-smoker donors in the test"), ("n0", "never-smoker donors in the test"),
                  ("group1", "ever"), ("group0", "never"), ("median1", "median over ever-smoker donors"),
                  ("median0", "median over never-smoker donors"),
                  ("estimate", "weighted mean of within-study differences in medians, ever minus never"),
                  ("ci_lo", "bootstrap 2.5th percentile"), ("ci_hi", "bootstrap 97.5th percentile"),
                  ("stat", "van Elteren Z"), ("p", "permutation p (printed)"), ("p_asym", "asymptotic p"),
                  ("note", "studies with both groups"), ("cell_set", "fib, mac or tnk"),
                  ("n_studies", "studies with both groups (the strata)")])
    main = T[T.family == "B"]
    N2 = pd.DataFrame([("tests_p_below_005", int((main.p < 0.05).sum()),
                        "printed tests (O/E, three cell sets) with permutation p < 0.05 ('did not differ significantly')")],
                      columns=["item", "value", "meaning"])
    o = L.out("smoking_atlas_numbers.csv")
    N2.to_csv(o, index=False)
    L.sidecar(o, "Counts behind the smoking sentence of the thioredoxin section.",
              [("item", "key"), ("value", "value"), ("meaning", "what it counts")])
    C = pd.DataFrame(count_rows)
    o = L.out("smoking_atlas_counts.csv")
    C.to_csv(o, index=False)
    L.sidecar(o, "Donors of the working sets by cell set, histology and smoking group, all donors (min_cells 0) "
              f"and donors with at least {MIN_CELLS} cells ('nearly all donors with LUSC had smoked').", [
                  ("cell_set", "fib, mac or tnk"), ("cell_type", "cell set as printed"), ("histology", "LUAD or LUSC"),
                  ("min_cells", "least cells of the type per donor"), ("donors", "donors"),
                  ("ever", "ever smokers"), ("never", "never smokers"), ("unknown", "smoking not recorded"),
                  ("left_out_inconsistent", "donors left out for inconsistent study or smoking")])
    LOG(f"\nwrote smoking_atlas_donors.csv, smoking_atlas_tests.csv, smoking_atlas_numbers.csv, "
        f"smoking_atlas_counts.csv   [{LOG.elapsed()}]")


# ====================================================================== s2
def s2_panel(a):
    """Supplementary Figure 2b: one dot per donor (O/E), never and ever smokers, median bars and the permutation p."""
    RP.style()
    plt.rcParams["axes.edgecolor"] = AXIS
    D = pd.read_csv(L.out("smoking_atlas_donors.csv"), float_precision="round_trip", dtype={"patient": str})
    Tt = pd.read_csv(L.out("smoking_atlas_tests.csv"), float_precision="round_trip").query("kind == 'main'")
    D = D[D.smoking.isin(["ever", "never"]) & (D.n_cells >= MIN_CELLS)]
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.2))
    fig.subplots_adjust(left=0.08, right=0.99, bottom=0.17, top=0.72, wspace=0.35)
    rng = np.random.default_rng(SEED_JITTER)
    for ax, (cell, spec) in zip(axs, L.CELL_SETS.items()):
        g = D[D.cell_set == cell]
        both = [s for s, x in g.groupby("study") if x.smoking.nunique() == 2]
        g = g[g.study.isin(both)]
        t = Tt[Tt.cell_set == cell].iloc[0]
        if (int(t.n1), int(t.n0)) != (int((g.smoking == "ever").sum()), int((g.smoking == "never").sum())):
            raise SystemExit(f"{cell}: the donors drawn do not match the test table")
        for k, grp in enumerate(["never", "ever"]):
            v = g.loc[g.smoking == grp, "OE"].to_numpy()
            ax.scatter(k + rng.uniform(-0.16, 0.16, len(v)), v, s=9, facecolor="white", edgecolor=RP.INK2, lw=0.6,
                       zorder=2)
            ax.plot([k - 0.25, k + 0.25], [np.median(v)] * 2, color=RP.INK, lw=1.4, zorder=3)
        ax.axhline(1, color=AXIS, lw=0.5, ls=":")
        ax.set_xticks([0, 1])
        ax.set_xticklabels([f"never\n({int(t.n0)})", f"ever\n({int(t.n1)})"], fontsize=6.4)
        ax.set_xlim(-0.6, 1.6)
        ax.set_title(f"{spec['label']}\np = {t.p:.3f}", fontsize=7, pad=3)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(labelsize=6, length=2)
    axs[0].set_ylabel("TXNRD1-positive cells,\nobserved / expected", fontsize=6.6)
    fig.text(0.03, 0.965, "TXNRD1 expression and smoking, LUAD patients from studies with smokers and never-smokers",
             fontsize=7.2, fontweight="bold", va="top")
    RP.save(fig, 2, "suppfig2d_smoking_atlas")
    LOG("wrote suppfig2d_smoking_atlas.png and .pdf")


def main():
    """Run s1 and/or s2 and write the log."""
    ap = argparse.ArgumentParser(description=__doc__.split("RUN:")[0].strip(),
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", nargs="*", default=["s1", "s2"], choices=["s1", "s2"])
    a = ap.parse_args()
    LOGS.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    if "s1" in a.only:
        s1_tests(a)
    if "s2" in a.only:
        s2_panel(a)
    LOG.write(LOGS / "07_smoking_atlas.log")
    print(f"07_smoking_atlas.py done in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
