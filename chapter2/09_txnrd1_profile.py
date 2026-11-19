#!/usr/bin/env python3
"""
Chapter 2, step 09: TXNRD1 across the stromal cell types of the atlas (Figure 2a).

WHAT: s1, one pass over the count matrix of the atlas object: for every cell of a primary LUSC or LUAD
          tumor in a block that passed QC (the cells of 02_working_sets.py, all compartments), its TXNRD1
          count and its total count. The object holds three compartments, T and NK cells, myeloid cells and
          fibroblasts; no cancer, endothelial or B cells.
      s2, per patient with at least 20 cells of a type:
        - the fraction of the patient's cells of the type with a positive TXNRD1 count (detection);
        - the total counts per cell, the mean over the patient's cells of the cell's total count (the
          count depth behind "lowest in neutrophils and in T and NK cells");
        for two sets of types: the nine types of the panel (CD4 T cells, CD8 T cells, regulatory T cells,
        NK cells; fibroblasts; macrophages, monocytes, dendritic cells, neutrophils, from level2_celltype;
        plasmacytoid dendritic cells, other T and other myeloid cells and unlabeled cells are left out),
        and the three cell sets of Figure 3 (T and NK cells, fibroblasts, macrophages). The per-patient
        values of the three cell sets are checked against the working sets of 02_working_sets.py.
        Summary per type of the panel and histology (patients, cells, median and quartiles over patients of
        detection and of total counts, and the rank of the median total counts among the nine types).
      s3, Figure 2a: horizontal boxplots of detection over patients (box, interquartile range; line,
          median; whiskers, the most extreme patients within 1.5 interquartile ranges of the box; patients beyond
          them are not drawn); colors are the lineage colors of redox.plotting in shades per subtype. Title and
          patient counts (from the summary) are written on the panel.
ITEMS: Figure 2a (file suppfig2a_txnrd1_profile) and the numbers of the Results sentence on TXNRD1 detection
      by cell type (txnrd1_profile_summary.csv). The leave-out run that backs the paper is
      --leave-out-study Jimenez_Bruno_2026 (the ablation study GSE325414; Methods sentence on it), with --only s2 s3.
FROM: the atlas object named by the configuration key "atlas_object" (X/indices and X/data streamed in
      windows, obs columns read with h5py); for the check, work_root/intermediate/chapter2/work_<cell set>/.
TO:   work_root/intermediate/chapter2/txnrd1_percell.npz (per cell: row, patient, histology, level 1 and
      level 2 labels, TXNRD1 count, total count), and work_root/output/chapter2/, every table with a .md:
        s2  txnrd1_profile_patients.csv   one row per patient and type (public atlas donor identifiers;
                                          stays under work_root)
            txnrd1_profile_summary.csv    per type of the panel and histology: patients, cells, median, quartiles
        s3  suppfig2a_txnrd1_profile.png/.pdf   (Figure 2a)
      and the log work_root/logs/chapter2/09_txnrd1_profile.log.
WHY:  Figure 2 asks what TXNRD1 marks within each cell type (panel c) and first how often each stromal cell
      type detects it (panel a); detection is shown per patient, so that a patient with many cells does not
      dominate.
RUN:  python3 chapter2/09_txnrd1_profile.py [--only s1 s2 s3] [--window 64000000] [--leave-out-study PREFIX ...]
      s1 reads the whole count matrix once (about 22 GB); s2 and s3 read its output. No random numbers.
      --leave-out-study (with --only s2 s3) leaves out the patients of one study or several (patient
      identifiers starting with the study label, which is also the prefix of its batch names, as in
      03_program_tests.py) and writes the tables, the panel and the log with the suffix _without_<prefix> beside
      the ordinary ones, which it leaves untouched; txnrd1_percell.npz is read, not rewritten.
"""
import argparse
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ch2lib as L  # noqa: E402
from redox import atlas_fixes, paths  # noqa: E402
from redox import plotting as RP  # noqa: E402  (chooses the Agg backend)
import matplotlib.pyplot as plt  # noqa: E402

import importlib.util  # noqa: E402
_spec = importlib.util.spec_from_file_location("smoking07", Path(__file__).resolve().parent / "07_smoking_atlas.py")
S07 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(S07)    # atlas_object and npz_member (the memory-mapped reader of the working sets)

MIN_CELLS = 20                   # cells of the type per patient
HISTS = ["LUSC", "LUAD"]
COMPARTMENTS = ["T_NK", "fibroblast", "myeloid"]
# The nine types of the panel, top to bottom: key, label, level 1, level 2 (None: the whole compartment),
# lineage color of redox.plotting, shade (negative darker, positive lighter).
TYPES = [("cd4", "CD4 T cells", "T_NK", "CD4_T", "T", -0.25),
         ("cd8", "CD8 T cells", "T_NK", "CD8_T", "T", 0.0),
         ("treg", "regulatory T cells", "T_NK", "Treg", "T", 0.3),
         ("nk", "NK cells", "T_NK", "NK", "T", 0.55),
         ("fib_all", "fibroblasts", "fibroblast", None, "Fib", 0.0),
         ("macro", "macrophages", "myeloid", "macrophage", "Mac", 0.0),
         ("mono", "monocytes", "myeloid", "monocyte", "Mac", 0.25),
         ("dc", "dendritic cells", "myeloid", "DC", "Mac", 0.45),
         ("neut", "neutrophils", "myeloid", "neutrophil", "Mac", 0.62)]
# The three cell sets of Figure 3 (ch2lib.CELL_SETS), for the text and the check against the working sets.
SETS = [("tnk", "T and NK cells", "T_NK", None), ("fib", "fibroblasts", "fibroblast", None),
        ("mac", "macrophages", "myeloid", "macrophage")]
LOGS = paths.WORK_ROOT / "logs" / "chapter2"
AXIS = "#c3c2b7"
TITLE = "TXNRD1 detection by stromal cell type, primary LUSC and LUAD tumors"
LOG = L.Log()


def percell_path():
    """The per-cell file written by s1 (work_root/intermediate/chapter2/txnrd1_percell.npz)."""
    return paths.intermediate("chapter2") / "txnrd1_percell.npz"


def s1_stream(a):
    """TXNRD1 count and total count of every selected cell, in one pass over X."""
    LOG("09_txnrd1_profile.py, s1: TXNRD1 and total counts per cell")
    with h5py.File(S07.atlas_object(), "r") as f:
        def col(k):
            """One obs column of the atlas object as a str array."""
            return f["obs"][k][:].astype(str)
        gk = f["var"].attrs["_index"]
        genes = f[f"var/{gk.decode() if isinstance(gk, bytes) else gk}"][:].astype(str)
        gidx = int(np.nonzero(genes == L.SPLITTER)[0][0])
        l1, hist = col("level1_compartment"), col("histology")
        keep = (col("tissue_class") == "primary_tumour") & f["obs/block_qc_pass"][:].astype(bool)
        keep &= np.isin(hist, HISTS) & np.isin(l1, COMPARTMENTS)
        keep &= ~atlas_fixes.duplicate(col("patient_id"))     # the redistributed Lambrechts 2018 copies
        rows = np.nonzero(keep)[0]
        LOG(f"  {len(rows):,} cells selected; TXNRD1 is gene column {gidx}")
        indptr = f["X/indptr"][:].astype(np.int64)
        nnz = int(indptr[-1])
        tx = np.zeros(len(indptr) - 1)
        tot = np.zeros(len(indptr) - 1)
        Xi, Xd = f["X/indices"], f["X/data"]
        t0 = time.time()
        for w0 in range(0, nnz, a.window):
            w1 = min(w0 + a.window, nnz)
            ind = Xi[w0:w1]
            dat = Xd[w0:w1].astype(np.float64)
            hit = np.nonzero(ind == gidx)[0]
            tx[np.searchsorted(indptr, hit + w0, side="right") - 1] = dat[hit]
            r0 = int(np.searchsorted(indptr, w0, side="right") - 1)
            r1 = int(np.searchsorted(indptr, w1 - 1, side="right") - 1)
            b = np.clip(indptr[r0:r1 + 2], w0, w1) - w0
            cs = np.concatenate([[0.0], np.cumsum(dat)])
            tot[r0:r1 + 1] += cs[b[1:]] - cs[b[:-1]]
        LOG(f"  streamed {nnz:,} stored values in {time.time() - t0:.0f}s")
        ncounts = f["obs/n_counts"][:][rows]
        np.savez(percell_path(), row=rows, patient=atlas_fixes.relabel(col("patient_id")[rows]), hist=hist[rows], l1=l1[rows],
                 l2=col("level2_celltype")[rows], tx=tx[rows], tot=tot[rows])
    gap = np.abs(ncounts - tot[rows]).max()
    LOG(f"  total counts against obs n_counts: largest difference {gap:.3f} (non-integer estimated counts)")
    if gap > 1:
        raise SystemExit("total counts disagree with obs n_counts")
    LOG(f"  -> {percell_path()}")


def per_patient(C, key, l1, l2):
    """Per patient and histology with at least MIN_CELLS cells of the type: cells, detection fraction, mean total counts."""
    s = (C.l1 == l1) if l2 is None else ((C.l1 == l1) & (C.l2 == l2))
    g = C[s].groupby(["hist", "patient"])
    d = pd.DataFrame(dict(n_cells=g.size(), detection=g.det.mean(), tot=g.tot.mean()))
    d = d[d.n_cells >= MIN_CELLS].reset_index().rename(columns={"hist": "histology"})
    d.insert(0, "type", key)
    return d


def check_sets(C, P, leave_out=None):
    """The per-patient values of the three cell sets equal those computed from the working sets (without the
    left-out patients when a study is left out)."""
    for key, _, _, _ in SETS:
        W = L.work_dir(key)
        rows = []
        for h in HISTS:
            m = np.load(W / f"meta_{h}.npz", allow_pickle=True)
            raw = W / f"raw_{h}.npz"
            sub = np.asarray(S07.npz_member(raw, "sub"))
            rows.append(pd.DataFrame(dict(row=sub, patient=m["patient"].astype(str), hist=h)))
        R = pd.concat(rows)
        R = R[~L.left_out(R.patient, leave_out)]
        mine = C[C.row.isin(R.row)]
        if len(mine) != len(R) or set(mine.row) != set(R.row):
            raise SystemExit(f"{key}: the streamed cells are not the cells of work_{key}")
        wanted = P[P.type == key].set_index(["histology", "patient"]).n_cells
        again = R.groupby(["hist", "patient"]).size()
        again = again[again >= MIN_CELLS]
        again.index.names = ["histology", "patient"]
        if not wanted.sort_index().equals(again.sort_index()):
            raise SystemExit(f"{key}: patients or cell counts differ from work_{key}")
    LOG("  check: the three cell sets hold exactly the cells and patients of the working sets")


def s2_tables(a):
    """The per-patient table and the summary per type and histology (patients with at least 20 cells of the type)."""
    LOG("09_txnrd1_profile.py, s2: per patient and type")
    sfx = L.leave_out_suffix(a.leave_out_study)
    z = np.load(percell_path(), allow_pickle=True)
    C = pd.DataFrame({k: z[k] for k in ("row", "patient", "hist", "l1", "l2", "tx", "tot")})
    C["patient"] = C.patient.astype(str)
    if sfx:
        out_ = L.left_out(C.patient, a.leave_out_study)
        LOG(f"  leaving out {' and '.join(L.leave_out_prefixes(a.leave_out_study))}: {int(out_.sum()):,} cells of "
            f"{C.patient[out_].nunique()} patients")
        C = C[~out_].reset_index(drop=True)
    C["det"] = (C.tx > 0).astype(float)
    P = pd.concat([per_patient(C, k, l1, l2) for k, _, l1, l2, _, _ in TYPES]
                  + [per_patient(C, k, l1, l2) for k, _, l1, l2 in SETS], ignore_index=True)
    check_sets(C, P, a.leave_out_study)
    labels = {k: lab for k, lab, *_ in TYPES}
    summary = []
    for key in labels:
        for h, hh in (("LUSC", ["LUSC"]), ("LUAD", ["LUAD"]), ("both", HISTS)):
            x = P[(P.type == key) & P.histology.isin(hh)]
            row = dict(type=key, label=labels[key], histology=h, n=len(x),
                       cells=int((C["hist"].isin(hh) & _mask(C, key)).sum()))
            for k in ("detection", "tot"):
                q1, med, q3 = np.quantile(x[k], [.25, .5, .75])
                row.update({f"{k}_median": med, f"{k}_q1": q1, f"{k}_q3": q3})
            summary.append(row)
        LOG(f"  {labels[key]:20s} {summary[-1]['n']:4d} patients, detection {summary[-1]['detection_median']:.3f}, "
            f"total counts {summary[-1]['tot_median']:,.0f}")
    summary = pd.DataFrame(summary)
    summary["tot_rank"] = summary.groupby("histology").tot_median.rank(method="min").astype(int)
    nine = ("cd4, cd8, treg, nk, fib_all, macro, mono, dc, neut: the nine types of the panel (level2_celltype; "
            "fib_all the whole fibroblast compartment)")
    kinds = nine + "; tnk, fib, mac: the three cell sets of the tests (checked against the working sets)"
    write(P, f"txnrd1_profile_patients{sfx}",
          "One row per atlas patient and type with at least 20 cells of the type: TXNRD1 detection and total "
          "counts. Written by chapter2/09_txnrd1_profile.py s2.",
          [("type", kinds), ("histology", "LUSC or LUAD"), ("patient", "atlas patient_id"),
           ("n_cells", "cells of the type of the patient"),
           ("detection", "share of the cells with a positive TXNRD1 count"),
           ("tot", "mean over the cells of the total count (all genes) of a cell")])
    write(summary, f"txnrd1_profile_summary{sfx}",
          "Per type of the panel and histology: patients, cells and the median and quartiles over patients of the "
          "two measures (one value per patient with at least 20 cells of the type). Written by "
          "chapter2/09_txnrd1_profile.py s2.",
          [("type", nine), ("label", "the name shown"), ("histology", "LUSC, LUAD or both"),
           ("n", "patients with at least 20 cells of the type"), ("cells", "cells of the type, all patients"),
           ("detection_median", "median detection"), ("detection_q1", "25th percentile"),
           ("detection_q3", "75th percentile"),
           ("tot_median", "median total counts per cell (counts, per patient the mean over the cells)"),
           ("tot_q1", "25th percentile"), ("tot_q3", "75th percentile"),
           ("tot_rank", "rank of tot_median among the nine types of the same histology row, 1 = fewest counts")])


def _mask(c, key):
    """Boolean mask of the cells of a type (key of TYPES or SETS)."""
    spec = {k: (l1, l2) for k, _, l1, l2, *_ in TYPES} | {k: (l1, l2) for k, _, l1, l2 in SETS}
    l1, l2 = spec[key]
    return (c.l1 == l1) if l2 is None else ((c.l1 == l1) & (c.l2 == l2))


def write(frame, name, what, columns):
    """Write a table to the chapter 2 output folder with its sidecar; every column must be described."""
    path = L.out(f"{name}.csv")
    missing = [c for c in frame.columns if c not in dict(columns)]
    if missing:
        raise ValueError(f"{name}: the sidecar does not describe {missing}")
    frame.to_csv(path, index=False)
    L.sidecar(path, what, columns)
    LOG(f"  -> {path.name} ({len(frame)} rows)")


def shade(hex_color, f):
    """A lighter (f > 0, toward white) or darker (f < 0, toward black) shade of a color."""
    c = np.array([int(hex_color[i:i + 2], 16) for i in (1, 3, 5)], float)
    c = c + (255 - c) * f if f >= 0 else c * (1 + f)
    return "#" + "".join(f"{int(round(v)):02x}" for v in c)


def s3_panel(a):
    """Figure 2a: TXNRD1 detection per patient, boxplots, with the title and the patients per type; the total
    counts are in the summary table and the legend."""
    RP.style()
    sfx = L.leave_out_suffix(a.leave_out_study)
    P = pd.read_csv(L.out(f"txnrd1_profile_patients{sfx}.csv"), dtype={"patient": str}, float_precision="round_trip")
    Sm = pd.read_csv(L.out(f"txnrd1_profile_summary{sfx}.csv")).query("histology == 'both'").set_index("type")
    ys, y = [], 0.0
    for k, (key, _, l1, *_r) in enumerate(TYPES):
        if k and l1 != TYPES[k - 1][2]:
            y += 0.55                                   # a gap between compartments
        ys.append(y)
        y += 1.0
    fig = plt.figure(figsize=(7.2, 2.05))
    fig.text(0.03, 0.975, TITLE, fontsize=7.4, fontweight="bold", va="top")
    ax = fig.add_axes([0.22, 0.17, 0.36, 0.70])
    data = [P.loc[P["type"] == key, "detection"].to_numpy() * 100 for key, *_ in TYPES]
    bp = ax.boxplot(data, positions=ys, vert=False, widths=0.62, showfliers=False, patch_artist=True,
                    medianprops=dict(color=RP.INK, lw=1.0), whiskerprops=dict(color=RP.INK2, lw=0.6),
                    capprops=dict(color=RP.INK2, lw=0.6), boxprops=dict(lw=0))
    for patch, (_, _, _, _, lin, f) in zip(bp["boxes"], TYPES):
        patch.set_facecolor(shade(RP.LC[lin], f))
    ax.set_ylim(ys[-1] + 0.6, -0.6)
    ax.set_xlim(0, None)
    ax.set_xlabel("cells with TXNRD1 detected, % per patient", fontsize=6.4, labelpad=1)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(AXIS)
    ax.tick_params(labelsize=6.2, length=2)
    ax.set_yticks(ys)
    ax.set_yticklabels([lab for _, lab, *_ in TYPES], fontsize=6.4)
    groups = (TYPES[:2], TYPES[2:4], TYPES[4:6], TYPES[6:])   # the lines of the patient counts (four T-cell types on two lines)
    fig.text(0.64, 0.80, "patients with at least 20 cells of the type\n" + "\n".join(
        ", ".join(f"{lab} {int(Sm.loc[key, 'n'])}" for key, lab, *_ in g) for g in groups),
        fontsize=5.9, color=RP.INK2, va="top", linespacing=1.4)
    RP.save(fig, 2, f"suppfig2a_txnrd1_profile{sfx}")
    LOG(f"  -> suppfig2a_txnrd1_profile{sfx}.png and .pdf")


def main():
    """Run the steps named with --only (default s1 to s3)."""
    ap = argparse.ArgumentParser(description=__doc__.split("RUN:")[0].strip(),
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", nargs="*", default=["s1", "s2", "s3"], choices=["s1", "s2", "s3"])
    ap.add_argument("--window", type=int, default=64_000_000)
    ap.add_argument("--leave-out-study", nargs="+", default=None, metavar="PREFIX",
                    help="study label (patient-identifier prefix) to leave out, or several; s2 and s3 only")
    a = ap.parse_args()
    if a.leave_out_study and "s1" in a.only:
        raise SystemExit("--leave-out-study runs s2 and s3 only: add --only s2 s3")
    LOGS.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    for s, fn in (("s1", s1_stream), ("s2", s2_tables), ("s3", s3_panel)):
        if s in a.only:
            fn(a)
    LOG.write(LOGS / f"09_txnrd1_profile{L.leave_out_suffix(a.leave_out_study)}.log")
    print(f"09_txnrd1_profile.py done in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
