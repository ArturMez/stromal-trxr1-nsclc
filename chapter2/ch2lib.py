"""
Shared pieces of the chapter 2 scripts (single-cell RNA atlas, TXNRD1 split).

WHAT: a library, not run by itself. It provides the three cell sets (fibroblasts, macrophages, T and
      NK cells) and the cancer-cell set, the location of their working sets, the depth-matched
      pairing of TXNRD1-detected and TXNRD1-negative cells (match_fast, with the reference matcher
      match_random_ties), the working-set loader (load_set), the program scorer (score_programs and
      Scorer), the Benjamini-Hochberg q, batch-to-study labels, and a small logger.
USED BY: 01 to 12 of this folder (03 and 06 use the matcher and the scorer; 02 uses the cell sets).
FROM: work_root/intermediate/chapter2/work_<cell set>/ (written by 02_working_sets.py) and
      chapter2/data/ (figure_panels.json).
TO:   nothing on disk; the numbered scripts write the tables and panels.
WHY:  every script must score the same cells in the same way, so the setup is written once here.

The matcher. Cells that detect TXNRD1 are paired one to one with cells that do not, within a stratum
(batch and histology), choosing the unexposed cell with the nearest number of detected genes (ties
broken at random). match_random_ties is the reference implementation. match_fast is a faster
implementation that draws the same random numbers in the same order and returns the same pairs; it
is checked against match_random_ties by check_matcher() (03 calls it at the start of every run). The
random seed is SEED = 6 and the detection tolerance of the comparison genes is DET_TOL = 0.20. SEED,
DET_TOL and match_random_ties reproduce the original analysis code exactly.
"""
import bisect
import json
import re
import sys
import time
import resource
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, vstack as spvstack

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from redox import paths  # noqa: E402

DATA = HERE / "data"

# ---------------------------------------------------------------- seed and tolerance
SEED = 6            # global seed: the observed pairing, the cap per batch and histology, and the permutation stream
DET_TOL = 0.20      # comparison genes are detected within 20 % of TXNRD1's own detection rate

# ---------------------------------------------------------------- thresholds of the working set and the scores
# PRESENCE: a gene must be measured in at least 60 % of the cells; DETFLOOR: detected in at least 2 % of them;
# MINPAIR: least pairs for an effect; MINBLK: least cells of a batch for its z score; MINBLKPAIR: least pairs of
# a batch; MINSET: least tested genes of a program; Z95: two-sided 5 % normal quantile.
PRESENCE, DETFLOOR, MINPAIR, MINBLK, MINBLKPAIR, MINSET, Z95 = 0.60, 0.02, 50, 10, 20, 3, 1.96
SPLITTER = "TXNRD1"
HISTS = ("LUSC", "LUAD")

# the three cell sets of the gene-program analysis (Figure 3, Supplementary Figures 3 and 4)
CELL_SETS = {
    "fib": dict(label="Fibroblasts", l1=["fibroblast"], l2=None, short="fibroblasts"),
    "mac": dict(label="Macrophages", l1=["myeloid"], l2="macrophage", short="macrophages"),
    "tnk": dict(label="T and NK cells", l1=["T_NK"], l2=None, short="T and NK cells"),
}

# The cancer-cell set (Supplementary Figure 2c). It is built and tested only when named (--cells can). "can": the
# epithelial cells of primary LUSC and LUAD tumors that the original authors labeled as cancer cells (LuCA "Cancer
# cells", E-MTAB-6149 "tumor"), read from the epithelial object built by atlas/06_cancer_cell_object.py (configuration
# key "epi_object"), which leaves the atlas object untouched.
EXTRA_CELL_SETS = {
    "can": dict(label="Cancer cells", l1=["epithelial_malignant"], l2=None, l3_native=["Cancer cells", "tumor"],
                object="epi_object", short="cancer cells"),
}


def cell_spec(cell):
    """The definition of a cell set of CELL_SETS or EXTRA_CELL_SETS."""
    return CELL_SETS[cell] if cell in CELL_SETS else EXTRA_CELL_SETS[cell]


def work_dir(cell):
    """work_root/intermediate/chapter2/work_<cell set>, created on demand."""
    d = paths.intermediate("chapter2") / f"work_{cell}"
    d.mkdir(parents=True, exist_ok=True)
    return d


def out(name=None):
    """work_root/output/chapter2[/name]."""
    return paths.output(2, name)


def panels():
    """The figure panel configuration (chapter2/data/figure_panels.json)."""
    return json.loads((DATA / "figure_panels.json").read_text())


def stem(set_id):
    """The set id without its version suffix (the part before the first dot)."""
    return str(set_id).split(".")[0]


# ---------------------------------------------------------------- logging
class Log:
    """Print and keep every line; write() puts them in a file."""

    def __init__(self):
        self.lines, self.t0 = [], time.time()

    def __call__(self, s=""):
        print(s, flush=True)
        self.lines.append(str(s))

    def elapsed(self):
        """Seconds since the logger was made and the peak memory, as text."""
        return f"{time.time() - self.t0:.0f}s, peak {rss_gb():.2f} GB"

    def write(self, path):
        """Write all kept lines to a text file."""
        Path(path).write_text("\n".join(self.lines) + "\n")


def rss_gb():
    """Peak resident memory of this process in GB."""
    v = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return v / 1e9 if sys.platform == "darwin" else v / 1e6


def sidecar(path, what, columns):
    """The sibling .md of an output table: what it holds, one line per column."""
    path = Path(path)
    lines = [f"# {path.name}", "", what, "", "| column | meaning |", "|---|---|"]
    lines += [f"| {c} | {m} |" for c, m in columns]
    path.with_suffix(".md").write_text("\n".join(lines) + "\n")


# ---------------------------------------------------------------- the matcher
def match_random_ties(depth, pos, strat, r):
    """within-stratum nearest-depth 1:1 matching of exposed to unexposed cells, ties broken at
    random. depth is genes detected per cell; strat is the stratum label per cell.
    This is the reference matcher; match_fast gives the same pairs faster."""
    a, b = [], []
    for st in np.unique(strat):
        m = strat == st
        p = np.where(m & pos)[0]
        q = np.where(m & ~pos)[0]
        if len(p) == 0 or len(q) == 0:
            continue
        jit = r.random(len(q)) * 1e-6
        dq = depth[q] + jit
        o = np.argsort(dq, kind="stable")
        qs, dqs = q[o], dq[o]
        used = np.zeros(len(qs), bool)
        for i in p[r.permutation(len(p))]:
            free = np.where(~used)[0]
            if not len(free):
                break
            j = free[np.argmin(np.abs(dqs[free] - depth[i]))]
            used[j] = True
            a.append(i)
            b.append(qs[j])
    return np.asarray(a, int), np.asarray(b, int)


class Strata:
    """The stratum membership of every cell, computed once, in the order match_random_ties visits it."""

    def __init__(self, strat):
        self.labels, inv = np.unique(strat, return_inverse=True)
        order = np.argsort(inv, kind="stable")
        bounds = np.searchsorted(inv[order], np.arange(len(self.labels) + 1))
        self.members = [order[bounds[k]:bounds[k + 1]] for k in range(len(self.labels))]


def match_fast(depth, pos, strata, r):
    """The same pairs as match_random_ties, drawing the same random numbers in the same order.

    For each exposed cell (in the random order) match_random_ties takes, among the free unexposed
    cells in sorted-depth order, the first one with the smallest |depth difference|. Along the
    sorted axis the distance falls to the exposed cell's depth and rises after it, so that cell is
    the nearest free one on the left or the right; on a tie the left one, and on the left the
    leftmost of equal distances (np.argmin returns the lowest index). The nearest free cell on
    either side is found with two path-compressed pointer arrays instead of a scan of every free
    cell, which is what made the permutations slow.
    """
    A, B = [], []
    for mem in strata.members:
        pm = pos[mem]
        p, q = mem[pm], mem[~pm]
        if len(p) == 0 or len(q) == 0:
            continue
        jit = r.random(len(q)) * 1e-6
        dq = depth[q] + jit
        o = np.argsort(dq, kind="stable")
        qs, dqs = q[o], dq[o]
        n = len(qs)
        dl = dqs.tolist()
        nxt = list(range(n + 1))          # nxt[k] == k: index k free (k == n: none to the right)
        prv = list(range(n + 1))          # prv[k + 1] == k + 1: index k free (prv[0]: none to the left)

        def find(par, k):
            root = k
            while par[root] != root:
                root = par[root]
            while par[k] != root:
                par[k], k = root, par[k]
            return root

        nfree = n
        for i in p[r.permutation(len(p))]:
            if nfree == 0:
                break
            d = float(depth[i])
            k = bisect.bisect_left(dl, d)
            jr = find(nxt, k)                     # smallest free index >= k, n when none
            jl = find(prv, k) - 1                 # largest free index <= k - 1, -1 when none
            if jl >= 0:
                dist = abs(dl[jl] - d)
                while jl > 0:
                    j2 = find(prv, jl) - 1        # next free index to the left
                    if j2 >= 0 and abs(dl[j2] - d) == dist:
                        jl = j2
                    else:
                        break
            if jl < 0:
                j = jr
            elif jr >= n:
                j = jl
            else:
                j = jl if abs(dl[jl] - d) <= abs(dl[jr] - d) else jr
            nxt[j] = j + 1
            prv[j + 1] = j
            nfree -= 1
            A.append(i)
            B.append(qs[j])
    return np.asarray(A, int), np.asarray(B, int)


def check_matcher(depth, pos, strat, n=3, seed=SEED):
    """match_fast against match_random_ties on n successive draws of one generator each; raises if
    any pair differs. Called by 03 before the permutations start."""
    st = Strata(strat)
    r1, r2 = np.random.default_rng(seed), np.random.default_rng(seed)
    for _ in range(n):
        a1, b1 = match_random_ties(depth, pos, strat, r1)
        a2, b2 = match_fast(depth, pos, st, r2)
        if not (np.array_equal(a1, a2) and np.array_equal(b1, b2)):
            raise SystemExit("match_fast differs from match_random_ties; run with --slow-matcher")
    if r1.random() != r2.random():
        raise SystemExit("match_fast consumed a different number of random draws")
    return True


# ---------------------------------------------------------------- statistics
def bh(p):
    """Benjamini-Hochberg q, missing values kept missing."""
    p = np.asarray(p, float)
    q = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return q
    pv = p[ok]
    o = np.argsort(pv)
    m = len(pv)
    qq = pv[o] * m / np.arange(1, m + 1)
    qq = np.minimum.accumulate(qq[::-1])[::-1]
    res = np.empty(m)
    res[o] = np.minimum(qq, 1.0)
    q[ok] = res
    return q


def bh_all(p):
    """Benjamini-Hochberg q for a vector without missing values, clipped to [0, 1]."""
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    q = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    res = np.empty(n)
    res[o] = np.clip(q, 0, 1)
    return res


# ---------------------------------------------------------------- working set
def load_meta(W):
    """The metadata arrays of the LUSC and LUAD working sets (meta_LUSC.npz, meta_LUAD.npz)."""
    mA = np.load(W / "meta_LUSC.npz", allow_pickle=True)
    mB = np.load(W / "meta_LUAD.npz", allow_pickle=True)
    return mA, mB


def load_rows(W, h, rows):
    """The count rows of one histology (raw_<histology>.npz, sparse) for the given row numbers."""
    r = np.load(W / f"raw_{h}.npz", allow_pickle=True)
    return csr_matrix((r["data"], r["indices"], r["indptr"]),
                      shape=(int(r["n"][0]), int(r["ngenes"][0])))[rows]


def leave_out_prefixes(leave_out):
    """The left-out study prefixes as a tuple: None -> (), one string -> (it,), a list -> the list."""
    if leave_out is None:
        return ()
    return (leave_out,) if isinstance(leave_out, str) else tuple(leave_out)


def left_out(ids, leave_out):
    """True for the identifiers (batch or patient names) that start with one of the left-out prefixes."""
    ids = np.asarray(ids).astype(str)
    hit = np.zeros(len(ids), bool)
    for p in leave_out_prefixes(leave_out):
        hit |= np.char.startswith(ids, p)
    return hit


def leave_out_suffix(leave_out):
    """The file-name suffix of a leave-out run: "" for none, "_without_<prefix>" for one study,
    "_without_<prefix 1>_and_<prefix 2>" for several."""
    p = leave_out_prefixes(leave_out)
    return f"_without_{'_and_'.join(p)}" if p else ""


def load_set(cell, cap=2000, log=print, leave_out=None, histology=None):
    """The working set of one cell set: the cap per batch and histology (at most `cap` cells, drawn with
    seed SEED = 6), counts, TXNRD1 exposure, arm presence and the tested gene universe. Returns a dict.

    leave_out: a batch-name prefix (a study label) whose cells are removed before the cap, for the
    leave-one-study-out runs of 03, 04, 06 and 09 (--leave-out-study), or a list of such prefixes. The runs
    that back the paper leave out Jimenez_Bruno_2026 (the ablation study GSE325414) and, for T and NK cells in
    03, Guo_Zhang_2018. With None the set and the random draws are those of the ordinary run.
    histology: "LUSC" or "LUAD" keeps the cells of that histology only (06 --histology). The cap is drawn
    first, exactly as without the argument, so the two histologies together are the cells of the ordinary run.
    The detection rule is then applied to the cells kept. With None nothing changes."""
    W = work_dir(cell)
    mA, mB = load_meta(W)
    nA, nB = mA["block"].shape[0], mB["block"].shape[0]
    hist = np.concatenate([np.full(nA, "LUSC"), np.full(nB, "LUAD")])
    block = np.concatenate([mA["block"].astype(str), mB["block"].astype(str)])
    pat = np.concatenate([mA["patient"].astype(str), mB["patient"].astype(str)])
    arm = np.concatenate([mA["arm"].astype(str), mB["arm"].astype(str)])
    depth = np.concatenate([mA["ngenes_obs"], mB["ngenes_obs"]]).astype(float)
    genes = np.array([str(g) for g in np.load(W / "genes.npy", allow_pickle=True)])
    present = pd.read_pickle(W / "present.pkl")
    strat_all = np.char.add(np.char.add(block, "|"), hist)
    base = np.arange(len(block)) if leave_out is None else np.where(~left_out(block, leave_out))[0]
    rng = np.random.default_rng(SEED)
    sel = []
    for st in np.unique(strat_all[base]):
        ix = base[strat_all[base] == st]
        sel.append(ix if len(ix) <= cap else rng.choice(ix, cap, replace=False))
    sel = np.sort(np.concatenate(sel))
    if histology is not None:                   # the cells of this histology within the sample of the ordinary run
        sel = sel[hist[sel] == histology]
    n_before = nA + nB
    X = spvstack([load_rows(W, "LUSC", sel[sel < nA]), load_rows(W, "LUAD", sel[sel >= nA] - nA)]).tocsc()
    block, pat, arm = block[sel], pat[sel], arm[sel]
    depth, hist, strat = depth[sel], hist[sel], strat_all[sel]
    N = X.shape[0]
    tot = np.asarray(X.sum(1)).ravel().astype(np.float32)
    GI = {g: i for i, g in enumerate(genes)}
    TX = np.asarray(X[:, GI[SPLITTER]].todense()).ravel() > 0
    arms_u = np.array(sorted(set(arm)))
    acode = np.searchsorted(arms_u, arm)
    PRES = np.zeros((len(arms_u), len(genes)), bool)
    for i, a in enumerate(arms_u):
        if a in present.columns:
            PRES[i] = present[a].values
    nmeas = np.bincount(acode, minlength=len(arms_u)).astype(float) @ PRES
    colnnz = np.diff(X.indptr)
    det = np.where(nmeas > 0, colnnz / np.maximum(nmeas, 1), 0.0)
    keepg = np.where((nmeas >= PRESENCE * N) & (det >= DETFLOOR))[0]
    blks = np.unique(block)
    blk_of = np.searchsorted(blks, block)
    blk_arm = np.array([np.unique(acode[blk_of == b])[0] for b in range(len(blks))])
    pidx = [np.where(pat == p)[0] for p in np.unique(pat)]
    studies = np.unique(block)
    study_of = np.searchsorted(studies, block)
    return dict(cell=cell, W=W, n_before=n_before, sel=sel, nA=nA, X=X, N=N, block=block, pat=pat, arm=arm, depth=depth,
                hist=hist, strat=strat, tot=tot, genes=genes, GI=GI, TX=TX, PRES=PRES, acode=acode,
                nmeas=nmeas, det=det, keepg=keepg, UNI=set(genes[keepg]), blks=blks, blk_of=blk_of,
                blk_arm=blk_arm, pidx=pidx, studies=studies, study_of=study_of)


def perm_expo(e, pidx, r):
    """TXNRD1 detection permuted among the cells within each patient (pidx lists the cells of each patient)."""
    o = e.copy()
    for m in pidx:
        o[m] = r.permutation(e[m])
    return o


def control_a(S):
    """Depth control: within each stratum, the residual of log10 total counts on log10 genes detected;
    the cells with the highest residuals, as many as there are TXNRD1-detected cells, are marked."""
    tot, depth, strat, TX = S["tot"], S["depth"], S["strat"], S["TX"]
    N = len(tot)
    lg_t, lg_g = np.log10(np.maximum(tot, 1)), np.log10(np.maximum(depth, 1))
    resid = np.zeros(N)
    for st in np.unique(strat):
        m = strat == st
        if m.sum() < 20 or lg_g[m].std() == 0:
            continue
        b = np.polyfit(lg_g[m], lg_t[m], 1)
        resid[m] = lg_t[m] - np.polyval(b, lg_g[m])
    ctrl = np.zeros(N, bool)
    for st in np.unique(strat):
        m = np.where(strat == st)[0]
        k = int(round(TX[m].mean() * len(m)))
        if k > 0:
            ctrl[m[np.argsort(-resid[m])[:k]]] = True
    return ctrl


def read_gmt(path):
    """set id -> (description, genes in file order, duplicates dropped)."""
    sets = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 3:
            continue
        sets[f[0]] = (f[1], [g for g in dict.fromkeys(f[2:]) if g])
    return sets


# ---------------------------------------------------------------- batch and study names
PLATFORM = r"_(?:10x|Smart-seq2|BD Rhapsody|inDrop|GEXSCOPE)"


def study_of_batch(batch):
    """First author and year of a batch name ("Leader_Merad_2021_10x 3' v2" -> "Leader 2021");
    the batch name without its platform part when it has no year."""
    base = re.split(PLATFORM, batch)[0]
    p = base.split("_")
    yr = [x for x in p if re.fullmatch(r"(19|20)\d\d", x)]
    return f"{p[0]} {yr[0]}" if yr else base


def batch_labels(batches):
    """Printed column labels of Supplementary Figure 4 (study and platform; repeated labels numbered)."""
    def plat(full):
        m = re.search(r"_(10x[^_]*|Smart-seq2|BD Rhapsody|inDrop|GEXSCOPE)", full)
        return (m.group(1) if m else "").replace("10x 3' ", "10x ")
    a = [study_of_batch(s) for s in batches]
    cnt = Counter(a)
    lab = [f"{x} ({plat(s)})" if cnt[x] > 1 else x for x, s in zip(a, batches)]
    c2 = Counter(lab)
    seen = Counter()
    res = []
    for lb in lab:
        if c2[lb] > 1:
            seen[lb] += 1
            lb = f"{lb} #{seen[lb]}"
        res.append(lb)
    return dict(zip(batches, res))


# US spelling for labels taken from the catalog
US = [("Programmes", "Programs"), ("programmes", "programs"), ("Programme", "Program"), ("programme", "program"),
      ("Colour", "Color"), ("colour", "color"), ("Tumour", "Tumor"), ("tumour", "tumor"),
      ("signalling", "signaling"), ("fibre", "fiber"), ("labelled", "labeled")]


def us(s):
    """Replace UK spellings in a label by US spellings."""
    for a, b in US:
        s = s.replace(a, b)
    return s


# ---------------------------------------------------------------- program scores and per-pairing statistics
def score_programs(S, gmt_path, log=print, minset=None):
    """Program setup: tested members of every set in the GMT file (TXNRD1 removed), the within-batch z of
    every member gene's count fraction (count divided by the cell's total counts), and the program score of
    every cell (mean z of its tested members).

    minset: the least number of tested genes a set needs to be scored (default MINSET, 3). The
    thioredoxin step (06) passes 1 so that one-gene sets are scored, and applies the floor of 3 to
    its multi-gene sets itself."""
    minset = MINSET if minset is None else minset
    X, N, genes, GI, UNI = S["X"], S["N"], S["genes"], S["GI"], S["UNI"]
    tot, PRES, blks, blk_of, blk_arm = S["tot"], S["PRES"], S["blks"], S["blk_of"], S["blk_arm"]
    PROG = {}
    for sid, (_, listed) in read_gmt(gmt_path).items():
        PROG[sid] = dict(listed=listed, tested=[g for g in listed if g in UNI and g != SPLITTER],
                         not_testable=[g for g in listed if g not in UNI and g != SPLITTER],
                         splitter_removed=SPLITTER in listed)
    names = [n for n, p in PROG.items() if len(p["tested"]) >= minset]
    R = len(names)
    mem_genes = sorted({g for n in names for g in PROG[n]["tested"]})
    mpos = {g: i for i, g in enumerate(mem_genes)}
    M_of = [np.array([mpos[g] for g in PROG[n]["tested"]]) for n in names]
    G = len(mem_genes)
    ZF = np.zeros((N, G), np.float32)
    FIN = np.zeros((N, G), bool)
    gidx = np.array([GI[g] for g in mem_genes])
    CH = 1000
    for c0 in range(0, G, CH):
        idx = gidx[c0:c0 + CH]
        A = X[:, idx].toarray().astype(np.float32)
        A /= tot[:, None]
        Z = np.full(A.shape, np.nan, np.float32)
        for b in range(len(blks)):
            mb = blk_of == b
            if mb.sum() < MINBLK:
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
        ZF[:, c0:c0 + len(idx)] = np.where(fin, Z, 0.0)
        FIN[:, c0:c0 + len(idx)] = fin
        del Z, fin
    SUM = np.zeros((N, R), np.float32)
    CNT = np.zeros((N, R), np.float32)
    for r in range(R):
        SUM[:, r] = ZF[:, M_of[r]].sum(1)
        CNT[:, r] = FIN[:, M_of[r]].sum(1)
    SCORE = np.where(CNT > 0, SUM / np.maximum(CNT, 1), np.nan).astype(np.float32)
    off = np.concatenate([[0], np.cumsum([len(m) for m in M_of])]).astype(int)
    return dict(PROG=PROG, names=names, R=R, mem_genes=mem_genes, M_of=M_of, G=G, ZF=ZF, FIN=FIN,
                SUM=SUM, CNT=CNT, SCORE=SCORE, off=off)


class Scorer:
    """The statistics of one pairing (a = exposed cells, b = their partners): base (effect per program,
    the mean paired score difference), logo (leave one gene out, used by 06) and hist (effects per
    histology over all pairs, and their difference)."""

    ARRAYS = ("SCORE", "ZF", "FIN", "SUM", "CNT", "hist_code", "off", "mof_flat")

    def __init__(self, arrays):
        for k in self.ARRAYS:
            setattr(self, k, arrays[k])
        self.R = self.SCORE.shape[1]
        self.K = int(self.off[-1])
        self.M_of = [np.asarray(self.mof_flat[self.off[r]:self.off[r + 1]]) for r in range(self.R)]

    @classmethod
    def build(cls, P, S):
        """A Scorer from the output of score_programs (P) and the working set (S)."""
        arrays = dict(SCORE=P["SCORE"], ZF=P["ZF"], FIN=P["FIN"], SUM=P["SUM"], CNT=P["CNT"],
                      hist_code=(S["hist"] == "LUSC").astype(np.int8), off=P["off"],
                      mof_flat=np.concatenate(P["M_of"]).astype(int))
        return cls(arrays)

    def save(self, folder):
        """Write the arrays to .npy files in a folder."""
        folder = Path(folder)
        folder.mkdir(parents=True, exist_ok=True)
        for k in self.ARRAYS:
            np.save(folder / f"{k}.npy", getattr(self, k))

    @classmethod
    def load(cls, folder):
        """Read the arrays back (memory-mapped)."""
        folder = Path(folder)
        return cls({k: np.load(folder / f"{k}.npy", mmap_mode="r") for k in cls.ARRAYS})

    def base(self, a, b):
        """Effect per program: the mean of (score of a minus score of b) over valid pairs, missing below MINPAIR pairs."""
        Sa, Sb = self.SCORE[a], self.SCORE[b]
        V = np.isfinite(Sa) & np.isfinite(Sb)
        D = np.where(V, Sa - Sb, 0.0).astype(np.float32)
        c = V.sum(0)
        e = np.where(c >= MINPAIR, D.sum(0) / np.maximum(c, 1), np.nan)
        return e, D, V, c

    def logo(self, a, b):
        """Effect per program and member gene with that gene left out of the program score."""
        res = np.full(self.K, np.nan, np.float32)
        ZA, ZB, FA, FB = self.ZF[a], self.ZF[b], self.FIN[a], self.FIN[b]
        for r in range(self.R):
            m = self.M_of[r]
            da = self.CNT[a, r][:, None] - FA[:, m]
            db = self.CNT[b, r][:, None] - FB[:, m]
            v = (da > 0) & (db > 0)
            d = ((self.SUM[a, r][:, None] - ZA[:, m]) / np.maximum(da, 1)
                 - (self.SUM[b, r][:, None] - ZB[:, m]) / np.maximum(db, 1))
            c = v.sum(0)
            res[self.off[r]:self.off[r + 1]] = np.where(c >= MINPAIR, np.where(v, d, 0.0).sum(0) / np.maximum(c, 1), np.nan)
        return res

    def hist(self, a, D, V):
        """The effects per histology (and over both) and their LUSC minus LUAD difference, over all pairs
        (the keys carry an "all_" prefix)."""
        res = {}
        ha = self.hist_code[a]
        for h in ("ALL", "LUSC", "LUAD"):
            m = np.ones(len(a), bool) if h == "ALL" else ha == (1 if h == "LUSC" else 0)
            c = V[m].sum(0)
            res[f"all_{h}"] = np.where(c >= MINPAIR, D[m].sum(0) / np.maximum(c, 1), np.nan)
            res[f"all_n_{h}"] = c
            res[f"all_npairs_{h}"] = int(m.sum())
        res["all_DIFF"] = res["all_LUSC"] - res["all_LUAD"]
        return res

    def real(self, a, b):
        """Everything the REAL arm needs from one pairing: the effect and the histology arm."""
        e, D, V, c = self.base(a, b)
        return dict(e=e, hist=self.hist(a, D, V))


_WORKER = {}


def worker_init(folder):
    """Pool initializer: load the scorer arrays once in each worker process."""
    _WORKER["scorer"] = Scorer.load(folder)


def worker_real(pair):
    """Pool task: the statistics of one pairing (a, b)."""
    a, b = pair
    return _WORKER["scorer"].real(a, b)
