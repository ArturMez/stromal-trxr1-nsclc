#!/usr/bin/env python3
"""
Atlas stage 04, step verify: verify the merged object before it is analyzed.

Paper items: none directly; it checks fullgene_object.h5ad, the atlas behind Figure 2, Figure 3,
Supplementary Figures 2 to 4 and Supplementary Data 1. Random numbers: a numpy generator with seed 31
(default_rng(31)) draws the sampled cells of checks 2, 3, 4 and 6; N_ROUNDTRIP = 300 cells per arm,
N_MASK = 4,000 cells per arm for the mask check, 600 cells for the estimates flag.

A 21 GB object should not be trusted because the script that wrote it exited zero.
Six checks, each of which can fail independently, and each printed as a pass or a fail rather than
summarized. A single FAIL means the object is not fit to analyze.

  1  STRUCTURE. indptr is monotone, its last entry equals the nonzero count, its length equals cells
     plus one, and the declared shape matches obs and var.
  2  COUNTS ROUND-TRIP. Random cells are read from the merged object and compared gene by gene to the
     arm they came from, on the union gene names rather than on positions, so a remapping error cannot
     hide. This is the check that a gene-index remap silently shifting by one would fail.
  3  THE PRESENCE MASK IS TRUE. For each arm, genes the mask marks ABSENT must carry zero counts in
     every cell of that arm. If a gene marked absent has counts, the mask is wrong and every gene-set
     score that trusts it is wrong. Sampled over cells, exhaustive over the absent genes.
  4  ROW ORDER IS ASCENDING. Within each sampled row the gene indices must be strictly increasing,
     which CSR consumers assume and which the per-arm remap could have broken.
  5  ANNOTATIONS SURVIVED. Cell ids unique; `block_qc_pass`, `counts_are_estimates`, `histology`,
     `tissue_class`, `donor_id` all present and non-empty; the compartment counts and the patient set
     match `panel_object.h5ad` exactly, because losing a patient silently would change the cohort.
  6  THE ESTIMATES FLAG IS RIGHT. Cells flagged `counts_are_estimates` should be the non-integer ones
     and unflagged cells should be integral. Checked on a sample, reported as a rate.

RUN: python3 atlas/steps/fullbuild_04_verify.py   (normally through atlas/04_build_objects.py)
FROM: fullgene_object.h5ad, the seven arm_*_full.h5ad files and stage1a/panel_object.h5ad.
OUT -> output/RNA/meta_dataset/fullbuild/verify_report.csv
"""
import os
import numpy as np
import pandas as pd
import h5py

import os as _os  # noqa: E402
import sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # atlas/
from common import ROOT, SCRNA, SCRATCH, DATA, RSCRIPT  # noqa: E402,F401  paths from config.json
D = f"{ROOT}/output/RNA/meta_dataset/fullbuild"
OBJ = os.environ.get("FULLGENE_DST", f"{D}/fullgene_object.h5ad")
PANEL = f"{ROOT}/output/RNA/meta_dataset/stage1a/panel_object.h5ad"
ARMS = {"LuCA": f"{D}/arm_luca_full.h5ad",
        "E-MTAB-13526": f"{D}/arm_emtab13526_full.h5ad",
        "E-MTAB-6653": f"{D}/arm_emtab6653_full.h5ad",
        "E-MTAB-6149": f"{D}/arm_emtab6149_full.h5ad",
        "STTT": f"{D}/arm_sttt_full.h5ad",
        "GSE325414": f"{D}/arm_gse325414_full.h5ad",
        "GSE246215": f"{D}/arm_gse246215_full.h5ad"}
N_ROUNDTRIP, N_MASK = 300, 4000
rng = np.random.default_rng(31)   # seed 31: the sampled cells of the checks
res = []


def dec(a):
    """Decode an array of bytes or strings to an array of str."""
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in a])


def col(h, path):
    """Read one obs column of an h5ad file (categorical, nullable or plain) as an array."""
    g = h[path]
    if isinstance(g, h5py.Group):
        if "categories" in g:
            return dec(g["categories"][:])[g["codes"][:]]
        if "values" in g:
            return g["values"][:]
    v = g[:]
    return dec(v) if v.dtype.kind in "SO" else v


def record(check, ok, detail):
    """Append one check result (PASS or FAIL with detail) to the report and print it."""
    res.append(dict(check=check, result="PASS" if ok else "FAIL", detail=detail))
    print(f"  [{'PASS' if ok else 'FAIL'}] {check}: {detail}", flush=True)


if not os.path.exists(OBJ):
    raise SystemExit("fullgene_object.h5ad not found")

print("loading axes ...", flush=True)
with h5py.File(OBJ, "r") as h:
    genes = dec(h["var/gene_symbol"][:])
    ids = dec(h["obs/_index"][:])
    ip = h["X/indptr"][:]
    shape = tuple(int(x) for x in h["X"].attrs["shape"])
    nnz = int(h["X/data"].shape[0])
    obs = {k: col(h, f"obs/{k}") for k in ["level1_compartment", "source_key", "block_qc_pass",
                                           "counts_are_estimates", "histology", "tissue_class",
                                           "donor_id", "patient_id", "n_genes",
                                           "block_qc_assessed"] if k in h["obs"]}
    armnames = [k for k in h["varm/gene_present"].keys() if k != "gene_symbol"]
    mask = {a: h[f"varm/gene_present/{a}"][:].astype(bool) for a in armnames}
gpos = {g: i for i, g in enumerate(genes)}

# ---- 1 structure
ok = (len(ip) == shape[0] + 1 and int(ip[-1]) == nnz and np.all(np.diff(ip) >= 0)
      and len(ids) == shape[0] and len(genes) == shape[1])
record("1 structure", ok, f"{shape[0]:,} cells x {shape[1]:,} genes, {nnz:,} nonzeros, "
                          f"indptr[-1]={int(ip[-1]):,}")

# ---- 2 counts round-trip + 4 ascending order
src = pd.Series(obs["source_key"]).astype(str).values
bad_rt, checked_rt, bad_order = 0, 0, 0
pos_in_obj = pd.Series(np.arange(len(ids)), index=ids)
with h5py.File(OBJ, "r") as ho:
    for aname, apath in ARMS.items():
        if not os.path.exists(apath):
            continue
        with h5py.File(apath, "r") as ha:
            aid = dec(ha["obs/_index"][:])
            asym = dec(ha["var/gene_symbol"][:])
            aip = ha["X/indptr"][:]
            pick = rng.choice(len(aid), min(N_ROUNDTRIP, len(aid)), replace=False)
            for r in np.sort(pick):
                a0, a1 = int(aip[r]), int(aip[r + 1])
                agi = ha["X/indices"][a0:a1]
                adv = ha["X/data"][a0:a1]
                o = int(pos_in_obj.get(aid[r], -1))
                if o < 0:
                    bad_rt += 1
                    continue
                b0, b1 = int(ip[o]), int(ip[o + 1])
                ogi = ho["X/indices"][b0:b1]
                odv = ho["X/data"][b0:b1]
                if len(ogi) > 1 and np.any(np.diff(ogi) <= 0):
                    bad_order += 1
                m_arm = {asym[g]: v for g, v in zip(agi, adv)}
                m_obj = {genes[g]: v for g, v in zip(ogi, odv)}
                if m_arm != m_obj:
                    bad_rt += 1
                checked_rt += 1
record("2 counts round-trip", bad_rt == 0,
       f"{checked_rt:,} cells compared gene-by-gene against their arm, {bad_rt} mismatched")
record("4 row indices ascending", bad_order == 0,
       f"{bad_order} of {checked_rt:,} sampled rows had non-ascending gene indices")

# ---- 3 presence mask is true
KEYMAP = {"LuCA": "LuCA_census", "E-MTAB-13526": "EMTAB13526",
          "E-MTAB-6653": "EMTAB6653", "E-MTAB-6149": "EMTAB6149",
          "STTT": "STTT", "GSE325414": "GSE325414", "GSE246215": "GSE246215"}
viol = 0
detail = []
with h5py.File(OBJ, "r") as ho:
    for aname in armnames:
        sk = KEYMAP.get(aname)
        rows = np.where(src == sk)[0] if sk else np.array([], int)
        if not len(rows):
            continue
        absent = np.where(~mask[aname])[0]
        aset = set(absent.tolist())
        pick = rng.choice(rows, min(N_MASK, len(rows)), replace=False)
        v = 0
        for r in np.sort(pick):
            b0, b1 = int(ip[r]), int(ip[r + 1])
            gidx = ho["X/indices"][b0:b1]
            v += int(len(aset.intersection(gidx.tolist())))
        viol += v
        detail.append(f"{aname}: {len(absent):,} genes marked absent, {v} counts found in "
                      f"{len(pick):,} sampled cells")
record("3 presence mask is true", viol == 0, "; ".join(detail))

# ---- 5 annotations
with h5py.File(PANEL, "r") as h:
    pl1 = col(h, "obs/level1_compartment").astype(str)
    ppat = col(h, "obs/patient_id").astype(str)
T = ["fibroblast", "myeloid", "T_NK"]
pm = np.isin(pl1, T)
# The wave-1 arms ADD cells and patients, so 5b and 5c compare only the panel-derived subset.
PANEL_KEYS = {"LuCA_census", "EMTAB13526", "EMTAB6653", "EMTAB6149"}
from_panel = np.isin(pd.Series(obs["source_key"]).astype(str).values, list(PANEL_KEYS))
want_counts = pd.Series(pl1[pm]).value_counts().to_dict()
got_counts = pd.Series(obs["level1_compartment"]).astype(str)[from_panel].value_counts().to_dict()
want_pat = set(ppat[pm])
got_pat = set(pd.Series(obs["patient_id"]).astype(str)[from_panel])
new_cells = int((~from_panel).sum())
new_pat = set(pd.Series(obs["patient_id"]).astype(str)[~from_panel]) - want_pat
print(f"  (new arms contribute {new_cells:,} cells and {len(new_pat)} patients not in panel_object)",
      flush=True)
lost = want_pat - got_pat
uniq_ids = len(set(ids)) == len(ids)
record("5a cell ids unique", uniq_ids, f"{len(ids):,} ids, {len(set(ids)):,} distinct")
record("5b panel-derived compartment counts match panel_object", want_counts == got_counts,
       f"panel {want_counts} vs merged-subset {got_counts}")
record("5c no patient lost", len(lost) == 0,
       f"panel has {len(want_pat)} patients in these compartments, merged has {len(got_pat)}, "
       f"lost {sorted(lost)[:5] if lost else 'none'}")
# A REQUIRED column that is absent must FAIL, not vanish from the report.
REQUIRED = ["block_qc_pass", "histology", "tissue_class", "donor_id", "patient_id",
            "level1_compartment", "source_key", "n_genes", "counts_are_estimates"]
absent = [k for k in REQUIRED if k not in obs]
record("5d required annotations present", not absent,
       f"missing: {absent}" if absent else f"all {len(REQUIRED)} present")
empt = {k: int((pd.Series(obs[k]).astype(str).isin(["", "nan", "None"])).sum())
        for k in REQUIRED if k in obs}
record("5e required annotations non-empty", all(x == 0 for x in empt.values()), str(empt))

# ---- 6 estimates flag
est = pd.Series(obs["counts_are_estimates"]).astype(str).isin(["True", "1", "true"]).values
agree = dis = 0
with h5py.File(OBJ, "r") as ho:
    pick = rng.choice(len(ids), 600, replace=False)
    for r in np.sort(pick):
        b0, b1 = int(ip[r]), int(ip[r + 1])
        if b1 <= b0:
            continue
        d = ho["X/data"][b0:b1]
        nonint = not np.allclose(d, np.round(d))
        if nonint == bool(est[r]):
            agree += 1
        else:
            dis += 1
record("6 estimates flag matches integrality", dis == 0,
       f"{agree} agree, {dis} disagree over {agree+dis} sampled cells; "
       f"{int(est.sum()):,} cells flagged in total")

R = pd.DataFrame(res)
R.to_csv(f"{D}/verify_report.csv", index=False)
pd.set_option("display.width", 200)
print("\n" + "=" * 104)
print("VERIFICATION OF fullgene_object.h5ad")
print("UNIT: one check. A single FAIL means the object should not be analysed.")
print("=" * 104)
print(R.to_string(index=False))
print(f"\n{int((R.result=='PASS').sum())} of {len(R)} checks passed")
print(f"written -> {D}/verify_report.csv")
