#!/usr/bin/env python
"""
Atlas stage 02, step merge_add: add arm_emtab13526 (De Zuani / Cvejic 2024) to panel_object.h5ad.

Paper items: part of the atlas behind Figure 2, Figure 3, Supplementary Figures 2 to 4 and Supplementary Data 1.
Produces the 2,043,367-cell, 21-block draft object. No random numbers are drawn.
Run: python atlas/steps/stage2_merge_01_add_emtab13526.py   (normally through atlas/02_arm_emtab13526.py)

Ordering, deliberately:
  1. back up stage1a/panel_object.h5ad  ->  stage1a/panel_object_prestage2.h5ad
  2. build the merged object into stage2/panel_object_draft.h5ad
  3. (merge_blockqc) recompute the block QC dropout test on the DRAFT, all blocks
  4. (merge_verify) verify the draft against the backup, then move it into place

Nothing is written to stage1a/panel_object.h5ad by this script.

Gene axis: the merged object keeps its own 1,906-gene axis in its own order. The arm is
subset to that axis by gene symbol, using only genes with var.present_in_source True, because
an all-zero column is a fabricated non-detection. The arm gene order is the panel_genes.csv
order and the merged gene order is a subsequence of it, so the per-row column indices stay
sorted after remapping; that is asserted, not assumed.

X of the existing 1,146,594 cells is copied byte-for-byte value-for-value, chunked, never
re-derived. Its indptr is reused verbatim.

Every categorical remap is an explicit dict with a completeness assertion. No substring rule.
"""
import datetime as dt
import json
import os
import shutil
import sys
import time

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem, write_elem

T0 = time.time()
import os as _os  # noqa: E402
import sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # atlas/
from common import MD  # noqa: E402  paths from config.json, see atlas/common.py
BASE = MD + "/"
S1A = BASE + "stage1a/"
S2 = BASE + "stage2/"
CUR = S1A + "panel_object.h5ad"
BAK = S1A + "panel_object_prestage2.h5ad"
ARM = S2 + "arm_emtab13526.h5ad"
DRAFT = S2 + "panel_object_draft.h5ad"
SRC_KEY = "EMTAB13526"

CELL_CHUNK = 100_000
NNZ_CHUNK = 20_000_000


def log(*a):
    """Print a message with the elapsed time."""
    print(f"[{time.time()-T0:7.1f}s]", *a, flush=True)


def die(msg):
    """Print a stop message and exit with status 1."""
    print("\n*** STOP: " + msg, flush=True)
    sys.exit(1)


# ------------------------------------------------------------------ 1. backup
if not os.path.exists(CUR):
    die("current panel_object.h5ad not found")
if os.path.exists(BAK):
    if os.path.getsize(BAK) != os.path.getsize(CUR):
        die(f"{BAK} already exists and differs in size from the current object. Refusing.")
    log(f"backup already present and same size, keeping it: {BAK}")
else:
    log(f"backing up {CUR} -> {BAK}")
    shutil.copy2(CUR, BAK)
    log(f"backup written, {os.path.getsize(BAK):,} bytes")

# ------------------------------------------------------------------ 2. gene axis
fc = h5py.File(CUR, "r")
fa = h5py.File(ARM, "r")

pg = np.array([g.decode() for g in fc["var/gene_symbol"][:]])
var_cur = read_elem(fc["var"])
var_arm = read_elem(fa["var"])
arm_genes = np.array(var_arm.index.astype(str))
arm_present = var_arm["present_in_source"].astype(bool).values

n_cur = int(fc["X"].attrs["shape"][0])
n_arm = int(fa["X"].attrs["shape"][0])
n_genes = len(pg)
log(f"current object {n_cur:,} x {n_genes}   arm {n_arm:,} x {len(arm_genes)}")

if pd.Index(arm_genes).duplicated().any():
    die("duplicate gene symbols in the arm var")
if pd.Index(pg).duplicated().any():
    die("duplicate gene symbols in the current object var")

arm_pos = {g: i for i, g in enumerate(arm_genes)}
arm_present_set = set(arm_genes[arm_present])

keep_genes = [g for g in pg if g in arm_present_set]
lost_genes = [g for g in pg if g not in arm_present_set]
log(f"INTERSECTION = {len(keep_genes)} genes; lost {len(lost_genes)}: {lost_genes}")

panel = pd.read_csv(S1A + "panel_genes.csv")
panel_ss = panel.set_index("gene_symbol")["source_set"].astype(str)
if lost_genes:
    for g in lost_genes:
        ss = panel_ss.get(g, "")
        curated = ss.replace("hvg_luca_2000", "").strip(",").strip()
        log(f"    LOST {g}  source_set={ss}  curated={'YES' if curated else 'no (HVG only)'}")

# curated redox sets, named explicitly so a loss cannot pass unremarked
REDOX_TOKENS = ["redox", "nrf2", "thioredoxin", "glutathione", "peroxiredoxin", "oxidative"]
redox_genes = set(panel.loc[
    panel["source_set"].astype(str).str.lower().str.contains("|".join(REDOX_TOKENS)),
    "gene_symbol"].astype(str))
redox_lost = sorted(set(lost_genes) & redox_genes)
log(f"curated redox genes in the panel: {len(redox_genes)}; lost by this merge: {len(redox_lost)} {redox_lost}")

if len(keep_genes) != n_genes:
    log("!!! the gene axis SHRINKS. Every downstream number changes. Reporting and continuing "
        "is not appropriate without a decision, so this script stops.")
    die(f"gene axis would shrink from {n_genes} to {len(keep_genes)}")

# column lookup: arm column -> merged column, -1 to drop
lut = np.full(len(arm_genes), -1, dtype=np.int32)
for j, g in enumerate(pg):
    lut[arm_pos[g]] = j
kept_arm_cols = np.where(lut >= 0)[0]
if not np.all(np.diff(lut[kept_arm_cols]) > 0):
    die("the merged gene order is NOT a monotonic subsequence of the arm gene order; "
        "per-row indices would need re-sorting")
log(f"column map built: {len(kept_arm_cols)} arm columns kept, "
    f"{len(arm_genes)-len(kept_arm_cols)} dropped, order monotonic PASS")

# ------------------------------------------------------------------ 3. obs
log("reading obs of both")
obs_cur = read_elem(fc["obs"])
obs_arm = read_elem(fa["obs"])
CUR_COLS = list(obs_cur.columns)
log(f"current obs columns ({len(CUR_COLS)}): {CUR_COLS}")

LEVEL2_MAP = {
    "AT1": "AT1", "AT2": "AT2", "B_cell": "B_cell", "CD4_T": "CD4_T", "CD8_T": "CD8_T",
    "DC": "DC", "NK": "NK", "SMC": "SMC", "Treg": "Treg", "ciliated": "ciliated",
    "club": "club", "endothelial": "endothelial",
    "endothelial_lymphatic": "endothelial_lymphatic", "epithelial_other": "epithelial_other",
    "fibroblast": "fibroblast", "macrophage": "macrophage", "malignant": "malignant",
    "mast": "mast", "mesothelial": "mesothelial", "monocyte": "monocyte",
    "myeloid_other": "myeloid_other", "neutrophil": "neutrophil", "pDC": "pDC",
    "pericyte": "pericyte", "plasma_cell": "plasma_cell",
    "B": "B_cell", "plasma": "plasma_cell", "lymphatic_EC": "endothelial_lymphatic",
    "epithelial": "epithelial_other", "basal": "basal", "erythroid": "erythroid",
    # NEW with this arm: T cells the authors did not resolve to CD4/CD8/Treg/NK
    "T_other": "T_other",
    "": "",
}
COMPARTMENTS = ["epithelial_malignant", "T_NK", "myeloid", "B_plasma", "fibroblast",
                "mural", "endothelial", "other"]

seen = set(obs_arm["level2_celltype"].astype(str).unique())
missing = seen - set(LEVEL2_MAP)
if missing:
    die(f"arm level2_celltype values with no explicit map entry: {sorted(missing)}")
bad = set(obs_arm["level1_compartment"].astype(str).unique()) - set(COMPARTMENTS) - {""}
if bad:
    die(f"arm level1_compartment outside the allowed vocabulary: {sorted(bad)}")

o = pd.DataFrame(index=obs_arm.index.astype(str))
for c in CUR_COLS:
    if c in ("block_id_native", "level2_celltype_native", "n_counts_intersect",
             "n_genes_intersect", "block_qc_pass"):
        continue
    if c not in obs_arm.columns:
        die(f"arm obs is missing the required column {c}")
    o[c] = obs_arm[c].values

o["level2_celltype"] = obs_arm["level2_celltype"].astype(str).map(LEVEL2_MAP).values
o["block_id_native"] = obs_arm["block_id"].astype(str).values
o["level2_celltype_native"] = obs_arm["level2_celltype_map_native"].astype(str).values
o["block_qc_pass"] = False          # provisional; merge_blockqc recomputes it for every block
o["n_counts_intersect"] = np.float32(0)   # filled during the X pass
o["n_genes_intersect"] = np.int32(0)

# block_id must already be the canonical <study>_<platform>; this study has one source_key
canon = obs_arm["study"].astype(str) + "_" + obs_arm["platform"].astype(str)
if not (canon.values == obs_arm["block_id"].astype(str).values).all():
    die("arm block_id does not equal the canonical <study>_<platform>")
log("arm block_id == canonical <study>_<platform> PASS")

o = o[CUR_COLS]
log(f"arm obs harmonised to the merged schema {o.shape}")

obs = pd.concat([obs_cur, o], axis=0)
del o
log(f"obs concatenated {obs.shape}")

CAT_COLS = ["source_key", "study", "donor_id", "library_id", "histology", "histology_source",
            "histology_confidence", "tissue_class", "sampling_site", "platform", "modality",
            "enrichment", "treatment", "block_id", "level1_compartment", "level2_celltype",
            "celltype_source", "donor_id_native", "block_id_native", "level2_celltype_native",
            "patient_id", "technical_replicate_group", "tnm", "sex", "smoking"]
for c in CAT_COLS:
    obs[c] = obs[c].astype(str).astype("category")
for c in ["level3_celltype_native", "stage"]:
    obs[c] = obs[c].astype(str)
for c in ["composition_safe", "counts_are_estimates", "block_qc_pass"]:
    obs[c] = obs[c].astype(bool)
for c in ["is_malignant", "macrophage_strict"]:
    obs[c] = pd.array(obs[c].values, dtype="boolean")
obs["n_counts"] = obs["n_counts"].astype("int64")
obs["n_genes"] = obs["n_genes"].astype("int32")
for c in ["celltype_confidence", "malignancy_score", "doublet_score", "age"]:
    obs[c] = obs[c].astype("float64")
obs["n_counts_intersect"] = obs["n_counts_intersect"].astype("float32")
obs["n_genes_intersect"] = obs["n_genes_intersect"].astype("int32")

n_tot = len(obs)

# --- assertions on obs that do not need X ------------------------------------------
log("OBS ASSERTIONS")
ids = obs.index
dup = ids.duplicated(keep=False)
if dup.any():
    d = pd.Series(ids[dup]).value_counts()
    die(f"cell_id collisions: {int(dup.sum())} rows, e.g. {d.head(10).to_dict()}")
log(f"  cell_id unique: {n_tot:,} / {n_tot:,} PASS")
pipes = pd.Series(ids).str.count(r"\|")
log(f"  cell_id has exactly two pipes on {int(pipes.eq(2).sum()):,} of {n_tot:,} "
    f"({'PASS' if pipes.eq(2).all() else 'FAIL'})")
pref = pd.Series(ids).str.split("|").str[0].value_counts().to_dict()
log(f"  cell_id source prefixes: {pref}")

chk = obs.groupby("block_id", observed=True)["source_key"].nunique()
log(f"  no block spans two source_keys: {'PASS' if (chk<=1).all() else 'FAIL ' + str(chk[chk>1].to_dict())}")
chk = obs.groupby("donor_id", observed=True)["histology"].nunique()
log(f"  no donor spans two histologies: {'PASS' if (chk<=1).all() else 'FAIL ' + str(chk[chk>1].to_dict())}")
cs = obs["composition_safe"].values
wt = (obs["enrichment"].astype(str) == "whole_tissue").values
log(f"  composition_safe == (enrichment == whole_tissue): {'PASS' if (cs==wt).all() else 'FAIL'} "
    f"({int(cs.sum()):,} safe cells)")
new = (obs["source_key"].astype(str) == SRC_KEY).values
log(f"  composition_safe False on every {SRC_KEY} cell: "
    f"{'PASS' if not cs[new].any() else 'FAIL'} ({int(new.sum()):,} cells)")
mast = obs["level2_celltype"].astype(str) == "mast"
log(f"  mast cells all in compartment other: "
    f"{'PASS' if (obs.loc[mast,'level1_compartment'].astype(str)=='other').all() else 'FAIL'} "
    f"({int(mast.sum()):,} cells)")
sub_t = obs["level3_celltype_native"].str.contains("t cell", case=False, regex=False)
sub_t_not_tnk = sub_t & (obs["level1_compartment"].astype(str) != "T_NK")
log(f"  substring trap live check: {int(sub_t.sum()):,} native labels contain the substring "
    f"'t cell', of which {int(sub_t_not_tnk.sum()):,} are correctly NOT T_NK "
    f"({obs.loc[sub_t_not_tnk,'level1_compartment'].astype(str).value_counts().to_dict()})")
trail = (obs["level3_celltype_native"] == "Cycling anti-inflammatory alveolar macrophages ")
log(f"  trailing-space label preserved: {int(trail.sum()):,} cells")
log(f"  n_counts > 0 on every cell: {'PASS' if (obs.n_counts.values>0).all() else 'FAIL'}")

viol = []
REQUIRED = [
    "source_key", "study", "donor_id", "library_id", "histology", "histology_source",
    "histology_confidence", "tissue_class", "platform", "modality", "enrichment",
    "composition_safe", "treatment", "block_id", "level1_compartment", "level2_celltype",
    "level3_celltype_native", "celltype_source", "n_counts", "n_genes", "counts_are_estimates",
]
for c in REQUIRED:
    s = obs[c]
    n_na = int(pd.isna(s).sum())
    n_empty = int((s.astype(str).str.strip() == "").sum()) if s.dtype.name in ("object", "category") else 0
    if n_na or n_empty:
        by = obs.loc[s.astype(str).str.strip() == "", "source_key"].astype(str).value_counts().to_dict()
        viol.append(f"{c}: {n_na+n_empty:,} cells empty/NA; by source_key {by}")
if viol:
    for v in viol:
        log(f"  NOTE required column {v}")
else:
    log("  all required columns present and non-empty on every cell PASS")

# ------------------------------------------------------------------ 4. var and uns
var = var_cur.copy()
var["ensembl_id_emtab13526"] = var_arm.loc[pg, "ensembl_id_in_source"].astype(str).values
var["gene_id_resolution_emtab13526"] = var_arm.loc[pg, "gene_id_resolution"].astype(str).values
var["source_symbol_in_reference_emtab13526"] = var_arm.loc[pg, "source_symbol_in_reference"].astype(str).values
var["symbol_agrees_emtab13526"] = var_arm.loc[pg, "symbol_agrees_with_reference"].astype(bool).values
agree_all = (var["ensembl_id_luca"].astype(str) == var["ensembl_id_emtab13526"].astype(str))
log(f"Ensembl id agrees between LuCA and EMTAB13526 on {int(agree_all.sum())} of {len(var)} genes; "
    f"symbol differs in the 2018 reference for {int((~var['symbol_agrees_emtab13526']).sum())} genes")

uns = read_elem(fc["uns"])
uns_arm = read_elem(fa["uns"])
uns["stage"] = "2"
uns["built_utc"] = dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")
uns["script"] = ("atlas/steps/stage2_merge_01_add_emtab13526.py, "
                 "stage2_merge_02_blockqc.py, stage2_merge_03_verify.py")
uns["arms"] = list(uns.get("arms", [])) + ["arm_emtab13526.h5ad"]
uns["level2_map"] = json.dumps(LEVEL2_MAP)
uns["level2_new_levels_stage2"] = (
    "T_other is new with EMTAB13526: 35,716 cells, from the author labels 'Downregulated T "
    "cells' (30,177) and 'Gamma delta T cells' (5,539). No earlier arm uses it.")
uns["emtab13526_level2_caveat"] = (
    "The arm maps the author label 'Naive T cells' (52,665 cells) to CD4_T. The author label "
    "does not resolve CD4 from CD8, so those cells are naive T cells of unknown lineage counted "
    "as CD4_T. uns['emtab13526']['level2_new_levels'] as written by the arm build says they went "
    "to T_other; the object says CD4_T. The object is what is true. Any CD4-specific analysis "
    "must decide what to do with them; obs.level3_celltype_native identifies them exactly.")
uns["emtab13526"] = {k: uns_arm[k] for k in [
    "accession", "study_full", "X_is", "gene_space", "gene_id_convention", "histology_note",
    "tissue_class_note", "enrichment_note", "platform_note", "depth_note", "celltype_note",
    "celltype_low_confidence", "macrophage_strict_note", "immune_contrast",
    "fibroblast_limitation", "structural_compartment_limitation", "sorted_fraction_limitation",
    "level2_new_levels"] if k in uns_arm}
uns["emtab13526"]["obs_columns_not_carried"] = (
    "The arm file keeps 20 native columns this object does not carry, joined on cell_id: "
    "barcode, author_object, author_barcode, author_sample, author_sorting, author_environment, "
    "sangerID, level2_celltype_map_native, celltype_confidence_word, author_flagged_remove, "
    "histology_disputed, tissue_class_disputed, enrichment_disputed, platform_disputed, "
    "n_counts_panel, n_genes_panel, author_total_counts, author_n_genes, pct_counts_mt. "
    "They live in output/RNA/meta_dataset/stage2/arm_emtab13526.h5ad.")
uns["X_is"] = (str(uns["X_is"]) + " EMTAB13526 is integer UMI counts stored as float32.")
uns["panel_genes_lost_stage2"] = [str(g) for g in lost_genes]

# ------------------------------------------------------------------ 5. write draft
log(f"writing {DRAFT}")
if os.path.exists(DRAFT):
    os.remove(DRAFT)

with h5py.File(DRAFT, "w") as fo:
    fo.attrs["encoding-type"] = "anndata"
    fo.attrs["encoding-version"] = "0.1.0"
    write_elem(fo, "var", var)
    for k in ["layers", "obsm", "obsp", "varm", "varp"]:
        write_elem(fo, k, {})

    g = fo.create_group("X")
    g.attrs["encoding-type"] = "csr_matrix"
    g.attrs["encoding-version"] = "0.1.0"
    g.attrs["shape"] = np.array([n_tot, n_genes], dtype="int64")

    nnz_cur = fc["X/data"].shape[0]
    d_out = g.create_dataset("data", shape=(nnz_cur,), maxshape=(None,), dtype="float32",
                             chunks=(40049,), compression="gzip", compression_opts=4)
    i_out = g.create_dataset("indices", shape=(nnz_cur,), maxshape=(None,), dtype="int32",
                             chunks=(40049,), compression="gzip", compression_opts=4)

    log("  copying the existing X, chunked")
    for s in range(0, nnz_cur, NNZ_CHUNK):
        e = min(s + NNZ_CHUNK, nnz_cur)
        d_out[s:e] = fc["X/data"][s:e]
        i_out[s:e] = fc["X/indices"][s:e]
        log(f"    {e:,} / {nnz_cur:,}")
    indptr_cur = fc["X/indptr"][:].astype(np.int64)
    assert int(indptr_cur[-1]) == nnz_cur

    log("  appending the arm, remapped to the merged gene axis")
    arm_indptr = fa["X/indptr"][:].astype(np.int64)
    new_ptr = np.zeros(n_arm + 1, dtype=np.int64)
    nc = np.zeros(n_arm, dtype=np.float64)
    ng = np.zeros(n_arm, dtype=np.int64)
    write_at = nnz_cur
    total_arm_nnz_in = 0
    for start in range(0, n_arm, CELL_CHUNK):
        stop = min(start + CELL_CHUNK, n_arm)
        lo, hi = int(arm_indptr[start]), int(arm_indptr[stop])
        d = fa["X/data"][lo:hi].astype(np.float32)
        idx = fa["X/indices"][lo:hi]
        total_arm_nnz_in += hi - lo
        rows = np.repeat(np.arange(start, stop, dtype=np.int64),
                         np.diff(arm_indptr[start:stop + 1]))
        newcol = lut[idx]
        keep = newcol >= 0
        d_k, c_k, r_k = d[keep], newcol[keep].astype(np.int32), rows[keep]
        cnt = np.bincount(r_k - start, minlength=stop - start)
        new_ptr[start + 1:stop + 1] = cnt
        nc[start:stop] = np.bincount(r_k - start, weights=d_k.astype(np.float64),
                                     minlength=stop - start)
        ng[start:stop] = cnt
        d_out.resize((write_at + d_k.size,))
        i_out.resize((write_at + d_k.size,))
        d_out[write_at:write_at + d_k.size] = d_k
        i_out[write_at:write_at + d_k.size] = c_k
        write_at += d_k.size
        if (start // CELL_CHUNK) % 2 == 0:
            log(f"    cells {stop:,} / {n_arm:,}   nnz written {write_at:,}")
    assert total_arm_nnz_in == fa["X/data"].shape[0], "arm nnz lost while reading"
    new_ptr = np.cumsum(new_ptr) + nnz_cur
    assert int(new_ptr[-1]) == write_at, "arm indptr does not close"
    indptr = np.concatenate([indptr_cur, new_ptr[1:]])
    assert indptr.shape[0] == n_tot + 1
    assert indptr[-1] == write_at
    if write_at > np.iinfo(np.int32).max:
        die("nnz exceeds int32; indptr dtype must change")
    ip = g.create_dataset("indptr", shape=(n_tot + 1,), dtype="int32", chunks=(8958,),
                          compression="gzip", compression_opts=4)
    ip[:] = indptr.astype(np.int32)
    log(f"  X written: {n_tot:,} x {n_genes}, nnz {write_at:,} "
        f"(existing {nnz_cur:,} + arm {write_at-nnz_cur:,})")

    obs.loc[obs.index[n_cur:], "n_counts_intersect"] = nc.astype(np.float32)
    obs.loc[obs.index[n_cur:], "n_genes_intersect"] = ng.astype(np.int32)
    obs["n_counts_intersect"] = obs["n_counts_intersect"].astype("float32")
    obs["n_genes_intersect"] = obs["n_genes_intersect"].astype("int32")
    write_elem(fo, "obs", obs)
    write_elem(fo, "uns", uns)

fc.close()
fa.close()
log(f"draft written {DRAFT} ({os.path.getsize(DRAFT):,} bytes)")

# sanity on the new depth columns
sub = obs.iloc[n_cur:]
log(f"  new cells: n_counts_intersect / n_counts median ratio "
    f"{float(np.median(sub.n_counts_intersect.values / sub.n_counts.values)):.4f}; "
    f"cells with zero panel counts {int((sub.n_counts_intersect.values==0).sum()):,}")
log("step 1 done")
