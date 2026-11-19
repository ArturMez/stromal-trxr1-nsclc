"""
Atlas stage 01, steps lt_build, lt_validate, lt_apply, lt_write: label transfer of level1_compartment
onto the 32,341 E-MTAB-6653 cells (a k-nearest-neighbor vote from a labelled LuCA reference).

Paper items: gives the E-MTAB-6653 cells their compartment labels in the atlas behind Figure 2,
Figure 3, Supplementary Figures 2 to 4 and Supplementary Data 1.
Run: python atlas/steps/stage1b_label_transfer_emtab6653_v2.py {build,validate,apply,write}
         [--variant perbatch] [--floor 0.5]
(normally through atlas/01_arms_luca_lambrechts.py, which runs build, validate, apply and write with
--variant perbatch --floor 0.5). Stages: build caches the log-normalized matrix, validate compares the
two scaling variants on held-out data, apply transfers labels onto the query and writes the tables,
write puts the labels into panel_object.h5ad (after saving panel_object_prelabel.h5ad).
Random numbers: numpy default_rng(SEED), SEED = 0, draws the leave-one-study-out test subsample (at most
25,000 cells per held-out study); the PCA and the exact kNN are deterministic.

Design
------
QUERY      EMTAB6653, 32,341 cells, 10x 3' v2, whole_tissue, donors P6/P7/P8.
           level1_compartment empty (the deposit carries no author cell types).

REFERENCE  LuCA census cells, platform 10x 3' v2, enrichment whole_tissue, labelled.
           Chemistry-matched to the query. Sorted libraries are excluded on purpose
           (Goveia stromal_sort, the Leader CD45pos subset): a kNN vote share is a
           local-density statistic, so a compositionally distorted reference would
           distort the confidence of every ambiguous query cell.

VALIDATION (external, headline) EMTAB6149, 10x 3' v2, whole_tissue, donors P3/P4/P5.
           Same laboratory, same paper, same dissociation protocol and same chemistry
           as the query, but annotated independently by the original authors and by
           disjoint donors. This is the closest available analogue of the real task.
VALIDATION (internal) leave-one-study-out across the five reference studies.

METHOD     normalize_total(1e4) -> log1p -> gene scaling -> PCA(50) fitted on the
           reference -> exact kNN (k=25) vote in PC space.
           Two scaling variants are compared on held-out data and the winner is used:
             shared    reference mean/sd applied to reference and query alike
                       (this is what scanpy ingest does)
             perbatch  every batch scaled by its own mean/sd
           Per-cell confidence = fraction of the 25 nearest reference cells carrying
           the winning label. Nothing below the floor is assigned.

Outputs go to output/RNA/meta_dataset/stage1a/ (labeltransfer_validation_summary.csv,
labeltransfer_emtab6653_percell.csv.gz and the by-patient, by-site and epithelial-fraction tables); the
cache is in scratch/labeltransfer_emtab6653/.
"""

import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd
import scipy.sparse as sp

import os as _os  # noqa: E402
import sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # atlas/
from common import ROOT, SCRNA, SCRATCH, DATA, RSCRIPT  # noqa: E402,F401  paths from config.json
OBJ = os.path.join(ROOT, "output/RNA/meta_dataset/stage1a/panel_object.h5ad")
OUT = os.path.join(ROOT, "output/RNA/meta_dataset/stage1a")
CACHE = f"{SCRATCH}/labeltransfer_emtab6653"

N_PCS = 50
K = 25
CLIP = 10.0
SEED = 0

REF_STUDIES = [
    
    "Chen_Zhang_2020",
    "He_Fan_2021",
    "Kim_Lee_2020",
    "Leader_Merad_2021",
]


def log(msg):
    """Print a time-stamped message."""
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def check(cond, msg, detail=""):
    """Print PASS or FAIL for a condition and exit when it fails."""
    tag = "PASS" if cond else "FAIL"
    print(f"[{tag}] {msg}" + (f" :: {detail}" if detail else ""), flush=True)
    if not cond:
        sys.exit(f"ASSERTION FAILED: {msg} :: {detail}")


# ----------------------------------------------------------------------------- load
def build_groups():
    """Select the reference (LuCA, block QC pass), validation (E-MTAB-6149) and query (E-MTAB-6653) cells, log-normalize and cache them."""
    import anndata as ad

    log("reading panel object")
    A = ad.read_h5ad(OBJ)
    o = A.obs
    plat = o.platform.astype(str)
    src = o.source_key.astype(str)
    enr = o.enrichment.astype(str)
    study = o.study.astype(str)
    lab = o.level1_compartment.astype(str)

    qc = o["block_qc_pass"].values.astype(bool) if "block_qc_pass" in o.columns else np.ones(len(o), bool)
    # The reference also requires block_qc_pass: the sorted libraries are excluded because a
    # compositionally distorted reference distorts every ambiguous cell's confidence, and
    # Adams_Kaminski_2020 (19.5 pct of the reference, 70 pct myeloid, CD68 and IGKC at exactly
    # zero counts) fails the block QC and is excluded for the same reason.
    is_ref = (src == "LuCA_census") & (plat == "10x 3' v2") & (enr == "whole_tissue") & (lab != "") & qc
    is_val = (src == "EMTAB6149") & (plat == "10x 3' v2") & (lab != "")
    is_qry = src == "EMTAB6653"

    check(int(is_qry.sum()) == 32341, "query is the 32,341 EMTAB6653 cells", str(int(is_qry.sum())))
    check(int(is_val.sum()) == 48384, "external validation is 48,384 EMTAB6149 v2 cells", str(int(is_val.sum())))
    check(True, "query labels will be overwritten by this v2 run")
    check(
        sorted(study[is_ref].unique()) == sorted(REF_STUDIES),
        "reference is the LuCA v2 whole_tissue studies passing block QC",
        str(sorted(study[is_ref].unique())),
    )
    qd = set(o.donor_id.astype(str)[is_qry])
    vd = set(o.donor_id.astype(str)[is_val])
    check(len(qd & vd) == 0, "validation donors are disjoint from query donors", f"{sorted(qd)} vs {sorted(vd)}")

    idx = np.concatenate([np.where(is_ref)[0], np.where(is_val)[0], np.where(is_qry)[0]])
    role = np.array(["ref"] * int(is_ref.sum()) + ["val"] * int(is_val.sum()) + ["qry"] * int(is_qry.sum()))
    # batch = the unit that shares a library prep / annotation provenance
    batch = np.where(role == "ref", study.values[idx], src.values[idx])

    log(f"ref={int(is_ref.sum()):,}  val={int(is_val.sum()):,}  qry={int(is_qry.sum()):,}")

    X = A.X[idx].astype(np.float32)
    check(sp.isspmatrix_csr(X), "X subset is CSR")
    check(float(X.min()) >= 0, "no negative values in X", str(float(X.min())))

    meta = pd.DataFrame(
        {
            "cell_id": o.index.values[idx],
            "role": role,
            "batch": batch,
            "study": study.values[idx],
            "label": lab.values[idx],
            "donor_id": o.donor_id.astype(str).values[idx],
            "sampling_site": o.sampling_site.astype(str).values[idx],
            "histology": o.histology.astype(str).values[idx],
            "n_counts": o.n_counts.values[idx],
            "n_genes": o.n_genes.values[idx],
        }
    )
    genes = A.var.index.values.copy()
    del A

    # normalize_total(1e4) + log1p, done manually to stay sparse
    log("normalize_total(1e4) + log1p")
    tot = np.asarray(X.sum(axis=1)).ravel()
    check((tot > 0).all(), "every cell has positive total signal", str(float(tot.min())))
    scale = (1e4 / tot).astype(np.float32)
    X = sp.diags(scale) @ X
    X = X.tocsr()
    X.data = np.log1p(X.data)

    os.makedirs(CACHE, exist_ok=True)
    sp.save_npz(os.path.join(CACHE, "lognorm.npz"), X)
    meta.to_parquet(os.path.join(CACHE, "meta.parquet"))
    np.save(os.path.join(CACHE, "genes.npy"), genes)
    log(f"cached lognorm matrix {X.shape} nnz={X.nnz:,}")


# ------------------------------------------------------------------------- embedding
def _stats(X, rows):
    """mean and sd of each gene over `rows` of sparse X, without densifying."""
    n = len(rows)
    sub = X[rows]
    mu = np.asarray(sub.mean(axis=0)).ravel()
    sq = np.asarray(sub.multiply(sub).mean(axis=0)).ravel()
    var = np.maximum(sq - mu**2, 0.0) * (n / max(n - 1, 1))
    sd = np.sqrt(var)
    sd[sd < 1e-8] = 1.0
    return mu.astype(np.float32), sd.astype(np.float32)


def _iter_chunks(rows, size=40000):
    """Yield consecutive slices of rows of at most the given size."""
    for i in range(0, len(rows), size):
        yield rows[i : i + size]


def embed(X, meta, variant, ref_mask, chunk=40000):
    """Scale, fit PCA on ref_mask rows, return PCs for every row."""
    n, g = X.shape
    batches = pd.unique(meta.batch.values)
    ref_rows = np.where(ref_mask)[0]

    if variant == "shared":
        mu_b, sd_b = _stats(X, ref_rows)
        stats = {b: (mu_b, sd_b) for b in batches}
    elif variant == "perbatch":
        stats = {}
        for b in batches:
            rows = np.where(meta.batch.values == b)[0]
            stats[b] = _stats(X, rows)
    else:
        raise ValueError(variant)

    # pass 1: covariance over reference rows only
    C = np.zeros((g, g), dtype=np.float64)
    n_ref = 0
    for b in batches:
        rows = np.where((meta.batch.values == b) & ref_mask)[0]
        if len(rows) == 0:
            continue
        mu, sd = stats[b]
        for ch in _iter_chunks(rows, chunk):
            Z = X[ch].toarray()
            Z -= mu
            Z /= sd
            np.clip(Z, -CLIP, CLIP, out=Z)
            C += (Z.T @ Z).astype(np.float64)
            n_ref += len(ch)
    C /= max(n_ref - 1, 1)

    w, V = np.linalg.eigh(C)
    order = np.argsort(w)[::-1][:N_PCS]
    comps = V[:, order].astype(np.float32)
    evr = float(w[order].sum() / max(w.sum(), 1e-12))

    # pass 2: project everything
    PC = np.zeros((n, N_PCS), dtype=np.float32)
    for b in batches:
        rows = np.where(meta.batch.values == b)[0]
        mu, sd = stats[b]
        for ch in _iter_chunks(rows, chunk):
            Z = X[ch].toarray()
            Z -= mu
            Z /= sd
            np.clip(Z, -CLIP, CLIP, out=Z)
            PC[ch] = Z @ comps
    log(f"  variant={variant} PCA fitted on {n_ref:,} ref cells, {N_PCS} PCs explain {evr:.3f} of scaled variance")
    return PC, evr


# ------------------------------------------------------------------------------- knn
def knn_vote(ref_pc, ref_lab_codes, qry_pc, n_classes, k=K, chunk=256):
    """Exact brute-force kNN. Returns winning class code and vote share."""
    ref_sq = (ref_pc**2).sum(axis=1)
    nq = qry_pc.shape[0]
    pred = np.zeros(nq, dtype=np.int16)
    conf = np.zeros(nq, dtype=np.float32)
    for i in range(0, nq, chunk):
        q = qry_pc[i : i + chunk]
        d = ref_sq[None, :] - 2.0 * (q @ ref_pc.T)  # |q|^2 omitted, constant per row
        nn = np.argpartition(d, k, axis=1)[:, :k]
        lab = ref_lab_codes[nn]
        counts = np.zeros((lab.shape[0], n_classes), dtype=np.int16)
        for c in range(n_classes):
            counts[:, c] = (lab == c).sum(axis=1)
        best = counts.argmax(axis=1)
        pred[i : i + chunk] = best
        conf[i : i + chunk] = counts[np.arange(len(best)), best] / float(k)
    return pred, conf


def confusion(truth, pred, classes):
    """Cross-tabulate true and predicted labels over a fixed class list."""
    cm = pd.crosstab(
        pd.Categorical(truth, categories=classes),
        pd.Categorical(pred, categories=classes),
        dropna=False,
    )
    cm.index.name = "truth"
    cm.columns.name = "predicted"
    return cm


# ---------------------------------------------------------------------------- stages
def stage_validate():
    """Compare the shared and per-batch scaling variants by external (E-MTAB-6149) and leave-one-study-out accuracy; writes the validation summary."""
    X = sp.load_npz(os.path.join(CACHE, "lognorm.npz"))
    meta = pd.read_parquet(os.path.join(CACHE, "meta.parquet"))
    role = meta.role.values
    classes = sorted([c for c in meta.label.unique() if c != ""])
    code = {c: i for i, c in enumerate(classes)}
    log(f"classes: {classes}")

    rng = np.random.default_rng(SEED)
    results = []
    for variant in ["shared", "perbatch"]:
        log(f"=== variant {variant} ===")
        PC, evr = embed(X, meta, variant, ref_mask=(role == "ref"))
        np.save(os.path.join(CACHE, f"pc_{variant}.npy"), PC)

        ref_i = np.where(role == "ref")[0]
        ref_codes = np.array([code[c] for c in meta.label.values[ref_i]], dtype=np.int16)

        # --- external validation: EMTAB6149, full reference
        val_i = np.where(role == "val")[0]
        t0 = time.time()
        pred, conf = knn_vote(PC[ref_i], ref_codes, PC[val_i], len(classes))
        truth = meta.label.values[val_i]
        pl = np.array(classes)[pred]
        acc = float((pl == truth).mean())
        log(f"  external EMTAB6149 accuracy {acc:.4f}  ({time.time()-t0:.0f}s)")
        results.append(dict(variant=variant, test="external_EMTAB6149", n=len(val_i), accuracy=acc, evr=evr))
        np.save(os.path.join(CACHE, f"valpred_{variant}.npy"), pred)
        np.save(os.path.join(CACHE, f"valconf_{variant}.npy"), conf)

        # --- leave-one-study-out inside the reference
        for held in REF_STUDIES:
            tr = ref_i[meta.study.values[ref_i] != held]
            te = ref_i[meta.study.values[ref_i] == held]
            if len(te) > 25000:
                te = rng.choice(te, 25000, replace=False)
            tr_codes = np.array([code[c] for c in meta.label.values[tr]], dtype=np.int16)
            p, c = knn_vote(PC[tr], tr_codes, PC[te], len(classes))
            a = float((np.array(classes)[p] == meta.label.values[te]).mean())
            log(f"  LOSO hold out {held:<26} n={len(te):>6,}  accuracy {a:.4f}")
            results.append(dict(variant=variant, test=f"LOSO_{held}", n=len(te), accuracy=a, evr=evr))

    res = pd.DataFrame(results)
    res.to_csv(os.path.join(OUT, "labeltransfer_validation_summary.csv"), index=False)
    print(res.to_string(index=False))


def stage_apply(variant, floor):
    """Transfer the labels onto the query cells for one variant and confidence floor and write the per-cell and summary tables."""
    meta = pd.read_parquet(os.path.join(CACHE, "meta.parquet"))
    classes = sorted([c for c in meta.label.unique() if c != ""])
    code = {c: i for i, c in enumerate(classes)}
    PC = np.load(os.path.join(CACHE, f"pc_{variant}.npy"))
    role = meta.role.values
    ref_i = np.where(role == "ref")[0]
    qry_i = np.where(role == "qry")[0]
    ref_codes = np.array([code[c] for c in meta.label.values[ref_i]], dtype=np.int16)

    log(f"transferring onto {len(qry_i):,} query cells (variant={variant}, k={K})")
    pred, conf = knn_vote(PC[ref_i], ref_codes, PC[qry_i], len(classes))
    lab = np.array(classes)[pred]
    assigned = conf >= floor
    lab_final = np.where(assigned, lab, "")

    q = meta.iloc[qry_i].copy()
    q["pred"] = lab
    q["conf"] = conf
    q["assigned"] = assigned
    q["level1_compartment"] = lab_final
    q.to_parquet(os.path.join(CACHE, "query_assignment.parquet"))
    q[["cell_id", "donor_id", "histology", "sampling_site", "pred", "conf", "assigned", "level1_compartment"]].to_csv(
        os.path.join(OUT, "labeltransfer_emtab6653_percell.csv.gz"), index=False, compression="gzip"
    )

    print(f"\nassigned {int(assigned.sum()):,} of {len(q):,} ({assigned.mean()*100:.2f} pct); "
          f"left unassigned {int((~assigned).sum()):,} ({(~assigned).mean()*100:.2f} pct)")
    print("\n=== overall compartment breakdown (assigned cells) ===")
    vc = q.loc[assigned, "level1_compartment"].value_counts()
    print(pd.DataFrame({"n": vc, "pct": (vc / vc.sum() * 100).round(2)}).to_string())

    print("\n=== per patient (row pct of assigned) ===")
    ct = pd.crosstab(q.loc[assigned, "donor_id"], q.loc[assigned, "level1_compartment"])
    print(ct.to_string())
    print((ct.div(ct.sum(axis=1), axis=0) * 100).round(2).to_string())
    ct.to_csv(os.path.join(OUT, "labeltransfer_emtab6653_by_patient.csv"))

    print("\n=== per sampling site (row pct of assigned) ===")
    site_order = ["core", "middle", "edge", "normal"]
    cs = pd.crosstab(q.loc[assigned, "sampling_site"], q.loc[assigned, "level1_compartment"]).reindex(site_order)
    print(cs.to_string())
    print((cs.div(cs.sum(axis=1), axis=0) * 100).round(2).to_string())
    cs.to_csv(os.path.join(OUT, "labeltransfer_emtab6653_by_site.csv"))

    print("\n=== per patient x site, epithelial pct ===")
    pt = pd.crosstab([q.loc[assigned, "donor_id"], q.loc[assigned, "sampling_site"]],
                     q.loc[assigned, "level1_compartment"])
    eps = (pt["epithelial_malignant"] / pt.sum(axis=1) * 100).round(2).unstack().reindex(columns=site_order)
    print(eps.to_string())
    eps.to_csv(os.path.join(OUT, "labeltransfer_emtab6653_epithelial_pct.csv"))

    print("\n=== unassigned rate per patient / site ===")
    print(q.groupby("donor_id")["assigned"].agg(["size", "mean"]).assign(unassigned_pct=lambda d: (1-d["mean"])*100).to_string())
    print(q.groupby("sampling_site")["assigned"].agg(["size", "mean"]).assign(unassigned_pct=lambda d: (1-d["mean"])*100).to_string())

    print("\n=== sanity: epithelial fraction core vs adjacent normal ===")
    for d in sorted(q.donor_id.unique()):
        sub = q[assigned & (q.donor_id == d)]
        core = (sub[sub.sampling_site == "core"].level1_compartment == "epithelial_malignant").mean() * 100
        norm = (sub[sub.sampling_site == "normal"].level1_compartment == "epithelial_malignant").mean() * 100
        print(f"  {d}: core {core:.2f} pct  vs  adjacent normal {norm:.2f} pct   {'OK' if core>norm else 'INVERTED'}")


def stage_write(variant, floor):
    """Write the transferred labels into panel_object.h5ad (backup first) and verify the written object."""
    import anndata as ad

    q = pd.read_parquet(os.path.join(CACHE, "query_assignment.parquet"))
    backup = os.path.join(OUT, "panel_object_prelabel.h5ad")
    if not os.path.exists(backup):
        log("copying panel_object.h5ad -> panel_object_prelabel.h5ad")
        import shutil

        shutil.copy2(OBJ, backup)
        check(os.path.getsize(backup) == os.path.getsize(OBJ), "backup byte size matches original")
    else:
        log("backup already exists, not overwriting it")

    log("reading object for update")
    A = ad.read_h5ad(OBJ)
    o = A.obs
    pos = pd.Index(o.index).get_indexer(q.cell_id.values)
    check((pos >= 0).all(), "every query cell id found in the object")
    check(len(set(pos)) == len(pos), "query cell ids map to distinct rows")

    l1 = o.level1_compartment.astype(str).values.copy()
    # The guard requires that every target row is empty or was written by this same transfer
    # process, never by an author or the census.
    _src = A.obs["celltype_source"].astype(str).values[pos]
    _ok = ((l1[pos] == "") | pd.Series(_src).isin(["transferred", "unassigned"]).values).all()
    check(bool(_ok), "target rows are empty or previously transferred, safe to overwrite")
    l1[pos] = q.level1_compartment.values
    cats = list(A.obs.level1_compartment.cat.categories)
    o["level1_compartment"] = pd.Categorical(l1, categories=cats)

    src = o.celltype_source.astype(str).values.copy()
    src[pos] = np.where(q.assigned.values, "transferred", "unassigned")
    o["celltype_source"] = pd.Categorical(src, categories=sorted(set(src)))

    cc = o.celltype_confidence.values.copy()
    cc[pos] = q.conf.values
    o["celltype_confidence"] = cc

    # level2 / level3 stay empty by design
    check(o.level2_celltype.astype(str).values[pos].tolist() == [""] * len(pos), "level2 left empty")
    check(o.level3_celltype_native.astype(str).values[pos].tolist() == [""] * len(pos), "level3 left empty")

    A.uns["label_transfer_emtab6653"] = json.dumps(
        {
            "target": "EMTAB6653, 32,341 cells, level1_compartment only",
            "reference": "LuCA census, 10x 3' v2, whole_tissue, labelled; studies "
            + ", ".join(REF_STUDIES),
            "method": f"normalize_total(1e4)+log1p, {variant} gene scaling clipped at {CLIP}, "
            f"PCA({N_PCS}) fitted on the reference, exact kNN k={K} vote in PC space",
            "confidence": f"fraction of the {K} nearest reference cells carrying the winning label",
            "confidence_floor": floor,
            "validation": "held out EMTAB6149 (10x 3' v2, same laboratory and protocol as the target, "
            "donors P3/P4/P5 disjoint from target donors P6/P7/P8, author labels)",
            "unassigned_encoding": "level1_compartment empty and celltype_source 'unassigned'; "
            "celltype_confidence is stored for every cell so the floor can be revisited",
            "level2_level3": "deliberately left empty, no fine labels were transferred",
            "script": "atlas/steps/stage1b_label_transfer_emtab6653.py",
        }
    )
    A.uns["unlabelled_cells"] = (
        "EMTAB6653 (32,341 cells) had no author cell types. level1_compartment was filled by kNN "
        "transfer from the LuCA 10x 3' v2 whole_tissue reference (see uns.label_transfer_emtab6653); "
        f"{int((~q.assigned).sum()):,} cells fell below the confidence floor and remain empty with "
        "celltype_source 'unassigned'. level2_celltype and level3_celltype_native are empty for the "
        "whole arm. Compartment queries must still filter level1_compartment != ''."
    )

    tmp = OBJ + ".tmp"
    log("writing object")
    A.write_h5ad(tmp, compression="gzip")
    os.replace(tmp, OBJ)
    log("wrote " + OBJ)

    B = ad.read_h5ad(OBJ, backed="r")
    ob = B.obs
    m = ob.source_key.astype(str) == "EMTAB6653"
    check(B.shape == (1146594, 1906), "round trip shape", str(B.shape))
    check(int(m.sum()) == 32341, "round trip arm size")
    n_assigned = int((ob.level1_compartment.astype(str)[m] != "").sum())
    check(n_assigned == int(q.assigned.sum()), "round trip assigned count", f"{n_assigned}")
    check(
        (ob.celltype_source.astype(str)[m & (ob.level1_compartment.astype(str) != "")] == "transferred").all(),
        "assigned cells carry celltype_source transferred",
    )
    check((ob.level2_celltype.astype(str)[m] == "").all(), "level2 still empty for the arm")
    check((ob.level3_celltype_native.astype(str)[m] == "").all(), "level3 still empty for the arm")
    check(
        int((ob.level1_compartment.astype(str) != "").sum())
        == 1146594 - 32341 + int(q.assigned.sum()),
        "total labelled cells increased by exactly the assigned count",
    )
    print("\nfinal level1_compartment:")
    print(ob.level1_compartment.value_counts(dropna=False).to_string())


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["build", "validate", "apply", "write"])
    ap.add_argument("--variant", default="perbatch")
    ap.add_argument("--floor", type=float, default=0.5)
    a = ap.parse_args()
    if a.stage == "build":
        build_groups()
    elif a.stage == "validate":
        stage_validate()
    elif a.stage == "apply":
        stage_apply(a.variant, a.floor)
    elif a.stage == "write":
        stage_write(a.variant, a.floor)
