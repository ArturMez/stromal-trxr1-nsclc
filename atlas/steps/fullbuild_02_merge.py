#!/usr/bin/env python3
"""
Atlas stage 04, step merge: merge the arms onto a union gene axis with a presence mask.

Paper items: writes fullgene_object.h5ad (1,888,910 cells x 50,552 genes), the atlas behind Figure 2,
Figure 3, Supplementary Figures 2 to 4 and Supplementary Data 1. No random numbers are drawn.

The full-transcriptome object covers every patient in `panel_object.h5ad`, for fibroblasts, myeloid and
T/NK cells. Seven arms feed it, built by the extract steps (extract_luca to extract_wave1).

THE THING THAT MAKES THIS HONEST: THE PRESENCE MASK. The arms measure different gene sets. LuCA has
17,811, E-MTAB-13526 has 33,514, E-MTAB-6653 has 23,968, E-MTAB-6149 has its own. On a union axis a
gene an arm never measured becomes a zero that looks exactly like a measured zero, and any analysis
that pools across arms will read "not measured" as "not expressed". So `varm/gene_present` is written:
a boolean gene-by-arm matrix saying which arm measured which gene, together with `uns/arm_gene_counts`.

The rule that makes this workable: every analysis is computed INSIDE a block and combined by block
votes, so an arm does not need to share the gene axis, it needs to carry the gene sets being scored.
The mask is what lets a scoring script check that instead of assuming it. **Any script scoring a gene
set must drop blocks whose arm does not measure the whole set, and say how many it dropped.**

TWO FLAGS TRAVEL WITH THE CELLS AND MUST NOT BE DROPPED.
  counts_are_estimates   True for every E-MTAB-6653 cell (Alevin EM estimates) and for LuCA's 17,891
                         Smart-seq2 cells. Anything needing true integer UMI counts must filter on it.
  block_qc_pass          the existing per-block gene-detection QC, carried unchanged from the panel
                         object so this object answers to the same QC as the panel object.

MEMORY. Arms are appended one at a time and each arm is streamed in row chunks, so peak memory is one
chunk rather than one arm. Gene indices are remapped per arm through a lookup into the union axis.

RUN: python3 atlas/steps/fullbuild_02_merge.py   (normally through atlas/04_build_objects.py)
FROM: the seven arm_*_full.h5ad files of the extract steps (the environment variable FULLGENE_DST
      overrides the output path).
OUT -> output/RNA/meta_dataset/fullbuild/fullgene_object.h5ad
       output/RNA/meta_dataset/fullbuild/fullgene_gene_presence.csv
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
DST = os.environ.get("FULLGENE_DST", f"{D}/fullgene_object.h5ad")
# Every arm built by the extract steps. The first four cover the patients in panel_object; the last
# three are the wave-1 sources. GSE136246 is not part of the atlas: normalized counts, shallow
# inDrop, and 7 of its 12 donors duplicate Zilionis_Klein_2019.
ARMS = [("LuCA", f"{D}/arm_luca_full.h5ad"),
        ("E-MTAB-13526", f"{D}/arm_emtab13526_full.h5ad"),
        ("E-MTAB-6653", f"{D}/arm_emtab6653_full.h5ad"),
        ("E-MTAB-6149", f"{D}/arm_emtab6149_full.h5ad"),
        ("STTT", f"{D}/arm_sttt_full.h5ad"),
        ("GSE325414", f"{D}/arm_gse325414_full.h5ad"),
        ("GSE246215", f"{D}/arm_gse246215_full.h5ad")]
CHUNK = 20000


def dec(a):
    """Decode an array of bytes or strings to an array of str."""
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in a])


def set_colorder(grp, cols):
    """Write the anndata `column-order` attribute (an empty list needs an explicit empty array)."""
    if cols:
        grp.attrs["column-order"] = np.array(list(cols), dtype=h5py.string_dtype())
    else:
        grp.attrs.create("column-order", np.array([], dtype="S1"))


present = [(n, p) for n, p in ARMS if os.path.exists(p)]
missing = [n for n, p in ARMS if not os.path.exists(p)]
if missing:
    print(f"NOTE: arms not yet built and therefore not merged: {missing}", flush=True)
if not present:
    raise SystemExit("no arms found")

print("reading arm gene axes ...", flush=True)
axes, ncells, obscols = {}, {}, None
for n, p in present:
    with h5py.File(p, "r") as h:
        axes[n] = dec(h["var/gene_symbol"][:])
        ncells[n] = int(h["X"].attrs["shape"][0])
        cols = set(h["obs"].keys()) - {"_index"}
        # UNION, not intersection. An intersection would drop block_qc_pass from the whole object,
        # because the three wave-1 arms do not carry it (block_qc_pass excludes 131,092
        # Adams_Kaminski cells, so losing it would change every analysis).
        obscols = cols if obscols is None else (obscols | cols)
    print(f"  {n:14s} {ncells[n]:>9,} cells  {len(axes[n]):>6,} genes", flush=True)
obscols = sorted(obscols)
union = np.array(sorted(set().union(*axes.values())))
gpos = pd.Series(np.arange(len(union)), index=union)
total = sum(ncells.values())
print(f"\n  union gene axis {len(union):,}   total cells {total:,}", flush=True)

mask = pd.DataFrame({"gene": union})
for n, _ in present:
    mask[n] = mask.gene.isin(set(axes[n]))
mask["n_arms"] = mask[[n for n, _ in present]].sum(axis=1)
mask.to_csv(f"{D}/fullgene_gene_presence.csv", index=False)
print("  genes measured by all arms: %d   by only one: %d"
      % (int((mask.n_arms == len(present)).sum()), int((mask.n_arms == 1).sum())), flush=True)

BOOLFILL = {"block_qc_pass", "composition_safe", "counts_are_estimates"}
PANEL_ARMS = {"LuCA", "E-MTAB-13526", "E-MTAB-6653", "E-MTAB-6149"}
obs_acc = {c: [] for c in obscols}
assessed = []
ids_acc = []
ptr = [0]
nnz = 0
with h5py.File(DST, "w") as o:
    g = o.create_group("X")
    d_out = g.create_dataset("data", (0,), maxshape=(None,), dtype="float32", chunks=(1 << 20,))
    i_out = g.create_dataset("indices", (0,), maxshape=(None,), dtype="int32", chunks=(1 << 20,))
    for n, p in present:
        remap = gpos.reindex(axes[n]).values.astype(np.int32)
        with h5py.File(p, "r") as h:
            ip = h["X/indptr"][:]
            nc = len(ip) - 1
            ids_acc.append(dec(h["obs/_index"][:]))
            for c in obscols:
                if c not in h["obs"]:
                    # An arm that lacks a column gets an explicit fill, and the fill for
                    # block_qc_pass is False: a block that has not been through the block QC
                    # has not passed it. block_qc_assessed records which is which (stage 04
                    # transfer_qc later writes the block QC of the 24 blocks onto these cells).
                    obs_acc[c].append(np.array(["" if c not in BOOLFILL else "False"] * nc,
                                               dtype=object))
                    continue
                x = h[f"obs/{c}"]
                if isinstance(x, h5py.Group):
                    if "categories" in x:
                        obs_acc[c].append(dec(x["categories"][:])[x["codes"][:]])
                    else:
                        obs_acc[c].append(x["values"][:])
                else:
                    v = x[:]
                    obs_acc[c].append(dec(v) if v.dtype.kind in "SO" else v)
            assessed.append(np.array([n in PANEL_ARMS] * nc, bool))
            for a in range(0, nc, CHUNK):
                b = min(a + CHUNK, nc)
                lo, hi = int(ip[a]), int(ip[b])
                dat = h["X/data"][lo:hi]
                ind = remap[h["X/indices"][lo:hi]]
                # each row's indices must stay ascending after the remap
                loc = ip[a:b + 1] - lo
                for r in range(b - a):
                    s0, s1 = int(loc[r]), int(loc[r + 1])
                    if s1 > s0:
                        srt = np.argsort(ind[s0:s1], kind="stable")
                        ind[s0:s1] = ind[s0:s1][srt]
                        dat[s0:s1] = dat[s0:s1][srt]
                    ptr.append(ptr[-1] + (s1 - s0))
                d_out.resize((nnz + len(dat),))
                i_out.resize((nnz + len(ind),))
                d_out[nnz:] = dat.astype("float32")
                i_out[nnz:] = ind.astype("int32")
                nnz += len(dat)
            print(f"  merged {n:14s} {nc:>9,} cells, running nonzeros {nnz:,}", flush=True)
    g.create_dataset("indptr", data=np.asarray(ptr, "int64"))
    g.attrs["encoding-type"] = "csr_matrix"
    g.attrs["encoding-version"] = "0.1.0"
    g.attrs["shape"] = np.array([total, len(union)], "int64")

    v = o.create_group("var")
    v.attrs["_index"] = "gene_symbol"
    v.attrs["encoding-type"] = "dataframe"
    v.attrs["encoding-version"] = "0.2.0"
    set_colorder(v, [])
    v.create_dataset("gene_symbol", data=np.array(union, dtype=object), dtype=h5py.string_dtype())

    vm = o.create_group("varm")
    gp = vm.create_group("gene_present")
    gp.attrs["encoding-type"] = "dataframe"
    gp.attrs["encoding-version"] = "0.2.0"
    gp.attrs["_index"] = "gene_symbol"
    set_colorder(gp, [n for n, _ in present])
    gp.create_dataset("gene_symbol", data=np.array(union, dtype=object), dtype=h5py.string_dtype())
    for n, _ in present:
        gp.create_dataset(n, data=mask[n].values.astype(bool))

    ob = o.create_group("obs")
    ob.attrs["_index"] = "_index"
    ob.attrs["encoding-type"] = "dataframe"
    ob.attrs["encoding-version"] = "0.2.0"
    set_colorder(ob, sorted(obscols) + ["block_qc_assessed"])
    ob.create_dataset("_index", data=np.array(np.concatenate(ids_acc), dtype=object),
                      dtype=h5py.string_dtype())
    ob.create_dataset("block_qc_assessed", data=np.concatenate(assessed))
    for c in obscols:
        x = np.concatenate([np.asarray(a, dtype=object) if np.asarray(a).dtype.kind in "OU"
                            else np.asarray(a) for a in obs_acc[c]])
        if x.dtype == object or x.dtype.kind in "OU":
            ob.create_dataset(c, data=np.array(x.astype(str), dtype=object),
                              dtype=h5py.string_dtype())
        else:
            ob.create_dataset(c, data=x)
    o.attrs["encoding-type"] = "anndata"
    o.attrs["encoding-version"] = "0.1.0"

print(f"\nwritten -> {DST}")
print(f"  cells {total:,}   genes {len(union):,}   nonzeros {nnz:,}   "
      f"{os.path.getsize(DST)/1e9:.1f} GB on disk")
print("\n" + "=" * 100)
print("GENE PRESENCE BY ARM. A zero for a gene an arm never measured is NOT a measured zero.")
print("VARIABLE: whether that arm's deposit contains the gene at all. UNIT: gene.")
print("DENOMINATOR: the %d genes on the union axis." % len(union))
print("=" * 100)
print(mask[[n for n, _ in present]].sum().to_string())
print(f"\nmeasured by every arm: {int((mask.n_arms == len(present)).sum()):,}")
print(f"measured by only one arm: {int((mask.n_arms == 1).sum()):,}")
