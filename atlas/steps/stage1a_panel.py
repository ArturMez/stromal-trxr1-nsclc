#!/usr/bin/env python
"""
Atlas stage 01, step panel: define the gene panel ONCE, so all three source arms subset identically.

Paper items: the 1,906-gene panel on which the atlas arms are first built (the atlas behind Figure 2,
Figure 3, Supplementary Figures 2 to 4 and Supplementary Data 1). Random numbers: the HVG cell sample
is drawn with numpy default_rng(SEED), SEED = 0 (contiguous runs of BLOCK = 25 cells, at most
CELL_CAP_PER_DATASET = 6,000 cells per LuCA sub-dataset); the 2,000 highly variable genes use scanpy
seurat_v3 with batch_key = dataset.

Run: python atlas/steps/stage1a_panel.py   (normally through atlas/01_arms_luca_lambrechts.py)

Inputs (read only)
  LuCA/luca_extended.h5ad                                       reference var index + HVG source
  output/RNA/meta_dataset/_features6149.csv                     E-MTAB-6149 symbol -> Ensembl table
  E-MTAB-6653 matrix rows (Ensembl ids)                         for the cross-source presence check

Output
  output/RNA/meta_dataset/stage1a/panel_genes.csv
      gene_symbol, ensembl_id, source_set, is_hvg

Two traps this script defends against
  1. LuCA feature_name is NOT a clean symbol. 22 genes carry a disambiguated
     "<SYMBOL>_<ENSG>" name because the symbol maps to more than one Ensembl id.
     RGS5 is one of them. A naive symbol match drops RGS5, which is the marker the
     level1 mural bucket and the vCAF signature both rest on, and which exposes the
     mislabelled Lambrechts Fibro2 cluster.
     Every symbol is therefore normalized and every gene carries its Ensembl id.
  2. The HVG sample is drawn from the cells that will actually enter the merged object:
     is_primary_data == True and NOT the three Lambrechts datasets (those 66,296 cells
     come in from the local E-MTAB files instead). HVGs computed on cells the object
     will not contain would be a different gene panel from the one the object supports.

No matrix is built here. Panel definition only.
"""

import os
import re
import sys
import numpy as np
import pandas as pd
import scipy.sparse as sp

import os as _os  # noqa: E402
import sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # atlas/
from common import ROOT, SCRNA, SCRATCH, DATA, RSCRIPT  # noqa: E402,F401  paths from config.json
LUCA = f"{SCRNA}/LuCA/luca_extended.h5ad"
OUT = f"{ROOT}/output/RNA/meta_dataset/stage1a"
os.makedirs(OUT, exist_ok=True)

N_HVG = 2000
CELL_CAP_PER_DATASET = 6000     # HVG sample cap, per LuCA sub-dataset
BLOCK = 25                      # sample in short contiguous runs so the CSR reads stay cheap
SEED = 0
LAMBRECHTS_DATASETS = (
    "Lambrechts_Thienpont_2018_6149v1",
    "Lambrechts_Thienpont_2018_6149v2",
    "Lambrechts_Thienpont_2018_6653",
)

# ------------------------------------------------------------------------------------
# 1. GENE SETS
#    The curated gene sets that the atlas analyses score (redox, T cell, myeloid, CAF,
#    histology identity and compartment markers); they are added to the 2,000 HVGs.
# ------------------------------------------------------------------------------------

GENE_SETS = {
    # --- redox -----------------------------------------------------------------------
    # core redox panel
    "redox_core": [
        "TXNRD1", "NFE2L2", "KEAP1", "CUL3", "SQSTM1", "NQO1", "GPX3", "HMOX1",
        "TXN", "TXNIP", "SLC7A11", "SOD3", "GCLC", "GCLM", "PRDX1",
    ],
    # thioredoxin program
    "redox_thioredoxin": [
        "TXN", "TXN2", "TXNRD1", "TXNRD2",
        "PRDX1", "PRDX2", "PRDX3", "PRDX5", "PRDX6", "SRXN1",
    ],
    # glutathione program
    "redox_glutathione": [
        "GCLC", "GCLM", "GSR", "GSS", "GPX1", "GPX4", "SLC7A11",
    ],
    # NRF2 target program
    "redox_nrf2_target": [
        "NQO1", "HMOX1", "SRXN1", "AKR1B10", "AKR1C1", "G6PD", "GCLM",
        "TXNRD1", "GPX2",
    ],
    # NRF2 pathway genes
    "redox_nrf2_pathway": ["NFE2L2", "KEAP1", "CUL3", "SQSTM1"],
    # superoxide dismutases
    "redox_superoxide_dismutase": ["SOD1", "SOD2", "SOD3"],
    # ROS-generating enzymes
    "redox_ros_source": ["NOX4", "CYBB", "NCF1", "NOS2"],

    # --- T cell programs ---------------------------------------------------
    "tcell_exhaustion": ["PDCD1", "HAVCR2", "LAG3", "TOX", "CTLA4", "TIGIT", "ENTPD1"],
    "tcell_activation": ["IL2RA", "CD38", "HLA-DRA", "TNFRSF9", "ENTPD1", "ICOS", "CD69"],
    "tcell_cytotoxicity": ["GZMB", "GZMA", "PRF1", "IFNG", "NKG7", "GNLY"],
    "tcell_naive_memory": ["TCF7", "SELL", "IL7R", "CCR7", "LEF1"],
    "treg_identity": ["FOXP3", "IKZF2", "TNFRSF18", "CTLA4", "TIGIT", "IL2RA"],
    # proliferation
    "proliferation": ["MKI67", "TOP2A", "PCNA", "BIRC5", "TYMS", "STMN1", "PCLAF"],

    # --- myeloid polarization ------------------------------------------------
    "m1_pro_oxidant": ["NOS2", "CYBB", "IL1B", "TNF", "CXCL9", "CXCL10", "SOD2", "NCF1"],
    "m2_reductive": ["CD163", "MRC1", "MSR1", "IL10", "HMOX1", "FTL", "FTH1", "SLC40A1"],

    # --- CAF subtypes --------------------------------------------------------
    "caf_myCAF": ["ACTA2", "TAGLN", "MYL9", "TPM1", "TPM2", "POSTN", "CTHRC1",
                  "MMP11", "COL1A1", "COL1A2"],
    "caf_iCAF": ["IL6", "CXCL12", "CXCL1", "CXCL2", "CXCL14", "CFD", "C3", "C7",
                 "PDGFRA", "LMNA", "HAS1", "IL1B"],
    "caf_apCAF": ["CD74", "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "SLPI", "CD83"],
    "caf_vCAF": ["RGS5", "NOTCH3", "PDGFRB", "MCAM", "ACTA2", "KCNJ8", "HIGD1B"],
    "myofibroblast": ["ACTA2", "TAGLN", "COL1A1", "COL1A2", "FAP", "PDGFRB", "POSTN"],
    # single genes tested against TXNRD1 in fibroblasts
    "fibroblast_single_genes": ["NOX4", "ACTA2", "CD74", "IL6", "CXCL12", "PDGFRA",
                                "RGS5", "FAP", "POSTN", "TGFB1"],
    # the four markers that expose Lambrechts Fibro2 as pericyte rather than fibroblast. 506 of 1,465 fibroblasts depend on this.
    "pericyte_disambiguation": ["COX4I2", "HIGD1B", "RGS5", "NOTCH3"],

    # --- histology identity ------------------------------------------------------
    "squamous_identity": ["KRT5", "KRT6A", "KRT14", "TP63", "SOX2", "PKP1", "DSG3"],
    "adeno_identity": ["NKX2-1", "NAPSA", "SFTPC", "SFTPB", "MUC1"],

    # --- compartment markers (level1_compartment of the obs schema) ------------
    "compartment_marker": [
        "PTPRC", "EPCAM", "COL1A1", "PDGFRA", "PDGFRB", "FAP", "ACTA2", "RGS5",
        "NOTCH3", "PECAM1", "VWF", "CD68", "CD163", "LYZ", "CD3E", "CD8A", "CD4",
        "FOXP3", "NKG7", "MS4A1", "JCHAIN", "TPSAB1",
    ],
}

# the core genes that must be confirmed against the LuCA var index
BRIEF_NAMED = (
    GENE_SETS["redox_core"]
    + GENE_SETS["squamous_identity"]
    + GENE_SETS["adeno_identity"]
    + GENE_SETS["compartment_marker"]
)

# Ensembl ids for genes absent from the LuCA var index, so the other two arms can still
# subset by a stable identifier. Resolved from the E-MTAB-6149 featuredata table
# (Associated.Gene.Name -> Ensembl rowname), not from memory.
FALLBACK_ENSEMBL = {
    "GPX1": "ENSG00000233276",
}

# gene tables of the other two arms, for the cross-source presence audit
F6149 = f"{SCRNA}/E-MTAB-6149/processed/AllSamplesCelllview.R"
F6653_ROWS = (f"{SCRNA}/E-MTAB-6653/raw_counts/matrix_market/"
              "E-MTAB-6653.aggregated_filtered_counts.mtx_rows")
FEAT6149_CACHE = f"{OUT}/_features6149.csv"   # symbol -> Ensembl, dumped once via R


# ------------------------------------------------------------------------------------
# 2. LuCA var index
# ------------------------------------------------------------------------------------

_MANGLED = re.compile(r"^(?P<sym>.+)_(?P<ens>ENSG\d{11})$")


def clean_symbol(name: str) -> str:
    """CELLxGENE appends _<ENSG> when a symbol maps to several gene ids. Strip it."""
    m = _MANGLED.match(name)
    return m.group("sym") if m else name


def read_luca_var():
    """Read the LuCA var index (Ensembl id, feature name) and add the cleaned gene symbol."""
    import h5py
    with h5py.File(LUCA, "r") as f:
        v = f["var"]
        gene_id = np.asarray(v["gene_id"]).astype(str)
        cats = np.asarray(v["feature_name"]["categories"]).astype(str)
        codes = np.asarray(v["feature_name"]["codes"])
        feature_name = cats[codes]
    df = pd.DataFrame({"gene_id": gene_id, "feature_name": feature_name})
    df["gene_symbol"] = [clean_symbol(x) for x in feature_name]
    dup = df["gene_symbol"].duplicated(keep=False)
    assert not dup.any(), (
        "cleaning the disambiguated feature_name collapsed two gene ids onto one "
        f"symbol: {df.loc[dup].to_dict('records')}")
    return df


# ------------------------------------------------------------------------------------
# 3. HVG on a sample of the LuCA atlas
# ------------------------------------------------------------------------------------

def sample_row_indices():
    """Cells eligible for the merged object, sampled in short contiguous runs."""
    import h5py
    with h5py.File(LUCA, "r") as f:
        o = f["obs"]
        prim = np.asarray(o["is_primary_data"])
        dcats = np.asarray(o["dataset"]["categories"]).astype(str)
        dcodes = np.asarray(o["dataset"]["codes"])
        dncats = np.asarray(o["donor_id"]["categories"]).astype(str)
        dncodes = np.asarray(o["donor_id"]["codes"])
        discats = np.asarray(o["disease"]["categories"]).astype(str)
        discodes = np.asarray(o["disease"]["codes"])
    lam = np.isin(dcats, LAMBRECHTS_DATASETS)
    eligible = prim & ~lam[dcodes]

    # the Census filter, asserted on the donor_id prefix
    by_donor = prim & ~np.char.startswith(dncats, "Lambrechts")[dncodes]
    assert np.array_equal(eligible, by_donor), (
        "the Lambrechts dataset-name filter and the donor_id-prefix filter disagree")
    assert int(prim.sum() - eligible.sum()) == 66296, (
        f"Lambrechts exclusion removed {int(prim.sum() - eligible.sum())}, expected 66296")
    assert int(eligible.sum()) == 1061555, (
        f"{int(eligible.sum())} cells remain, expected 1061555")
    d = pd.DataFrame({"donor": dncats[dncodes][eligible],
                      "disease": discats[discodes][eligible]})
    nd = d.groupby("disease")["donor"].nunique()
    assert int(nd.get("squamous cell lung carcinoma", 0)) == 42, nd.to_dict()
    assert int(nd.get("lung adenocarcinoma", 0)) == 155, nd.to_dict()
    rng = np.random.default_rng(SEED)
    picked = []
    for ci, name in enumerate(dcats):
        if lam[ci]:
            continue
        idx = np.where(eligible & (dcodes == ci))[0]
        if idx.size == 0:
            continue
        if idx.size <= CELL_CAP_PER_DATASET:
            picked.append(idx)
            continue
        n_blocks = int(np.ceil(CELL_CAP_PER_DATASET / BLOCK))
        starts = rng.choice(idx.size - BLOCK, size=n_blocks, replace=False)
        picked.append(np.concatenate([idx[s:s + BLOCK] for s in starts]))
    rows = np.unique(np.concatenate(picked))
    return rows, dcodes, dcats, int(prim.sum()), int(eligible.sum())


def read_counts_rows(rows, n_genes):
    """Read layers/count for the given rows, span by span, and rebuild a CSR."""
    import h5py
    with h5py.File(LUCA, "r") as f:
        g = f["layers"]["count"]
        indptr = np.asarray(g["indptr"])
        D, I = g["data"], g["indices"]
        breaks = np.where(np.diff(rows) != 1)[0]
        starts = np.concatenate([[0], breaks + 1])
        stops = np.concatenate([breaks + 1, [rows.size]])
        data_parts, ind_parts = [], []
        nnz_per_row = np.empty(rows.size, dtype=np.int64)
        for a, b in zip(starts, stops):
            r0, r1 = int(rows[a]), int(rows[b - 1]) + 1
            p0, p1 = int(indptr[r0]), int(indptr[r1])
            data_parts.append(np.asarray(D[p0:p1], dtype=np.float32))
            ind_parts.append(np.asarray(I[p0:p1], dtype=np.int32))
            nnz_per_row[a:b] = np.diff(indptr[r0:r1 + 1])
    ptr = np.zeros(rows.size + 1, dtype=np.int64)
    np.cumsum(nnz_per_row, out=ptr[1:])
    X = sp.csr_matrix(
        (np.concatenate(data_parts), np.concatenate(ind_parts), ptr),
        shape=(rows.size, n_genes),
    )
    return X


def compute_hvg(var):
    """Select the 2,000 highly variable genes (seurat_v3, batch_key = dataset) on the sampled LuCA cells; returns the Ensembl ids, cells and datasets used."""
    import anndata as ad
    import scanpy as sc

    rows, dcodes, dcats, n_prim, n_elig = sample_row_indices()
    print(f"LuCA cells 1,283,972 | is_primary_data True {n_prim:,} | "
          f"eligible after dropping Lambrechts {n_elig:,}", flush=True)
    print(f"HVG sample: {rows.size:,} cells from "
          f"{len(set(dcats[dcodes[rows]]))} datasets", flush=True)

    X = read_counts_rows(rows, len(var))

    # layers/count is NOT integer everywhere. The two Smart-seq2 sub-atlases
    # (Maynard_Bivona_2020, Guo_Zhang_2018) ship fractional EM abundance estimates,
    # the same trap the obs schema records for E-MTAB-6653. seurat_v3 fits a
    # variance-stabilising loess that assumes count-like input, so the fractional
    # entries are rounded here, for the HVG fit only. Nothing downstream reads this.
    ds_of_row = dcats[dcodes[rows]]
    nonint = X.data != np.floor(X.data)
    if nonint.any():
        per_row = np.add.reduceat(nonint, X.indptr[:-1]) if X.indptr[-1] else None
        aff = pd.Series(per_row, index=ds_of_row).groupby(level=0).sum()
        aff = aff[aff > 0]
        print(f"layers/count non-integer entries: {int(nonint.sum()):,} of "
              f"{X.data.size:,} ({100 * nonint.mean():.2f} pct), confined to "
              f"{', '.join(aff.index)}; rounded for the HVG fit only", flush=True)
        np.round(X.data, out=X.data)
        X.eliminate_zeros()

    a = ad.AnnData(X=X, obs=pd.DataFrame({"dataset": pd.Categorical(dcats[dcodes[rows]])},
                                         index=[f"c{i}" for i in rows]),
                   var=pd.DataFrame(index=var["gene_id"].values))
    sc.pp.highly_variable_genes(a, n_top_genes=N_HVG, flavor="seurat_v3",
                                batch_key="dataset")
    hv = a.var.index[a.var["highly_variable"].values].tolist()
    assert len(hv) == N_HVG, f"expected {N_HVG} HVGs, got {len(hv)}"
    return set(hv), rows.size, len(set(dcats[dcodes[rows]]))


# ------------------------------------------------------------------------------------
# 4. cross-source presence
# ------------------------------------------------------------------------------------

def presence_audit(panel):
    """Which of the three arms can actually supply each panel gene.

    E-MTAB-6149 is keyed on Ensembl with an Associated.Gene.Name column; only 22,180 of
    its 55,457 annotated features survive into the log2cpm matrix. E-MTAB-6653 carries
    Ensembl ids and no symbol column at all, so it can only be joined on ensembl_id.
    """
    f49 = pd.read_csv(FEAT6149_CACHE)
    in49 = set(f49.loc[f49["in_matrix"], "ensembl_id"])
    sym49 = dict(zip(f49.loc[f49["in_matrix"], "gene_symbol"],
                     f49.loc[f49["in_matrix"], "ensembl_id"]))
    r53 = pd.read_csv(F6653_ROWS, sep="\t", header=None, usecols=[0])
    in53 = set(r53[0].astype(str))

    out = panel.copy()
    out["in_LuCA"] = out["ensembl_id"].isin(set(read_luca_var()["gene_id"]))
    out["in_EMTAB6149"] = [
        (e in in49) or (s in sym49) for e, s in zip(out["ensembl_id"], out["gene_symbol"])
    ]
    out["in_EMTAB6653"] = out["ensembl_id"].isin(in53)
    out["n_arms"] = out[["in_LuCA", "in_EMTAB6149", "in_EMTAB6653"]].sum(axis=1)
    return out[["gene_symbol", "ensembl_id", "source_set", "is_hvg",
                "in_LuCA", "in_EMTAB6149", "in_EMTAB6653", "n_arms"]]


# ------------------------------------------------------------------------------------
# 5. assemble
# ------------------------------------------------------------------------------------

def main():
    """Assemble the panel (curated sets plus HVGs), check presence in the three arms, and write panel_genes.csv and panel_genes_source_presence.csv."""
    var = read_luca_var()
    sym2ens = dict(zip(var["gene_symbol"], var["gene_id"]))
    mangled = var[var["feature_name"] != var["gene_symbol"]]

    # membership check against the LuCA var index, on cleaned symbols
    curated = {}
    for setname, genes in GENE_SETS.items():
        for g in genes:
            curated.setdefault(g, []).append(setname)

    missing = sorted(g for g in curated if g not in sym2ens)
    missing_named = sorted(g for g in BRIEF_NAMED if g not in sym2ens)

    hvg_ids, n_sample_cells, n_sample_datasets = compute_hvg(var)

    ens2sym = dict(zip(var["gene_id"], var["gene_symbol"]))
    for e in hvg_ids:
        curated.setdefault(ens2sym[e], []).append("hvg_luca_2000")

    rows = []
    for gene, sets in curated.items():
        ens = sym2ens.get(gene, FALLBACK_ENSEMBL.get(gene, ""))
        rows.append({
            "gene_symbol": gene,
            "ensembl_id": ens,
            "source_set": ",".join(dict.fromkeys(sets)),
            "is_hvg": bool(ens in hvg_ids),
        })
    panel = pd.DataFrame(rows)
    panel["_curated"] = ~panel["source_set"].eq("hvg_luca_2000")
    panel = panel.sort_values(["_curated", "gene_symbol"],
                              ascending=[False, True]).drop(columns="_curated")
    panel.to_csv(f"{OUT}/panel_genes.csv", index=False)

    # cross-source presence audit: can all three arms subset to the same panel?
    audit = presence_audit(panel)
    audit.to_csv(f"{OUT}/panel_genes_source_presence.csv", index=False)

    assert not panel["gene_symbol"].duplicated().any(), "duplicate gene_symbol in panel"
    n_curated = int((panel["source_set"] != "hvg_luca_2000").sum())
    n_hvg_only = int((panel["source_set"] == "hvg_luca_2000").sum())
    both = int(((panel["source_set"] != "hvg_luca_2000") & panel["is_hvg"]).sum())
    print(f"\npanel_genes.csv: {len(panel)} genes "
          f"({n_curated} curated, {n_hvg_only} HVG-only, "
          f"{int(panel['is_hvg'].sum())} flagged is_hvg, "
          f"{both} curated genes that are also HVG)")
    print(f"ensembl_id resolved for {int((panel['ensembl_id'] != '').sum())} "
          f"of {len(panel)}")
    print(f"curated genes absent from LuCA var: {missing if missing else 'none'}")
    print(f"brief-named genes absent from LuCA var: "
          f"{missing_named if missing_named else 'none'}")
    print(f"\ndisambiguated feature_name entries in LuCA var: {len(mangled)}")
    hit = mangled[mangled["gene_symbol"].isin(curated)]
    for _, r in hit.iterrows():
        print(f"  panel gene affected: {r.gene_symbol} stored as {r.feature_name}")

    print("\ncross-source presence (panel_genes_source_presence.csv)")
    for col in ("in_LuCA", "in_EMTAB6149", "in_EMTAB6653"):
        print(f"  {col:14s} {int(audit[col].sum()):5d} / {len(audit)}")
    print(f"  in all three  {int((audit['n_arms'] == 3).sum()):5d} / {len(audit)}")
    cur = audit[audit["source_set"] != "hvg_luca_2000"]
    print(f"  curated genes in all three  {int((cur['n_arms'] == 3).sum())} / {len(cur)}")
    short = cur[cur["n_arms"] < 3].sort_values("gene_symbol")
    for _, r in short.iterrows():
        arms = [n for n, k in (("LuCA", "in_LuCA"), ("6149", "in_EMTAB6149"),
                               ("6653", "in_EMTAB6653")) if not r[k]]
        print(f"    {r.gene_symbol:12s} missing from {', '.join(arms)}  [{r.source_set}]")

    print("\nHVG sample: %d cells, %d LuCA sub-datasets, seurat_v3, batch_key=dataset"
          % (n_sample_cells, n_sample_datasets))
    return panel, audit


if __name__ == "__main__":
    main()
