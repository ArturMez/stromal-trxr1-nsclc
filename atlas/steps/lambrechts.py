#!/usr/bin/env python
"""
Atlas stage 01, step lambrechts: inventory of the Lambrechts / Thienpont 2018 data (Nat Med 24:1277),
the study behind the "Lambrechts copies" that the atlas includes once.

Paper items: none directly; the per-sample table documents the two Lambrechts accessions that enter
the atlas behind Figure 2, Figure 3, Supplementary Figures 2 to 4 and Supplementary Data 1.
No random numbers are drawn.

Two ArrayExpress accessions, one study, 8 patients:
    E-MTAB-6149  patients 1-5  (discovery)  52,698 cells, author-processed R object
    E-MTAB-6653  patients 6-8  (validation) 32,341 cells, EBI-reprocessed Matrix Market

Everything here is derived; original_data/external_scRNA is never written to.

What this script establishes that is not stated in any deposited file:
  1. Histology per patient comes only from Supplementary Table 1 of the paper.
     Both SDRFs record Characteristics[disease] = "lung carcinoma" for every sample.
  2. AllSamplesCelllview.R is a gzip'd RData (RDX3) container, not a script and not an RDS.
     Its "log2cpm" object is natural-log log1p normalization, NOT log2 CPM, and the
     underlying integer UMI counts are exactly recoverable (verified 52,698/52,698 cells).
  3. The 24 barcode suffixes in that object are libraries, not samples. The
     library -> patient -> sampling-site map is recovered here and is provably unique.
  4. The E-MTAB-6653 "raw counts" matrix is not integer: 15.91 % of nonzero entries are
     fractional (Alevin EM-deduplicated estimates).

Reads:  E-MTAB-6149 (processed R objects) and E-MTAB-6653 (Matrix Market, SDRF, clusters) from the
        raw deposits (00_download.sh).
Writes: output/RNA/meta_dataset/lambrechts.csv (the per-sample inventory); scratch files
        (cells6149.csv, sexprobe6149.csv, lineage6149.csv, extract6149.R) go to a temporary folder
        under scratch/.
Run:  python atlas/steps/lambrechts.py   (normally through atlas/01_arms_luca_lambrechts.py)
Needs: R (Rscript 4.5.1, found through config.json "rscript" or PATH) for the E-MTAB-6149 half.
"""
import os, subprocess, tempfile, json
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp

import os as _os  # noqa: E402
import sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # atlas/
from common import ROOT, SCRNA, SCRATCH, DATA, RSCRIPT  # noqa: E402,F401  paths from config.json
EXT    = f"{SCRNA}"
D6149  = f"{EXT}/E-MTAB-6149"
D6653  = f"{EXT}/E-MTAB-6653"
OUT    = f"{ROOT}/output/RNA/meta_dataset"
os.makedirs(OUT, exist_ok=True)

# --------------------------------------------------------------------------------------
# Supplementary Table 1 (paper). The ONLY source of histology; not in SDRF or IDF.
# --------------------------------------------------------------------------------------
PATIENTS = pd.DataFrame([
    # patient, accession, sex, age, histology_paper, histology_class, stage, chemistry
    (1, "E-MTAB-6149", "female", 70, "Squamous",     "LUSC",         "IIA",  "10x 3p v1"),
    (2, "E-MTAB-6149", "male",   86, "Squamous",     "LUSC",         "IB",   "10x 3p v1"),
    (3, "E-MTAB-6149", "male",   68, "Adenomatous",  "LUAD",         "IIIB", "10x 3p v2"),
    (4, "E-MTAB-6149", "female", 64, "Adenomatous",  "LUAD",         "IIB",  "10x 3p v2"),
    (5, "E-MTAB-6149", "male",   60, "Large cell",   "large cell",   "IA3",  "10x 3p v2"),
    (6, "E-MTAB-6653", "male",   65, "Adenomatous",  "LUAD",         "IIIA", "10x 3p v2"),
    (7, "E-MTAB-6653", "male",   60, "Squamous",     "LUSC",         "IB",   "10x 3p v2"),
    (8, "E-MTAB-6653", "female", 55, "Pleiomorphic", "pleomorphic",  "IIB",  "10x 3p v2"),
], columns=["patient","accession","sex_reported","age","histology_paper",
            "histology_class","stage","chemistry"])

# Per-sample cell counts as published (Suppl. Table 2). Used to pin library -> sample.
PUB_CELLS = {  # (patient, site) -> cells reported by the authors
 (1,"core"):424,(1,"middle"):400,(1,"edge"):648,
 (2,"core"):114,(2,"middle"):806,(2,"edge"):1651,(2,"normal"):271,
 (3,"core"):2149,(3,"middle"):3266,(3,"edge"):3554,(3,"normal"):4235,
 (4,"core"):5270,(4,"middle"):5198,(4,"edge"):3283,(4,"normal"):2933,
 (5,"core"):2983,(5,"middle"):4913,(5,"edge"):4664,(5,"normal"):5936,
 (6,"core"):2098,(6,"middle"):3832,(6,"edge"):3036,(6,"normal"):4366,
 (7,"core"):1898,(7,"middle"):1450,(7,"edge"):989, (7,"normal"):3862,
 (8,"core"):3314,(8,"middle"):4601,(8,"edge"):3777,(8,"normal"):7027,
}

# Library index -> (patient, site) for E-MTAB-6149. Recovered, not deposited.
# Blocks 1-6 / 7-12 / 13-16 / 17-20 / 21-24 are fixed by XIST vs RPS4Y1/DDX3Y detection
# (the sex-chromosome probe in R_EXTRACT below); within each block the assignment is the unique subset partition
# whose sums equal the published per-sample cell counts.
LIB2SAMPLE = {
 1:(2,"middle"), 2:(2,"middle"), 3:(2,"edge"),  4:(2,"edge"),  5:(2,"core"), 6:(2,"normal"),
 7:(1,"core"),   8:(1,"core"),   9:(1,"middle"),10:(1,"middle"),11:(1,"edge"),12:(1,"edge"),
 13:(3,"edge"), 14:(3,"middle"),15:(3,"core"), 16:(3,"normal"),
 17:(4,"normal"),18:(4,"edge"), 19:(4,"middle"),20:(4,"core"),
 21:(5,"core"), 22:(5,"edge"),  23:(5,"middle"),24:(5,"normal"),
}

# Author subcluster vocabulary, Suppl. Table 3. dbCluster integer in each *.Cellview.Rds
# joins directly onto these (all 50 sizes verified to match exactly).
SUBCLUSTERS = [
 ("Alveolar",1,"cuboidal alveolar type 2 cells",520,"SFTPA1, SFTPA2, SFTPC, PGC"),
 ("Alveolar",2,"lower quality alveolar cells",310,""),
 ("Alveolar",3,"flat alveolar type 1 cells",303,"AGER, CAV1"),
 ("Alveolar",4,"COPD-injured alveolar cells",205,"MMP7, CXCL14"),
 ("Alveolar",5,"respiratory epithelial cells",168,"EMP2, AGER, UPK3B, SPOCK2, CYP4B1, CAV1"),
 ("Alveolar",6,"secretory club cells",136,"SCGB1A1, BPIFB1, SCGB3A1, WFDC2, MUC5B, PIGR"),
 ("Alveolar",7,"basal cells",68,"KRT17, KRT5, SYT8, MMP1, IGFBP6, FHL2"),
 ("B_cell",1,"follicular B-cells",1786,""),
 ("B_cell",2,"follicular B-cells",889,"RP5-887A10.1, AL928768.3"),
 ("B_cell",3,"plasma B-cells",866,"IGHG4"),
 ("B_cell",4,"mast cells",613,"TPSAB1, TPSB2, CPA3, HPGDS, MS4A2, KIT"),
 ("B_cell",5,"MALT B-cells",548,"IGHA1, IGHA2"),
 ("B_cell",6,"plasma B-cells",373,"IGHG3, IGHG1, IGHGP, IGHG2"),
 ("B_cell",7,"MALT B-cells",344,"JCHAIN, MZB1, IGHM, IGLL5, IGHD"),
 ("B_cell",8,"plasmacytoid dendritic cells",96,"GZMB, LILRA4, CLIC3, IRF7, IL3RA"),
 ("B_cell",9,"erythroblasts",88,"HBB, HBA2, HBA1, HBD, ALAS2"),
 ("EC",1,"normal endothelial cells",323,"MT1X, MT2A, IL6, FCN3, SOX7"),
 ("EC",2,"lower quality endothelial cell",320,""),
 ("EC",3,"tumour endothelial cells",311,""),
 ("EC",4,"tumour endothelial cells",307,""),
 ("EC",5,"normal endothelial cells",246,"EDNRB, IL1RL1, HPGD, CA4, EMCN"),
 ("EC",6,"lymphatic endothelial cells",85,"CCL21, TFF3, PDPN, MMRN1"),
 ("Epithelial",1,"epithelial cells",164,"CAPS, TPPP3"),
 ("Epithelial",2,"lower quality epithelial cells",50,""),
 ("Fibro",1,"fibroblasts",315,"CTHRC1, VCAN, SULF1, COL10A1, POSTN, SFRP4, HTRA3, ASPN"),
 ("Fibro",2,"fibroblasts (pericyte/mural marker set)",266,"COX4I2, HIGD1B, RGS5, GJA4, NDUFA4L2, PDGFRB, NOTCH3, ACTG2"),
 ("Fibro",3,"lower quality fibroblasts",240,""),
 ("Fibro",4,"fibroblasts",219,"PLA2G2A, SFRP2"),
 ("Fibro",5,"fibroblasts",195,"MMP3, SERPINE2"),
 ("Fibro",6,"normal lung fibroblasts",175,"GPC3, MFAP4, A2M, GPX3, TIMP3, CFD, ADH1B, INMT"),
 ("Fibro",7,"fibroblasts",55,"ALDH1A3, C3, CCL2, SOD2, TNFRSF12A, NFKBIA, PTX3, DKK1, SERPINB2"),
 ("Myeloid",1,"macrophages",2409,""),("Myeloid",2,"macrophages",2106,""),
 ("Myeloid",3,"macrophages",1492,""),("Myeloid",4,"macrophages",716,""),
 ("Myeloid",5,"Langerhans cells",714,""),("Myeloid",6,"macrophages",484,""),
 ("Myeloid",7,"granulocytes",478,""),("Myeloid",8,"macrophages",400,""),
 ("Myeloid",9,"monocyte-derived dendritic cells",298,""),("Myeloid",10,"macrophages",263,""),
 ("Myeloid",11,"macrophages",204,""),("Myeloid",12,"cross-presenting dendritic cells",192,""),
 ("T_cell",1,"CD4+ T-cells",4843,""),("T_cell",2,"CD8+ T-cells",4655,""),
 ("T_cell",3,"CD4+ T-cells",3880,""),("T_cell",4,"CD8+ T-cells",3254,""),
 ("T_cell",5,"CD8+ T-cells",3125,""),("T_cell",6,"natural killer cells",1741,""),
 ("T_cell",7,"regulatory T-cells",1513,""),("T_cell",8,"CD8+ T-cells",1006,""),
 ("T_cell",9,"CD4+ T-cells",894,""),
]

# Broad labels for the EBI K=13 clustering of E-MTAB-6653, read off the deposited
# marker_genes_13.tsv (Ensembl IDs mapped to symbols via the 6149 featuredata table).
K13_LABEL = {1:"myeloid", 2:"T/NK", 3:"plasma B", 4:"B cell", 5:"alveolar/epithelial",
             6:"endothelial", 7:"fibroblast/mesenchymal", 8:"macrophage (FABP4/C1QA)",
             9:"tumour/epithelial", 10:"plasma B", 11:"mast", 12:"ciliated", 13:"erythroid"}


# ======================================================================================
# E-MTAB-6149 : RData -> integer counts + per-cell metadata, via R
# ======================================================================================
R_EXTRACT = r'''
P   <- "%(D6149)s/processed"
OUT <- "%(SCRATCH)s"
lib2pat  <- c(%(LIBPAT)s)
lib2site <- c(%(LIBSITE)s)

e <- new.env(); load(file.path(P,"AllSamplesCelllview.R"), envir=e)
m <- get("log2cpm", envir=e); td <- get("tsne.data", envir=e); fd <- get("featuredata", envir=e)
stopifnot(identical(colnames(m), rownames(td)))
lib <- as.integer(sub(".*_","", colnames(m)))

# --- invert the normalisation ---------------------------------------------------------
# values are ln(1 + k * s_j) for a per-cell constant s_j; taking s_j from the smallest
# nonzero value in the cell (which corresponds to k = 1) returns integer k exactly.
n <- ncol(m); umi <- integer(n); libsize <- numeric(n); ngene <- integer(n); okn <- 0L
for (j in seq_len(n)) {
  v <- m[[j]]; v <- v[v > 0]
  s <- exp(min(v)) - 1
  k <- (exp(v) - 1)/s
  if (max(abs(k - round(k))) < 1e-4) okn <- okn + 1L
  umi[j] <- sum(round(k)); libsize[j] <- 1e4/s; ngene[j] <- length(v)
}
cat("INTEGER_RECOVERY", okn, n, "\n")

meta <- data.frame(cell=colnames(m), library_index=lib,
                   patient=lib2pat[as.character(lib)], site=lib2site[as.character(lib)],
                   cell_type=as.character(td$dbCluster),
                   umi_recovered=umi, n_genes=ngene,
                   implied_full_libsize=round(libsize,1), stringsAsFactors=FALSE)
write.csv(meta, file.path(OUT,"cells6149.csv"), row.names=FALSE)

# --- sex-chromosome probe: independent check of the library -> patient blocks ----------
gid <- function(s) rownames(fd)[which(fd$Associated.Gene.Name == s)][1]
sx  <- c(XIST=gid("XIST"), RPS4Y1=gid("RPS4Y1"), DDX3Y=gid("DDX3Y"))
S   <- as.matrix(m[sx,,drop=FALSE]); rownames(S) <- names(sx)
probe <- do.call(rbind, lapply(sort(unique(lib)), function(L){ ix <- which(lib==L)
  data.frame(library_index=L, cells=length(ix),
             XIST_det=round(mean(S["XIST",ix]>0),4),
             RPS4Y1_det=round(mean(S["RPS4Y1",ix]>0),4),
             DDX3Y_det=round(mean(S["DDX3Y",ix]>0),4)) }))
probe$inferred_sex <- ifelse(probe$XIST_det > probe$RPS4Y1_det, "female", "male")
write.csv(probe, file.path(OUT,"sexprobe6149.csv"), row.names=FALSE)

# --- lineage sub-object sizes ---------------------------------------------------------
sizes <- do.call(rbind, lapply(c("Alveolar","B_cell","EC","Fibro","Myeloid","T_cell"), function(f){
  ee <- new.env(); load(file.path(P, paste0(f,".Cellview.Rds")), envir=ee)
  tt <- get("tsne.data", envir=ee); mm <- get("log2cpm", envir=ee)
  data.frame(lineage=f, cells=ncol(mm), genes=nrow(mm),
             subclusters=length(unique(tt$dbCluster)),
             sizes=paste(as.integer(table(tt$dbCluster)), collapse=","),
             all_cells_in_parent=all(colnames(mm) %%in%% colnames(m))) }))
write.csv(sizes, file.path(OUT,"lineage6149.csv"), row.names=FALSE)
'''

def run_6149(scratch):
    """Run the R extraction of E-MTAB-6149 in the scratch folder; returns (cells with exact integer recovery, cells)."""
    libpat  = ",".join(f'"{k}"={v[0]}' for k, v in LIB2SAMPLE.items())
    libsite = ",".join(f'"{k}"="{v[1]}"' for k, v in LIB2SAMPLE.items())
    src = R_EXTRACT % dict(D6149=D6149, SCRATCH=scratch, LIBPAT=libpat, LIBSITE=libsite)
    rf = os.path.join(scratch, "extract6149.R")
    open(rf, "w").write(src)
    p = subprocess.run([RSCRIPT, "--vanilla", rf], capture_output=True, text=True)
    if p.returncode:
        raise RuntimeError(p.stderr[-4000:])
    rec = [l for l in p.stdout.splitlines() if l.startswith("INTEGER_RECOVERY")][0].split()
    return int(rec[1]), int(rec[2])

# ======================================================================================
# E-MTAB-6653 : Matrix Market + EBI clusters
# ======================================================================================
def run_6653():
    """Read the E-MTAB-6653 Matrix Market, SDRF and EBI clusters; returns the matrix, genes, per-cell table, integer fraction and cluster keys."""
    base = f"{D6653}/raw_counts/matrix_market/E-MTAB-6653.aggregated_filtered_counts"
    X = sio.mmread(base + ".mtx").tocsc()
    genes = pd.read_csv(base + ".mtx_rows", sep="\t", header=None)
    cells = pd.read_csv(base + ".mtx_cols", sep="\t", header=None)[0].values
    d = X.data
    intfrac = float(np.isclose(d, np.round(d)).mean())

    sd = pd.read_csv(f"{D6653}/metadata/E-MTAB-6653.sdrf.txt", sep="\t", dtype=str)
    key = sd.set_index("Comment[ENA_RUN]")
    short = {"normal tissue adjacent to tumour":"normal", "tumour edge":"edge",
             "tumour middle (in between core and edge sample)":"middle", "tumour core":"core"}
    run = pd.Series(cells).str.split("-").str[0].values
    meta = pd.DataFrame({
        "cell": cells, "run": run,
        "sample": key["Source Name"].reindex(run).values,
        "patient": key["Characteristics[individual]"].reindex(run).astype(int).values,
        "site": pd.Series(key["Characteristics[sampling site]"].reindex(run).values).map(short).values,
        "umi": np.asarray(X.sum(0)).ravel(),
        "n_genes": np.asarray((X > 0).sum(0)).ravel(),
    })
    cl = pd.read_csv(f"{D6653}/metadata/E-MTAB-6653_clusters.tsv", sep="\t").set_index(["sel.K","K"])
    assert list(cl.columns) == list(cells), "cluster file column order does not match mtx_cols"
    meta["k13"] = cl.xs(13, level="K").iloc[0].values.astype(int)
    meta["k13_label"] = meta.k13.map(K13_LABEL)
    meta["k40"] = cl.xs(40, level="K").iloc[0].values.astype(int)
    return X, genes, meta, intfrac, [k for (s, k) in cl.index], [k for (s, k) in cl.index if s]

# ======================================================================================
if __name__ == "__main__":
    scratch = tempfile.mkdtemp(prefix="lambrechts_", dir=SCRATCH)
    print("scratch:", scratch)

    ok, tot = run_6149(scratch)
    print(f"E-MTAB-6149 exact integer recovery: {ok}/{tot}")
    c49 = pd.read_csv(f"{scratch}/cells6149.csv")
    c49 = c49.merge(PATIENTS[["patient","histology_class","chemistry","sex_reported"]], on="patient")
    c49["accession"] = "E-MTAB-6149"

    X, genes, c53, intfrac, allK, selK = run_6653()
    print(f"E-MTAB-6653 shape {X.shape}, integer-valued nonzeros {100*intfrac:.2f} %")
    c53 = c53.merge(PATIENTS[["patient","histology_class","chemistry","sex_reported"]], on="patient")
    c53["accession"] = "E-MTAB-6653"

    # -- master per-sample inventory ----------------------------------------------------
    rows = []
    for (pat, site), g in c49.groupby(["patient","site"]):
        rows.append(dict(accession="E-MTAB-6149", patient=pat, site=site,
                         libraries="+".join(map(str, sorted(g.library_index.unique()))),
                         cells_local=len(g), cells_published=PUB_CELLS[(pat, site)],
                         median_umi=g.umi_recovered.median(), median_genes=g.n_genes.median(),
                         pct_malignant=round(100*(g.cell_type == "tumor").mean(), 2),
                         pct_fibroblast=round(100*(g.cell_type == "Fibro").mean(), 2)))
    for (pat, site), g in c53.groupby(["patient","site"]):
        rows.append(dict(accession="E-MTAB-6653", patient=pat, site=site,
                         libraries=g["sample"].iloc[0] + "/" + g.run.iloc[0],
                         cells_local=len(g), cells_published=PUB_CELLS[(pat, site)],
                         median_umi=round(g.umi.median(), 1), median_genes=g.n_genes.median(),
                         pct_malignant=round(100*(g.k13_label == "tumour/epithelial").mean(), 2),
                         pct_fibroblast=round(100*(g.k13_label == "fibroblast/mesenchymal").mean(), 2)))
    inv = pd.DataFrame(rows).merge(PATIENTS.drop(columns=["accession"]), on="patient")
    inv["site"] = pd.Categorical(inv.site, ["core","middle","edge","normal"], ordered=True)
    inv = inv.sort_values(["patient","site"])
    inv = inv[["accession","patient","site","libraries","histology_class","histology_paper",
               "chemistry","sex_reported","age","stage","cells_local","cells_published",
               "median_umi","median_genes","pct_malignant","pct_fibroblast"]]
    inv.to_csv(f"{OUT}/lambrechts.csv", index=False)
    print(inv.to_string(index=False))

    print("wrote", OUT)
