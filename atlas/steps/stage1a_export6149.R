# =====================================================================================
# Atlas stage 01 (steps export6149) and stage 04 (export6149_full): export the E-MTAB-6149 arm from the
# author RData to something python can read. Part of the atlas behind Figure 2, Figure 3,
# Supplementary Figures 2 to 4 and Supplementary Data 1. No random numbers are drawn.
#
# AllSamplesCelllview.R is a gzip'd RData container (magic RDX2/RDX3), NOT an RDS and
# NOT a script.  It holds three objects:
#     log2cpm      data.frame 22,180 genes x 52,698 cells, rownames = Ensembl gene ids
#     tsne.data    52,698 x 4, dbCluster = the author top-level lineage label
#     featuredata  55,457 x 6 annotation table (Ensembl id -> symbol, biotype, locus)
#
# Despite the object name the values are NOT log2 CPM.  They are ln(1 + c * 1e4 / T),
# i.e. natural-log log1p of a 10,000-scaled library-size normalization.  That is exactly
# invertible: within one cell the smallest positive value is the value of c = 1, so
#     s = exp(min(v[v>0])) - 1  =  1e4 / T      and      c = (exp(v) - 1) / s
# returns the integer UMI count.  This script tests that per cell and reports the worst
# deviation from an integer, rather than asserting it.
#
# The six lineage *.Cellview.Rds files carry the author 52-subtype taxonomy as an integer
# dbCluster inside each lineage.  They are loaded for their cell names and dbCluster only.
#
# Writes to the output directory (argument 1); the deposit is only read.
# =====================================================================================
suppressPackageStartupMessages(library(Matrix))

# Arguments: 1 output directory, 2 row list (or ""), 3 deposit folder, 4 panel table.
#   Rscript stage1a_export6149.R <out dir> <row list or ""> <E-MTAB-6149/processed> <stage1a/panel_genes.csv>
P     <- commandArgs(trailingOnly = TRUE)[3]
PANEL <- commandArgs(trailingOnly = TRUE)[4]
OUT   <- commandArgs(trailingOnly = TRUE)[1]
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
CHUNK <- 2000L

t0 <- Sys.time()
e <- new.env(); load(file.path(P, "AllSamplesCelllview.R"), envir = e)
m  <- get("log2cpm",     envir = e)
td <- get("tsne.data",   envir = e)
fd <- get("featuredata", envir = e)
cat(sprintf("loaded in %.1f s : %d genes x %d cells\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs")), nrow(m), ncol(m)))
stopifnot(identical(colnames(m), rownames(td)))

genes_all <- rownames(m)
cells     <- colnames(m)
n         <- length(cells)

# --- panel rows -----------------------------------------------------------------------
# Second argument, if given, is a one-column file of matrix Ensembl ids to export.  It
# exists because 26 panel genes carry a different Ensembl id in this 2018 featuredata build
# than in the LuCA reference (TXNIP and MRC1 among them); python resolves those by unique
# symbol and hands the resolved matrix ids back here.  Without it a pure id join drops them.
rowfile <- commandArgs(trailingOnly = TRUE)[2]
if (!is.na(rowfile) && nzchar(rowfile)) {
  pid <- readLines(rowfile)
  cat("row list supplied:", length(pid), "ids\n")
} else {
  pid <- read.csv(PANEL, stringsAsFactors = FALSE)$ensembl_id
}
hit <- pid %in% genes_all
cat(sprintf("requested genes %d, present in 6149 matrix %d, absent %d\n",
            length(pid), sum(hit), sum(!hit)))
prow <- match(pid[hit], genes_all)          # row index in m, in requested order
stopifnot(!any(is.na(prow)), !any(duplicated(prow)))

# --- pass over cells in chunks: invert the normalization ------------------------------
umi     <- integer(n); ngene <- integer(n)
maxdev  <- numeric(n); sfac  <- numeric(n)
umi_pan <- integer(n); ngene_pan <- integer(n)
blocks  <- vector("list", ceiling(n / CHUNK))
b <- 0L
for (start in seq(1L, n, by = CHUNK)) {
  idx <- start:min(start + CHUNK - 1L, n)
  M   <- as.matrix(m[, idx, drop = FALSE])          # genes x chunk, double
  pos <- M > 0
  Mi  <- M; Mi[!pos] <- Inf
  s   <- exp(apply(Mi, 2, min)) - 1                 # = 1e4 / total UMI of that cell
  K   <- (exp(M) - 1) / rep(s, each = nrow(M))      # exp(0)-1 == 0, so zeros stay zero
  Kr  <- round(K)
  dev <- abs(K - Kr); dev[!pos] <- 0

  sfac[idx]   <- s
  maxdev[idx] <- apply(dev, 2, max)
  ngene[idx]  <- colSums(pos)
  umi[idx]    <- as.integer(colSums(Kr))

  Kp <- Kr[prow, , drop = FALSE]
  umi_pan[idx]   <- as.integer(colSums(Kp))
  ngene_pan[idx] <- as.integer(colSums(Kp > 0))
  b <- b + 1L
  blocks[[b]] <- as(as(Matrix(Kp, sparse = TRUE), "generalMatrix"), "dgCMatrix")
  if (b %% 5 == 0) cat("  chunk", b, "of", length(blocks), "\n")
}
C <- do.call(cbind, blocks)
dimnames(C) <- list(pid[hit], cells)
cat(sprintf("panel counts: %d x %d, nnz %d\n", nrow(C), ncol(C), length(C@x)))

# --- how exactly do the integers come back? -------------------------------------------
cat(sprintf("INTEGER_RECOVERY tol1e-4 %d / %d ; tol1e-6 %d ; worst deviation %.3e\n",
            sum(maxdev < 1e-4), n, sum(maxdev < 1e-6), max(maxdev)))
cat(sprintf("recovered counts are whole numbers in C: %s ; min %g max %g\n",
            all(C@x == round(C@x)), min(C@x), max(C@x)))

# --- per-cell metadata ----------------------------------------------------------------
lib <- as.integer(sub(".*_", "", cells))
bc  <- sub("_[^_]*$", "", cells)
meta <- data.frame(cell = cells, barcode = bc, library_index = lib,
                   lineage_native = as.character(td$dbCluster),
                   tsne1 = td$V1, tsne2 = td$V2,
                   n_counts = umi, n_genes = ngene,
                   n_counts_panel = umi_pan, n_genes_panel = ngene_pan,
                   implied_library_size = round(1e4 / sfac, 3),
                   max_integer_deviation = signif(maxdev, 4),
                   stringsAsFactors = FALSE)
write.table(meta, file.path(OUT, "cells.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# --- author 52-subtype taxonomy from the six lineage sub-objects ----------------------
subs <- list()
for (f in c("Alveolar", "B_cell", "EC", "Fibro", "Myeloid", "T_cell")) {
  ee <- new.env(); load(file.path(P, paste0(f, ".Cellview.Rds")), envir = ee)
  tt <- get("tsne.data", envir = ee)
  subs[[f]] <- data.frame(cell = rownames(tt), sub_lineage = f,
                          sub_cluster = as.integer(as.character(tt$dbCluster)),
                          stringsAsFactors = FALSE)
  rm(ee); gc(verbose = FALSE)
}
S <- do.call(rbind, subs)
cat(sprintf("sub-object cells %d, unique %d, all in parent %s\n",
            nrow(S), length(unique(S$cell)), all(S$cell %in% cells)))
write.table(S, file.path(OUT, "subclusters.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# --- sex-chromosome probe per library: independent check of library -> patient ---------
gid <- function(s) rownames(fd)[which(fd$Associated.Gene.Name == s)][1]
sx  <- c(XIST = gid("XIST"), RPS4Y1 = gid("RPS4Y1"), DDX3Y = gid("DDX3Y"))
sx  <- sx[!is.na(sx) & sx %in% genes_all]
SX  <- as.matrix(m[sx, , drop = FALSE]); rownames(SX) <- names(sx)
probe <- do.call(rbind, lapply(sort(unique(lib)), function(L) {
  ix <- which(lib == L)
  data.frame(library_index = L, cells = length(ix),
             XIST_det   = round(mean(SX["XIST", ix] > 0), 4),
             RPS4Y1_det = round(mean(SX["RPS4Y1", ix] > 0), 4),
             DDX3Y_det  = round(mean(SX["DDX3Y", ix] > 0), 4))
}))
probe$inferred_sex <- ifelse(probe$XIST_det > probe$RPS4Y1_det, "female", "male")
write.table(probe, file.path(OUT, "sexprobe.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(probe)

# --- matrix out -----------------------------------------------------------------------
writeMM(C, file.path(OUT, "panel_counts.mtx"))
writeLines(rownames(C), file.path(OUT, "panel_genes_present.txt"))
writeLines(colnames(C), file.path(OUT, "panel_cells.txt"))

# featuredata for the genes actually in the matrix, so python can audit symbols
fdm <- data.frame(ensembl_id = genes_all,
                  gene_symbol = as.character(fd[genes_all, "Associated.Gene.Name"]),
                  gene_biotype = as.character(fd[genes_all, "Gene.Biotype"]),
                  chromosome = as.character(fd[genes_all, "Chromosome.Name"]),
                  stringsAsFactors = FALSE)
write.table(fdm, file.path(OUT, "matrix_features.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("DONE", format(Sys.time()), "\n")
