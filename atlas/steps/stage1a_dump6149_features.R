# Atlas stage 01, step dump6149_features: dump the E-MTAB-6149 featuredata (Ensembl id, symbol, in
# matrix) for the panel step (stage1a_panel.py).
# Paper items: part of the atlas behind Figure 2, Figure 3, Supplementary Figures 2 to 4 and
# Supplementary Data 1. No random numbers are drawn.
# Args: <E-MTAB-6149/processed folder> <output csv, stage1a/_features6149.csv>
# Run (normally through atlas/01_arms_luca_lambrechts.py):
#   Rscript stage1a_dump6149_features.R <E-MTAB-6149/processed> <output csv>
args <- commandArgs(trailingOnly = TRUE)
P <- args[1]
OUTCSV <- args[2]
dir.create(dirname(OUTCSV), recursive = TRUE, showWarnings = FALSE)
e <- new.env()
load(file.path(P, "AllSamplesCelllview.R"), envir = e)
cat("objects:", paste(ls(e), collapse=", "), "\n")
fd <- get("featuredata", envir = e)
m  <- get("log2cpm", envir = e)
cat("featuredata cols:", paste(colnames(fd), collapse=", "), "\n")
cat("featuredata rows:", nrow(fd), " matrix rows:", nrow(m), "\n")
out <- data.frame(ensembl_id = rownames(fd),
                  gene_symbol = as.character(fd[["Associated.Gene.Name"]]),
                  in_matrix = rownames(fd) %in% rownames(m),
                  stringsAsFactors = FALSE)
write.csv(out, OUTCSV, row.names = FALSE)
cat("written\n")
