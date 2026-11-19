# Atlas stage 04, step dump6149_rows: dump every Ensembl id in the E-MTAB-6149 matrix, so the exporter
# stage1a_export6149.R can be reused with a FULL gene list instead of the 1,906-gene panel.
# The full-gene build needs all measured genes; the exporter takes a row list as argument 2.
# Paper items: part of the atlas behind Figure 2, Figure 3, Supplementary Figures 2 to 4 and
# Supplementary Data 1. No random numbers are drawn.
# Output: <out dir>/all_gene_ids.txt (one Ensembl id per line).
# Run (normally through atlas/04_build_objects.py):
#   Rscript fullbuild_01d_dump6149_rows.R <out dir> <E-MTAB-6149/processed>
OUT <- commandArgs(trailingOnly = TRUE)[1]
P <- commandArgs(trailingOnly = TRUE)[2]
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
e <- new.env(); load(file.path(P, "AllSamplesCelllview.R"), envir = e)
m <- get("log2cpm", envir = e)
cat(sprintf("matrix %d genes x %d cells\n", nrow(m), ncol(m)))
writeLines(rownames(m), file.path(OUT, "all_gene_ids.txt"))
cat("wrote", nrow(m), "gene ids\n")
