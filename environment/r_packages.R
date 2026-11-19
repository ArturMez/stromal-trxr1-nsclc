# Installs the R packages of mif/02_cells.R (and Matrix, used by the E-MTAB-6149 export of the atlas build) into the
# user library. The published outputs were made with R 4.5.1 and the versions below; a different version prints a
# warning.
#   Rscript environment/r_packages.R
versions <- c(data.table = "1.18.2.1", readxl = "1.5.0.1", autothresholdr = "1.4.3", jsonlite = "2.0.0",
              Matrix = "1.7.3")
new <- names(versions)[!sapply(names(versions), requireNamespace, quietly = TRUE)]
if (length(new)) install.packages(new, repos = "https://cloud.r-project.org")
for (p in names(versions)) {
  v <- as.character(packageVersion(p))
  if (v != versions[[p]]) warning(sprintf("%s %s is installed; the outputs were made with %s", p, v, versions[[p]]))
}
