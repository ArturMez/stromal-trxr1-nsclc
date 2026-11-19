#!/usr/bin/env Rscript
################################################################################
##  The per-cell table of a cohort, from the raw inForm export.
##
##  WHAT: reads the raw single-cell export, derives one positivity cut per marker, writes the
##        binary calls and the lineage of every Tumor or Stroma cell of a primary core, and
##        converts the per-core tissue areas to square millimeters.
##  FROM: raw_root/Raw_original_<cohort>/{cell_seg_data, tissue_seg_data_summary} and the core key
##        files under raw_root/Bomi2 and raw_root/Bomi1.
##  TO:   work_root/intermediate/<cohort>/{thresholds.csv, cells.tsv, areas.csv}, plus a sibling .md
##        for each that documents the columns.
##  ITEMS: no figure or table of the paper is made here; every later script starts from these
##        tables (cells.tsv: lineage and TrxR1 call of each cell; areas.csv: segmented areas).
##  WHY:  one script does both cohorts, with everything that differs between them in one block,
##        saves the cuts it used (thresholds.csv), converts the areas with the same micrometer
##        constant as the positions, and derives the reclassification of one patient from the
##        annotation sheet.
##
##  The step stays in R because the two cuts that are not fixed constants are defined by R's own
##  density estimator and by autothresholdr's Otsu implementation. Rewriting them in another
##  language would move the cut, and the cut decides the positivity of three million cells.
##
##  RUN: Rscript mif/02_cells.R --cohort BOMI2|BOMI1 [--phase-a]
##       The configuration comes from $REDOX_CONFIG, or from config.json beside the repository
##       root. --phase-a (BOMI2, Uppsala II, only) documents how the three fixed Uppsala II lineage
##       cuts (CD3, CD68 and alpha-SMA) were derived: it reruns the review of the hand-annotated
##       cores and writes the suggested cuts to thresholds_phase_a.csv; it changes nothing else.
##
##  Channels not analyzed in the study are not read. A cell positive for more than one lineage
##  marker is assigned to the first lineage in the priority order (CD3, CD68, cytokeratin,
##  alpha-SMA in Uppsala II; CD3, cytokeratin in Uppsala I).
##
##  Deps: data.table, jsonlite, readxl, autothresholdr.
################################################################################

suppressPackageStartupMessages({
  library(data.table); library(jsonlite); library(readxl); library(autothresholdr)
})

## ---------------------------------------------------------------- configuration and arguments

script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (!length(f)) stop("run this with Rscript, so that the repository root can be found")
  dirname(normalizePath(f[1]))
}
REPO <- dirname(script_dir())

read_config <- function() {
  named <- Sys.getenv("REDOX_CONFIG", unset = "")
  path <- if (nzchar(named)) named else file.path(REPO, "config.json")
  if (!file.exists(path))
    stop(sprintf("no configuration at %s. Copy config.example.json to config.json and set the roots.",
                 path))
  cfg <- jsonlite::fromJSON(path)
  for (k in c("raw_root", "work_root"))
    if (is.null(cfg[[k]]) || !nzchar(cfg[[k]])) stop(sprintf("%s has no %s", path, k))
  cfg
}

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default = "") {
  i <- match(name, args)
  if (is.na(i)) default else args[i + 1L]
}
COHORT  <- arg_value("--cohort", "")
PHASE_A <- "--phase-a" %in% args
if (!COHORT %in% c("BOMI2", "BOMI1"))
  stop("pass --cohort BOMI2 or --cohort BOMI1")
if (PHASE_A && COHORT != "BOMI2")
  stop("the annotated-core review exists for BOMI2 only")

CFG <- read_config()
raw_path <- function(...) file.path(CFG$raw_root, ...)
## the file name of a raw input that config.json names under "input_files"
input_file <- function(key) {
  name <- CFG$input_files[[key]]
  if (is.null(name) || !nzchar(name)) stop(sprintf("config.json has no input_files entry %s", key))
  name
}
out_dir  <- file.path(CFG$work_root, "intermediate", COHORT)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_path <- function(stem, ext) file.path(out_dir, paste0(stem, ".", ext))

## ------------------------------------------------------------------------ cohort configuration
##
##  Everything that differs between the two cohorts sits here and nowhere else. The panels differ:
##  BOMI1 was stained with a reduced panel whose Opal 620 and Opal 650 channels are empty, so that
##  cohort has no macrophage and no fibroblast marker, and its "Other" lineage holds the cells that
##  BOMI2 separates into those two. Channels not analyzed in the study are not read.

CHANNEL <- function(opal) sprintf("Cytoplasm Opal %d Mean (Normalized Counts, Total Weighting)", opal)

COHORTS <- list(
  BOMI2 = list(
    cells  = raw_path("Raw_original_BOMI2", "BOMI2_REDOX_merged_cell_seg_data.txt"),
    tissue = raw_path("Raw_original_BOMI2", "BOMI2_REDOX_merged_tissue_seg_data_summary.txt"),
    keys   = raw_path("Bomi2"),
    ## marker -> inForm cytoplasm channel, in lineage priority order
    channels = c(CD3 = CHANNEL(540), CD68 = CHANNEL(620), ASMA = CHANNEL(650), CK = CHANNEL(690)),
    trxr1    = CHANNEL(570),
    lineage  = c(CD3 = "T", CD68 = "Mac", CK = "CK", ASMA = "Fib"),
    ## the tissue summary gives pixels; one constant for the whole package (redox.coords.MPP)
    area_column = "Region Area (pixels)",
    mpp = 0.4964),
  BOMI1 = list(
    cells  = raw_path("Raw_original_BOMI1", input_file("bomi1_cell_export")),
    tissue = raw_path("Raw_original_BOMI1", input_file("bomi1_tissue_export")),
    keys   = raw_path("Bomi1"),
    channels = c(CD3 = CHANNEL(540), CK = CHANNEL(690)),
    trxr1    = CHANNEL(570),
    lineage  = c(CD3 = "T", CK = "CK"),
    ## the BOMI1 tissue summary gives square micrometers and the positions are micrometers already
    area_column = "Region Area (square microns)",
    mpp = 1.0))
CO <- COHORTS[[COHORT]]

## The three BOMI2 lineage cuts are fixed constants, not derived at run time. They were chosen in
## the Phase A review of 22 annotated cores (--phase-a below reruns it): the density crossover of
## annotated positive against annotated negative cells fell at 4.340 for CD3, 4.512 for CD68 and
## 9.243 for alpha-SMA, and the reviewed cuts were rounded up from those three numbers. They are
## absolute intensities on the BOMI2 scale and must not be carried to any other cohort, which is
## why BOMI1 derives its own CD3 cut from its own intensities.
BOMI2_FIXED_CUTS <- c(CD3 = 4.4, CD68 = 4.5, ASMA = 9.3)

## The BOMI1 CD3 cut is the 95.2nd percentile of CD3 inside cytokeratin-positive cells. Those cells
## are the cohort's own CD3-negative anchor, and the percentile transfers the BOMI2 operating point
## (about 95 % specificity) onto the BOMI1 intensity scale, which sits lower than the BOMI2 one.
BOMI1_CD3_PERCENTILE <- 0.952

## Cuts of the reference run of this script (the values the manuscript's tables were built with).
## The run stops when a reproduced cut moves from them by more than the relative tolerance below,
## because a moved cut changes the positivity of millions of cells silently.
EXPECTED <- list(BOMI2 = c(CK = 1.743312, TrxR1 = 4.651484),
                 BOMI1 = c(CK = 0.397389, CD3 = 5.079, TrxR1 = 2.776411))
EXPECTED_TOLERANCE <- 1e-3   # relative

## -------------------------------------------------------------------------- threshold methods

## Cytokeratin: an Otsu seed, then the lowest point of the intensity density between the two modes
## it separates. Otsu alone sits too high on a bimodal distribution with a heavy positive tail.
ck_cut <- function(x) {
  x  <- x[is.finite(x)]
  o  <- as.numeric(auto_thresh(as.integer(x * 1000), "Otsu")) / 1000
  p1 <- density(x[x < o]); p1 <- p1$x[which.max(p1$y)]
  p2 <- density(x[x > o]); p2 <- p2$x[which.max(p2$y)]
  dd <- density(x); reg <- dd$x > p1 & dd$x < p2
  dd$x[which(dd$y == min(dd$y[reg]))][1]
}

## TrxR1: the crossover of the Tumor and the Stroma intensity density, taken inside the window
## between the two medians. The two compartments carry different amounts of the protein, so the
## point where a cell stops looking like a stromal cell and starts looking like a tumor cell is a
## cut that does not depend on any annotation. The top percent of intensities is trimmed first
## because a handful of very bright cells otherwise sets the range of the density estimate.
## Fallbacks: when the window holds no crossing the whole range is searched, and when that fails
## too the midpoint of the medians is used.
crossover_cut <- function(a, b) {
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  if (length(a) < 5 || length(b) < 5) return(NA_real_)
  hi <- as.numeric(quantile(c(a, b), 0.99)); a <- a[a <= hi]; b <- b[b <= hi]
  da <- density(a, from = 0, to = hi, n = 1024); db <- density(b, from = 0, to = hi, n = 1024)
  x  <- da$x; diffy <- da$y - db$y
  pick <- function(lo, up) {
    reg <- which(x >= lo & x <= up)
    s <- which(diff(sign(diffy[reg])) != 0)
    if (length(s)) x[reg][s[1]] else NA_real_
  }
  m1 <- min(median(a), median(b)); m2 <- max(median(a), median(b))
  cut <- pick(m1, m2)
  if (is.na(cut)) cut <- pick(0, hi)
  if (is.na(cut)) cut <- (m1 + m2) / 2
  cut
}

## ------------------------------------------------------------- Phase A, the annotated-core review
##
##  Optional. It reproduces the review that set the three fixed BOMI2 lineage cuts, from the 22
##  cores that were phenotyped by hand: each core was exported twice, once with the lineage call
##  and once plain, and the two exports are joined on the cell. The suggestion for a marker is the
##  crossover of the density of the annotated positive cells with the density of the annotated
##  negative ones, inside the window between the two medians. It writes a table and changes no
##  other output, so that the reviewed constants above stay the cuts the package applies.

phase_a <- function() {
  folder <- raw_path("Raw_original_BOMI2", "outputPheno")
  if (!dir.exists(folder)) stop(sprintf("no annotated-core folder at %s", folder))
  keep <- c("Sample Name", "Cell ID", unname(CO$channels))
  files <- list.files(folder, pattern = "cell_seg_data\\.txt$", full.names = TRUE)
  base <- sub("(_[0-9]+)?_cell_seg_data\\.txt$", "", basename(files))
  cells <- rbindlist(lapply(split(files, base), function(fs) {
    numbered <- grepl("_[0-9]+_cell_seg_data\\.txt$", basename(fs))
    plain <- fs[!numbered]; lineage <- fs[numbered]
    if (!length(plain) || !length(lineage)) return(data.table())
    P <- fread(plain[1], sep = "\t", showProgress = FALSE)
    P <- P[, intersect(keep, names(P)), with = FALSE]
    setnames(P, unname(CO$channels), names(CO$channels), skip_absent = TRUE)
    N <- fread(lineage[1], sep = "\t", showProgress = FALSE)[, .(`Cell ID`, call = Phenotype)]
    merge(P, N, by = "Cell ID", all = FALSE)
  }), use.names = TRUE, fill = TRUE)
  cat(sprintf("[A] %d annotated cells across %d cores\n", nrow(cells), uniqueN(cells$`Sample Name`)))

  annotated <- c(CD3 = "CD3", CD68 = "CD68", ASMA = "ASMA", CK = "PanCK")
  negative <- cells[call == "negative"]
  suggest <- function(pos, neg) {
    pos <- pos[is.finite(pos)]; neg <- neg[is.finite(neg)]
    if (length(pos) < 5) return(NA_real_)
    rng <- range(c(pos, neg))
    dp <- density(pos, from = rng[1], to = rng[2], n = 512)
    dn <- density(neg, from = rng[1], to = rng[2], n = 512)
    reg <- which(dp$x >= min(median(neg), median(pos)) & dp$x <= max(median(neg), median(pos)))
    if (length(reg) < 2) return(NA_real_)
    s <- which(diff(sign(dp$y[reg] - dn$y[reg])) != 0)
    if (length(s)) dp$x[reg][s[1]] else NA_real_
  }
  rows <- rbindlist(lapply(names(annotated), function(mk) {
    pos <- cells[call == annotated[[mk]]][[mk]]
    data.table(marker = mk, suggested_cut = suggest(pos, negative[[mk]]),
               n_positive = sum(is.finite(pos)), n_negative = sum(is.finite(negative[[mk]])),
               median_positive = median(pos, na.rm = TRUE),
               median_negative = median(negative[[mk]], na.rm = TRUE))
  }))
  rows[, applied_cut := as.numeric(BOMI2_FIXED_CUTS[marker])]
  print(rows)
  target <- out_path("thresholds_phase_a", "csv")
  fwrite(rows, target)
  writeLines(sprintf(
"# thresholds_phase_a.csv (BOMI2)

The annotated-core review that set the three fixed BOMI2 lineage cuts, rerun by
Rscript mif/02_cells.R --cohort BOMI2 --phase-a. %d cells of %d hand-phenotyped cores.

| column | meaning |
|---|---|
| marker | the marker the review looked at |
| suggested_cut | crossover of the density of the annotated positive cells with the density of the annotated negative ones, inside the window between the two medians |
| n_positive | annotated positive cells with a measurement in that channel |
| n_negative | annotated negative cells with a measurement in that channel |
| median_positive | median intensity of the annotated positive cells |
| median_negative | median intensity of the annotated negative cells |
| applied_cut | the cut the package applies; empty for cytokeratin, whose cut is derived from the whole cohort instead, which is why its suggestion here is far from the cut in thresholds.csv |

The review is a decision aid and changes nothing else: the cuts the package applies are the named
constants at the top of mif/02_cells.R.
", nrow(cells), uniqueN(cells$`Sample Name`)),
    file.path(out_dir, "thresholds_phase_a.md"))
  cat("[A] wrote", target, "\n")
}

## ------------------------------------------------------------------------- core key files

## The BOMI2 core name, as "glass,row,letter". The annotation sheet names the same core in that
## form, which is what makes the two key files joinable.
bomi2_key <- function(name) {
  m <- regmatches(name, regexec("BOMI2_Redox_(\\d+)_Core\\[(\\d+),(\\d+),([A-Za-z])\\]", name))
  vapply(m, function(x) if (length(x) == 5) paste(x[2], x[4], toupper(x[5]), sep = ",")
                        else NA_character_, "")
}
## The BOMI1 core name, as "TMA|[column,row,letter]", the form the core map uses.
bomi1_key <- function(name) {
  m <- regmatches(name, regexec(
    "RedOx_BOMI_1_0*([0-9]+)_Core\\[\\s*([0-9]+)\\s*,\\s*([0-9]+)\\s*,\\s*([A-Za-z])\\s*\\]", name))
  vapply(m, function(x) if (length(x) == 5)
    sprintf("%d|[%s,%s,%s]", as.integer(x[2]), x[3], x[4], toupper(x[5])) else NA_character_, "")
}

##  The core table of a cohort: key, patient, and whether the core is kept.
##
##  BOMI2 excludes two cores that the core map itself marks and 52 that the annotation sheet flags,
##  and keeps only primary tissue. One patient is reclassified: the annotation sheet calls both of
##  that patient's surviving cores metastasis, a review of the tissue found them to be primary,
##  and the correction takes the cohort from 293 to 294 patients. The rule below derives the
##  correction from the sheet, because after the exclusions exactly one patient is left with
##  metastasis cores only. That keeps core identifiers out of the code. BOMI1 has one key file and
##  every core is primary; the cores of the patient its clinical record says to remove are
##  excluded (see the BOMI1 branch below).
core_table <- function() {
  if (COHORT == "BOMI2") {
    cores <- as.data.table(read_excel(file.path(CO$keys, "cores to cases_ Bomi2-Redox1.xlsx")))
    cores[, key := bomi2_key(core)]
    cores <- unique(cores[!is.na(key)], by = "key")
    cores[, dropped_by_map := as.character(Patient) == "exclude"]
    cores[, Case := fifelse(dropped_by_map, NA_character_, as.character(Patient))]

    sheet <- as.data.table(read_excel(
      file.path(CO$keys, "List of cores to exclude from Bomi2 20210602.xlsx"),
      sheet = "All annotations"))
    setnames(sheet, names(sheet)[4], "TO_EXCLUDE")
    sheet[, key := toupper(gsub(" ", "", as.character(TMA)))]
    flagged <- unique(sheet[toupper(as.character(TO_EXCLUDE)) == "YES", key])
    source_of <- unique(sheet[, .(key, Source)], by = "key")

    t <- merge(cores[, .(key, Case, dropped_by_map)], source_of, by = "key", all.x = TRUE)
    t[, TM := substr(as.character(Source), 1, 1)]
    t[, keep := !(dropped_by_map | key %in% flagged)]

    surviving <- t[keep == TRUE & !is.na(Case)]
    only_m <- surviving[, .(all_m = all(TM == "M")), by = Case][all_m == TRUE, Case]
    if (length(only_m) != 1L)
      stop(sprintf(paste("expected exactly one metastasis-only patient to reclassify, found %d.",
                         "The annotation sheet has changed; check it before trusting the cohort."),
                   length(only_m)))
    t[Case %in% only_m, TM := "T"]
    t[, keep := keep & !is.na(Case) & TM == "T"]
    cat(sprintf("[C] %d cores in the map, %d excluded, %d kept, %d patients (%d core%s of one\n",
                nrow(t), sum(!t$keep), sum(t$keep), uniqueN(t[keep == TRUE, Case]),
                nrow(surviving[Case %in% only_m]),
                if (nrow(surviving[Case %in% only_m]) == 1L) "" else "s"))
    cat("    metastasis-only patient reclassified as primary)\n")
    t[, .(key, Case, keep)]
  } else {
    cc <- as.data.table(read_excel(
      file.path(CO$keys, "Core_to_patient_matching_Bomi1_Redox1reducedpanel.xlsx")))
    setnames(cc, trimws(names(cc)))
    norm <- function(tma, id) {
      id <- toupper(gsub("\\s", "", as.character(id)))
      id <- ifelse(grepl("^\\[", id), id, paste0("[", id, "]"))
      sprintf("%d|%s", as.integer(tma), id)
    }
    cc[, key := norm(`TMA NUMBER`, `Core ID`)]
    t <- unique(cc[, .(key, Case = as.character(Case))], by = "key")
    ## The patients whose raw clinical record says to remove them are excluded with all their cores
    ## (one patient: "Remove due to mix-up" in the Comments column of the clinicopathological
    ## file). The rule reads the comment, so no identifier is written into the code.
    ## The thresholds are derived before this step, from every core, so they do not move.
    clin <- as.data.table(read_excel(file.path(CO$keys,
      "clinicopathological data of Bomi1 for the Redox1 panel cases 20200713.xls")))
    removed <- as.character(as.integer(clin[grepl("^\\s*remove", tolower(as.character(Comments))),
                                             `Project ID`]))
    t[, keep := !(Case %in% removed)]
    cat(sprintf("[C] %d cores in the map, %d excluded (clinical record: remove), %d patients\n",
                nrow(t), sum(!t$keep), uniqueN(t[keep == TRUE, Case])))
    t
  }
}

## -------------------------------------------------------------------------------- the cell table

cat(sprintf("[%s] reading %s\n", COHORT, basename(CO$cells)))
## na.strings is given explicitly because 361 BOMI1 cells carry the spreadsheet text for a missing
## measurement in every cytoplasm channel. Older versions of data.table read that text as a string
## and turned the whole channel into text, which would have made every positivity call wrong.
d <- fread(CO$cells, sep = "\t", showProgress = FALSE,
           na.strings = c("NA", "#N/A", ""),
           select = c("Sample Name", "Tissue Category", "Cell X Position", "Cell Y Position",
                      unname(CO$channels), CO$trxr1))
setnames(d, unname(CO$channels), names(CO$channels))
setnames(d, CO$trxr1, "TrxR1_intensity")
d <- d[`Tissue Category` %in% c("Tumor", "Stroma")]
cat(sprintf("[B] %d cells of the Tumor and Stroma categories, %d cores\n",
            nrow(d), uniqueN(d[["Sample Name"]])))

## ---- thresholds, on every core of the export, before any exclusion ---------------------------
##
##  The cuts are derived from all cores. Deriving them from the kept cores instead would make the
##  cut depend on the exclusion list, and the exclusion list was drawn up for reasons that have
##  nothing to do with staining intensity.

thresholds <- data.table(marker = character(), cut = numeric(), method = character(),
                         n_cells = integer())
add_cut <- function(marker, cut, method, n_cells) {
  thresholds <<- rbind(thresholds, data.table(marker = marker, cut = cut, method = method,
                                              n_cells = as.integer(n_cells)))
  cut
}

CK_cut <- add_cut("CK", ck_cut(d$CK),
                  "Otsu seed on the intensity times 1000, then the density valley between the two modes",
                  sum(is.finite(d$CK)))
if (COHORT == "BOMI2") {
  for (mk in names(BOMI2_FIXED_CUTS))
    add_cut(mk, as.numeric(BOMI2_FIXED_CUTS[[mk]]),
            "fixed constant from the Phase A review of the annotated cores", NA)
  cuts <- c(BOMI2_FIXED_CUTS, CK = CK_cut)
} else {
  CD3_cut <- add_cut("CD3",
                     as.numeric(quantile(d$CD3[d$CK > CK_cut], BOMI1_CD3_PERCENTILE, na.rm = TRUE)),
                     sprintf("percentile %.1f of CD3 within the cytokeratin-positive cells",
                             100 * BOMI1_CD3_PERCENTILE),
                     sum(is.finite(d$CD3[d$CK > CK_cut])))
  cuts <- c(CD3 = CD3_cut, CK = CK_cut)
}
TrxR1_cut <- add_cut("TrxR1",
                     crossover_cut(d[`Tissue Category` == "Tumor", TrxR1_intensity],
                                   d[`Tissue Category` == "Stroma", TrxR1_intensity]),
                     "crossover of the Tumor and Stroma intensity densities, inside the window between the two medians",
                     sum(is.finite(d$TrxR1_intensity)))
print(thresholds)

## ---- the recorded cuts must come back -------------------------------------------------------
for (mk in names(EXPECTED[[COHORT]])) {
  got <- thresholds[marker == mk, cut][1]
  want <- EXPECTED[[COHORT]][[mk]]
  rel <- abs(got - want) / want
  cat(sprintf("[B] %-5s cut %.6f, recorded %.6f, relative difference %.2e\n", mk, got, want, rel))
  if (!is.finite(rel) || rel > EXPECTED_TOLERANCE)
    stop(sprintf(paste("the %s cut of %s moved from the recorded %.6f to %.6f. Every positivity",
                       "call in the cohort depends on it; find out why before going on."),
                 mk, COHORT, want, got))
}

## ---- binary calls and lineage ----------------------------------------------------------------
##
##  A cell with no measurement in a channel is negative for that marker when assigning a lineage.
##  Its TrxR1 call stays missing, because the
##  paper averages TrxR1 over the cells that carry a call and counts the others in the cell count
##  only; writing a zero there would quietly move the exposure of every patient with such a cell.
for (mk in names(cuts))
  d[, (mk) := fifelse(is.na(get(mk)), 0L, as.integer(get(mk) > cuts[[mk]]))]
d[, TrxR1 := as.integer(TrxR1_intensity > TrxR1_cut)]

## Lineage by priority: the first marker in the cohort's order that the cell is positive for names
## it. BOMI2 runs CD3 before CD68 before CK before alpha-SMA; BOMI1 runs CD3 before CK. A cell
## positive for more than one lineage marker is assigned to the first lineage in this order.
d[, lineage := "Other"]
for (mk in rev(names(CO$lineage)))
  d[get(mk) == 1L, lineage := CO$lineage[[mk]]]
cat(sprintf("[B] positive fraction: %s\n",
            paste(sprintf("%s=%.3f", names(cuts),
                          vapply(names(cuts), function(m) mean(d[[m]]), 0)), collapse = "  ")))
cat(sprintf("[B] TrxR1 positive %.3f, %d cells without a TrxR1 measurement\n",
            mean(d$TrxR1, na.rm = TRUE), sum(is.na(d$TrxR1))))

## ---- patients, exclusions, primary tissue ----------------------------------------------------
keys <- core_table()
kept_keys <- keys[keep == TRUE]
d[, key := if (COHORT == "BOMI2") bomi2_key(`Sample Name`) else bomi1_key(`Sample Name`)]
d[kept_keys, Case := i.Case, on = "key"]
before <- nrow(d)
d <- d[!is.na(Case)]
cat(sprintf("[C] kept %d of %d cells, %d cores, %d patients\n",
            nrow(d), before, uniqueN(d$key), uniqueN(d$Case)))

## ---- the package schema ----------------------------------------------------------------------
## Positions are written as the export gives them (BOMI2 pixels, BOMI1 micrometers) and converted
## to micrometers in Python on load (redox.cells.load_cells, through redox.coords.to_um). A product
## written here as decimal text is one last bit away from pixel x 0.4964 computed in double for
## about 9 % of the BOMI2 values, and one bit is enough to break an exact tie between equidistant
## neighbors, which the six-neighbor graph and the distance rankings depend on.
## The raw values are short decimals; the check below makes sure fwrite, which writes at most 15
## significant digits, gives every one of them back exactly.
d[, x_raw := as.numeric(`Cell X Position`)]
d[, y_raw := as.numeric(`Cell Y Position`)]
for (v in c("x_raw", "y_raw")) {
  vals <- d[[v]][is.finite(d[[v]])]
  if (any(as.numeric(formatC(vals, digits = 15, format = "g")) != vals))
    stop(sprintf("some %s values need more than 15 significant digits and would not survive fwrite", v))
}
setnames(d, c("Sample Name", "Tissue Category"), c("core", "tissue"))
cells_out <- d[, c("core", "Case", "tissue", "x_raw", "y_raw", names(cuts), "lineage", "TrxR1"),
               with = FALSE]
target <- out_path("cells", "tsv")
fwrite(cells_out, target, sep = "\t")
cat(sprintf("[C] wrote %s (%d rows)\n", target, nrow(cells_out)))

## ---- tissue areas ----------------------------------------------------------------------------
ta <- fread(CO$tissue, sep = "\t", showProgress = FALSE,
            select = c("Sample Name", "Tissue Category", CO$area_column))
ta <- ta[`Tissue Category` %in% c("Tumor", "Stroma")]
## BOMI2 areas are in pixels and use the same micrometer constant as the positions; BOMI1 areas are
## in square micrometers. Both end in square millimeters.
ta[, area_mm2 := as.numeric(get(CO$area_column)) * (CO$mpp^2) / 1e6]
ta <- ta[is.finite(area_mm2) & area_mm2 > 0]
ta <- ta[`Sample Name` %in% unique(cells_out$core)]
setnames(ta, c("Sample Name", "Tissue Category"), c("core", "tissue"))
areas_out <- ta[, .(core, tissue, area_mm2)]
target <- out_path("areas", "csv")
fwrite(areas_out, target)
cat(sprintf("[D] wrote %s (%d rows, %d cores)\n", target, nrow(areas_out), uniqueN(areas_out$core)))

## ---- thresholds and the sibling notes --------------------------------------------------------
target <- out_path("thresholds", "csv")
fwrite(thresholds, target)
cat(sprintf("[D] wrote %s\n", target))

marker_rows <- paste0(vapply(names(cuts), function(m) sprintf(
  "| %s | 1 when the cell is positive for %s at the cohort cut, 0 otherwise or when the channel is missing |\n",
  m, m), ""), collapse = "")
writeLines(sprintf(
"# cells.tsv (%s)

One row per cell of the Tumor or Stroma tissue category of a primary core, after the published
exclusions. Written by mif/02_cells.R from the raw inForm export. %d rows, %d cores, %d patients.

| column | meaning |
|---|---|
| core | core identifier as in the raw export (inForm Sample Name) |
| Case | patient identifier, as a string |
| tissue | Tumor or Stroma, the inForm tissue category |
| x_raw | position along x as the export gives it, in %s |
| y_raw | position along y, same unit |
%s| lineage | priority %s; named %s |
| TrxR1 | 1 when the cytoplasmic TrxR1 intensity is above the cohort cut, 0 below, empty when the intensity is missing |

The cuts are in thresholds.csv beside this file.
", COHORT, nrow(cells_out), uniqueN(cells_out$core), uniqueN(cells_out$Case),
   if (CO$mpp == 1) "micrometers; redox.cells.load_cells returns them unchanged as x_um and y_um" else
     sprintf("pixels; redox.cells.load_cells returns x_um and y_um, the raw position times %s micrometers per pixel",
             format(CO$mpp)),
   marker_rows, paste(c(names(CO$lineage), "other"), collapse = " > "),
   paste(c(unname(CO$lineage), "Other"), collapse = ", ")),
  file.path(out_dir, "cells.md"))

writeLines(sprintf(
"# areas.csv (%s)

One row per core and tissue category. Written by mif/02_cells.R from the raw inForm tissue
summary. %d rows, %d cores.

| column | meaning |
|---|---|
| core | core identifier as in the raw export (inForm Sample Name) |
| tissue | Tumor or Stroma, the inForm tissue category |
| area_mm2 | segmented area in square millimeters |

%s
", COHORT, nrow(areas_out), uniqueN(areas_out$core),
   if (COHORT == "BOMI2")
     "The export gives the area in pixels; it is converted with 0.4964 micrometers per pixel, the one BOMI2 coordinate constant of the package, squared."
   else
     "The export gives the area in square micrometers; it is divided by a million and nothing else."),
  file.path(out_dir, "areas.md"))

writeLines(sprintf(
"# thresholds.csv (%s)

The positivity cuts this cohort uses, one row per marker. A cell is positive when its intensity is
strictly above the cut. Written by mif/02_cells.R.

| column | meaning |
|---|---|
| marker | the marker the cut belongs to |
| cut | the intensity, in normalized counts, above which a cell is positive |
| method | how the cut was obtained |
| n_cells | the cells the cut was derived from; empty for a fixed constant, whose derivation was the review of the annotated cores and not this run |

The cuts are derived from every core of the export, before any core is excluded, so that the
exclusion list cannot move a threshold.
", COHORT), file.path(out_dir, "thresholds.md"))

if (PHASE_A) {
  phase_a()
}

cat(sprintf("[%s] done\n", COHORT))
