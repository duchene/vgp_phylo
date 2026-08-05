#!/usr/bin/env Rscript
## Join every results/metrics/<name>.tsv into one wide table keyed on locus_id.
## Columns are prefixed with their metric name so two metrics may both emit
## e.g. "mean" without colliding. Loci missing from a metric get NA, so this is
## safe to re-run whenever a new metric finishes.
##
## Usage: scripts/merge_metrics.R [--out PATH] [--drop-errors]

args <- commandArgs(trailingOnly = TRUE)
getopt <- function(f, d = NULL) { i <- match(f, args); if (is.na(i)) d else args[i + 1L] }

root <- Sys.getenv("VGP_PHYLO_ROOT", unset = normalizePath(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "..")))
source(file.path(root, "config", "project.R"))
source(file.path(root, "scripts", "lib", "load.R"))

out_path    <- getopt("--out", PATHS$merged)
drop_errors <- "--drop-errors" %in% args

files <- list.files(PATHS$results, pattern = "\\.tsv$", full.names = TRUE)
if (!length(files)) stop("no metric outputs in ", PATHS$results)

base <- read_manifest()[, c("locus_id", "clade", "ntips", "nsites")]
merged <- base

for (f in files) {
  nm <- sub("\\.tsv$", "", basename(f))
  d <- utils::read.delim(f, stringsAsFactors = FALSE)
  if (!"locus_id" %in% colnames(d)) { warning("skipping ", basename(f), ": no locus_id"); next }
  nerr <- if ("error" %in% colnames(d)) sum(!is.na(d$error)) else 0L
  vals <- setdiff(colnames(d), c("locus_id", "clade", "error"))
  d <- d[, c("locus_id", vals), drop = FALSE]
  colnames(d)[-1] <- paste0(nm, ".", vals)
  merged <- merge(merged, d, by = "locus_id", all.x = TRUE, sort = FALSE)
  message(sprintf("%-26s %6d rows, %2d cols, %d errors", nm, nrow(d), length(vals), nerr))
}

## By default keep only loci that have at least one metric computed. Without this,
## a partial run (say one clade) yields all 112k manifest rows, almost all of them
## entirely NA. Pass --all-loci to get the full manifest regardless.
meta_cols <- c("locus_id", "clade", "ntips", "nsites")
if (!"--all-loci" %in% args) {
  val <- merged[, setdiff(colnames(merged), meta_cols), drop = FALSE]
  has_any <- rowSums(!is.na(val)) > 0L
  message(sprintf("keeping %d of %d loci with >= 1 metric computed",
                  sum(has_any), nrow(merged)))
  merged <- merged[has_any, , drop = FALSE]
}

if (drop_errors) {
  n0 <- nrow(merged)
  merged <- merged[stats::complete.cases(merged), , drop = FALSE]
  message(sprintf("dropped %d loci with any missing metric", n0 - nrow(merged)))
}

merged <- merged[order(merged$clade, merged$locus_id), , drop = FALSE]
utils::write.table(merged, out_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
message(sprintf("wrote %s (%d loci x %d columns)", out_path, nrow(merged), ncol(merged)))
