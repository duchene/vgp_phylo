#!/usr/bin/env Rscript
## Reproducibility check: compare metrics recomputed on a panel (e.g. after moving
## machines or changing IQ-TREE version) with a reference run, column by column.
##
## Usage: Rscript scripts/compare_panel.R [panel_dir] [reference_dir]
##   defaults: results/panel_187  results/metrics
## Writes <panel_dir>/comparison.tsv and prints it.
##
## Columns: n (loci compared), n_na_mismatch (NA in one run only), identical
## (bit-for-bit equal), max_abs, max_rel (relative to |reference|, zeros skipped),
## spearman. Deterministic metrics should be identical; IQ-TREE-based ones may
## differ in the last digits across versions, which max_rel quantifies.

root <- Sys.getenv("VGP_PHYLO_ROOT", unset = normalizePath(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "..")))
args <- commandArgs(trailingOnly = TRUE)
pdir <- if (length(args) >= 1) args[1] else file.path(root, "results", "panel_187")
rdir <- if (length(args) >= 2) args[2] else file.path(root, "results", "metrics")

rows <- list()
for (f in sort(list.files(pdir, pattern = "\\.tsv$", full.names = TRUE))) {
  metric <- sub("\\.tsv$", "", basename(f))
  if (metric == "comparison") next
  rf <- file.path(rdir, basename(f))
  if (!file.exists(rf)) { message("no reference for ", metric); next }
  p <- utils::read.delim(f, stringsAsFactors = FALSE)
  r <- utils::read.delim(rf, stringsAsFactors = FALSE)
  r <- r[match(p$locus_id, r$locus_id), ]
  if (anyNA(r$locus_id)) stop(metric, ": panel loci missing from the reference")
  for (col in setdiff(colnames(p), c("locus_id", "clade", "error"))) {
    a <- p[[col]]; b <- r[[col]]
    both <- !is.na(a) & !is.na(b)
    d <- abs(a[both] - b[both])
    rel <- d[b[both] != 0] / abs(b[both][b[both] != 0])
    rows[[length(rows) + 1L]] <- data.frame(
      metric = metric, column = col, n = sum(both),
      n_na_mismatch = sum(is.na(a) != is.na(b)),
      identical = isTRUE(all(a[both] == b[both])),
      max_abs = if (length(d)) max(d) else NA_real_,
      max_rel = if (length(rel)) max(rel) else NA_real_,
      spearman = if (sum(both) > 2 && stats::sd(b[both]) > 0)
        suppressWarnings(stats::cor(a[both], b[both], method = "spearman")) else NA_real_)
  }
  nerr <- sum(!is.na(p$error))
  if (nerr) message(metric, ": ", nerr, " loci with errors, e.g. ", p$error[!is.na(p$error)][1])
}
out <- do.call(rbind, rows)
utils::write.table(out, file.path(pdir, "comparison.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
num <- vapply(out, is.numeric, logical(1)); out[num] <- lapply(out[num], signif, 3)
print(out, row.names = FALSE)
