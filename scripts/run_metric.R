#!/usr/bin/env Rscript
## Compute ONE metric across loci and append tidy rows to results/metrics/<name>.tsv.
##
## Parallelism has two independent axes:
##   within a metric  -- this script forks --cores workers over loci
##   across metrics   -- run several copies of this script (see run_metrics.sh)
## Metrics never share state or output files, so the second axis is embarrassingly
## parallel and a new metric costs nothing to add.
##
## Usage:
##   scripts/run_metric.R --metric <name> [options]
##     --cores N        workers within this metric        (default DEFAULTS$cores)
##     --clade C        restrict to one clade             (repeatable, comma-separated)
##     --limit N        first N loci only (smoke testing)
##     --shard i/N      process shard i of N              (for splitting across hosts)
##     --block N        loci per flush / resume unit      (default DEFAULTS$block)
##     --out PATH       override output path
##     --no-resume      recompute loci already in the output
##     --seed N         RNG seed base, for quartet-sampling metrics (default 1)

suppressWarnings(suppressPackageStartupMessages({
  library(parallel)
}))

# ------------------------------------------------------------------- args ----

args <- commandArgs(trailingOnly = TRUE)
getopt <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop(flag, " needs a value")
  args[i + 1L]
}
hasflag <- function(flag) flag %in% args

metric_name <- getopt("--metric")
if (is.null(metric_name)) stop("--metric is required")

root <- Sys.getenv("VGP_PHYLO_ROOT", unset = normalizePath(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "..")))
source(file.path(root, "config", "project.R"))
source(file.path(root, "scripts", "lib", "load.R"))

cores  <- as.integer(getopt("--cores", DEFAULTS$cores))
block  <- as.integer(getopt("--block", DEFAULTS$block))
limit  <- getopt("--limit"); limit <- if (is.null(limit)) NA_integer_ else as.integer(limit)
clades <- getopt("--clade")
shard  <- getopt("--shard")
seed   <- as.integer(getopt("--seed", 1L))
resume <- !hasflag("--no-resume")

# ------------------------------------------------------- metric and loci ----

m <- load_metric(metric_name)
out_path <- getopt("--out", file.path(PATHS$results, paste0(m$name, ".tsv")))
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
dir.create(PATHS$logs, recursive = TRUE, showWarnings = FALSE)

loci <- read_manifest()
if (!is.null(clades)) {
  keep <- trimws(strsplit(clades, ",")[[1]])
  bad <- setdiff(keep, unique(loci$clade))
  if (length(bad)) stop("unknown clade(s): ", paste(bad, collapse = ", "))
  loci <- loci[loci$clade %in% keep, , drop = FALSE]
}
if (!is.null(shard)) {
  p <- as.integer(strsplit(shard, "/")[[1]])
  if (length(p) != 2L || is.na(p[1]) || is.na(p[2]) || p[1] < 1L || p[1] > p[2])
    stop("--shard must look like i/N with 1 <= i <= N")
  loci <- loci[(seq_len(nrow(loci)) - 1L) %% p[2] == (p[1] - 1L), , drop = FALSE]
}
if (!is.na(limit)) loci <- utils::head(loci, limit)

header <- c("locus_id", "clade", m$columns, "error")
if (resume && file.exists(out_path)) {
  done <- utils::read.delim(out_path, stringsAsFactors = FALSE)
  if (!identical(colnames(done), header)) {
    stop("existing ", basename(out_path), " has a different schema than the metric now",
         " declares.\n  on disk: ", paste(colnames(done), collapse = ","),
         "\n  metric:  ", paste(header, collapse = ","),
         "\nMove or delete the file, or pass --no-resume to overwrite.")
  }
  n0 <- nrow(loci)
  loci <- loci[!loci$locus_id %in% done$locus_id, , drop = FALSE]
  message(sprintf("[%s] resuming: %d of %d loci already done", m$name, n0 - nrow(loci), n0))
} else if (file.exists(out_path)) {
  file.remove(out_path)
}

if (nrow(loci) == 0L) {
  message(sprintf("[%s] nothing to do", m$name))
  quit(save = "no", status = 0)
}

# ----------------------------------------------------------------- compute ----

## Loaded once per process, not per locus: reference trees are shared and immutable.
if ("species_tree" %in% m$needs) sptree <- load_species_tree(m$species_tree)

na_row <- setNames(rep(list(NA_real_), length(m$columns)), m$columns)

one_locus <- function(i) {
  loc <- as.list(loci[i, , drop = FALSE])
  ## Deterministic per-locus seed: quartet-sampling metrics must be reproducible
  ## and must not depend on how loci were sharded across workers.
  set.seed(seed + sum(utils::head(utf8ToInt(loc$locus_id), 64L)))
  res <- tryCatch({
    if ("tree" %in% m$needs)         loc$tree        <- read_gene_tree(loc$treefile)
    if ("alignment" %in% m$needs)    loc$alignment   <- read_locus_alignment(loc$fasta)
    if ("species_tree" %in% m$needs) loc$species_tree <- sptree
    v <- m$compute(loc)
    miss <- setdiff(m$columns, names(v))
    if (length(miss))
      stop("compute() omitted declared column(s): ", paste(miss, collapse = ", "))
    c(v[m$columns], list(error = NA_character_))
  }, error = function(e) c(na_row, list(error = gsub("[\t\r\n]+", " ", conditionMessage(e)))))
  c(list(locus_id = loc$locus_id, clade = loc$clade), res)
}

write_block <- function(rows, append) {
  df <- do.call(rbind, lapply(rows, function(r) as.data.frame(r, stringsAsFactors = FALSE)))
  utils::write.table(df, out_path, sep = "\t", quote = FALSE, row.names = FALSE,
                     col.names = !append, append = append, na = "NA")
}

message(sprintf("[%s] %d loci, %d cores, blocks of %d -> %s",
                m$name, nrow(loci), cores, block, out_path))
t0 <- Sys.time()
appending <- file.exists(out_path)
starts <- seq(1L, nrow(loci), by = block)

for (s in starts) {
  idx <- s:min(s + block - 1L, nrow(loci))
  rows <- if (cores > 1L) mclapply(idx, one_locus, mc.cores = cores, mc.preschedule = TRUE)
          else lapply(idx, one_locus)
  ## mclapply returns a try-error object if a worker dies outright (e.g. killed by
  ## the OOM killer); those loci would otherwise vanish from the output silently.
  bad <- vapply(rows, inherits, logical(1), "try-error")
  if (any(bad)) rows[bad] <- lapply(which(bad), function(k) c(
    list(locus_id = loci$locus_id[idx[k]], clade = loci$clade[idx[k]]),
    na_row, list(error = "worker died")))
  write_block(rows, append = appending)
  appending <- TRUE
  nerr <- sum(!is.na(vapply(rows, function(r) r$error, character(1))))
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  message(sprintf("[%s] %d/%d  (%.1f loci/s, %d errors this block)",
                  m$name, max(idx), nrow(loci), length(idx) / max(dt, 1e-9), nerr))
  t0 <- Sys.time()
}

fin <- utils::read.delim(out_path, stringsAsFactors = FALSE)
message(sprintf("[%s] done: %d rows, %d with errors", m$name, nrow(fin), sum(!is.na(fin$error))))
