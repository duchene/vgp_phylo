#!/usr/bin/env Rscript
## Write the job table for the ModelFinder run: one row per locus with its
## alignment, fixed gene tree and deterministic IQ-TREE seed (locus_seed()).
## The seed is computed here, in R, so the bash array tasks use exactly the same
## per-locus seed as every other IQ-TREE run in the project.
##
## Usage: scripts/modelfinder_jobs.R <out.tsv> [loci_file]
##   loci_file  restrict to these locus_ids (first column, header ok), e.g.
##              data/panel_187.tsv; default all loci in data/manifest.tsv.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("usage: modelfinder_jobs.R <out.tsv> [loci_file]")

root <- Sys.getenv("VGP_PHYLO_ROOT", unset = normalizePath(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "..")))
source(file.path(root, "config", "project.R"))
source(file.path(root, "scripts", "lib", "load.R"))

m <- read_manifest()
if (length(args) >= 2L) {
  ids <- utils::read.delim(args[2], stringsAsFactors = FALSE)[[1]]
  miss <- setdiff(ids, m$locus_id)
  if (length(miss)) stop(length(miss), " locus_id(s) not in the manifest, e.g. ", miss[1])
  m <- m[match(ids, m$locus_id), ]
}
stopifnot(all(file.exists(m$fasta)), all(file.exists(m$treefile)))

m$seed <- vapply(m$locus_id, locus_seed, integer(1))
if (anyNA(m$seed)) stop("IQTREE_SEED is unset; this run must be seeded")

out <- m[, c("locus_id", "clade", "ntips", "nsites", "fasta", "treefile", "seed")]
utils::write.table(out, args[1], sep = "\t", quote = FALSE, row.names = FALSE)
cat(nrow(out), "jobs ->", args[1], "\n")
