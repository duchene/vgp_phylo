#!/usr/bin/env Rscript
## Build data/panel_187.tsv: the fixed pilot panel for porting checks and Stage 1
## timing (docs/PLAN_genome_landscapes.md).
##
## The original "187 loci spanning all clades" used for the measurements quoted in
## config/project.R and metrics/README.md was never written to disk. This panel
## REPLACES it, deterministically (no RNG), so it can be rebuilt exactly:
##
##   - 9 loci per clade, at evenly spaced quantiles of that clade's tip count
##     (all loci for the two 2-locus clades)                        -> 184
##   - plus the 3 largest remaining loci overall by ntips x nsites,
##     the worst cases for timing IQ-TREE model fits                 ->   3
##
## Usage: Rscript scripts/make_panel.R

root <- Sys.getenv("VGP_PHYLO_ROOT", unset = normalizePath(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "..")))
source(file.path(root, "config", "project.R"))
source(file.path(root, "scripts", "lib", "load.R"))

PER_CLADE <- 9L
N_EXTREME <- 3L

m <- read_manifest()
m <- m[order(m$clade, m$ntips, m$nsites, m$locus_id), ]

pick <- unlist(lapply(split(seq_len(nrow(m)), m$clade), function(ix) {
  if (length(ix) <= PER_CLADE) return(ix)
  ix[unique(round(seq(1, length(ix), length.out = PER_CLADE)))]
}))
rest <- setdiff(seq_len(nrow(m)), pick)
size <- m$ntips[rest] * m$nsites[rest]
pick <- c(pick, rest[order(-size, m$locus_id[rest])][seq_len(N_EXTREME)])

p <- m[sort(pick), c("locus_id", "clade", "ntips", "nsites")]
p$role <- ifelse(seq_len(nrow(p)) %in% match(m$locus_id[tail(pick, N_EXTREME)], p$locus_id),
                 "largest", "clade_quantile")
out <- file.path(PATHS$root, "data", "panel_187.tsv")
utils::write.table(p, out, sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("%d loci, %d clades, ntips %d-%d, nsites %d-%d -> %s",
                nrow(p), length(unique(p$clade)), min(p$ntips), max(p$ntips),
                min(p$nsites), max(p$nsites), out))
