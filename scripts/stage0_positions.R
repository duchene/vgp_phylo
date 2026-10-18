#!/usr/bin/env Rscript
## Stage 0 (minimal): where each locus sits in its clade's ANCESTRAL genome, read
## from the locus name <Anc>_<Anc>refChr<N>_<start>_<end>. Nothing else is used:
## the ancestral reconstructions may still change substantially, so this stays a
## thin, rebuildable layer (other frames are documented in docs/PLAN, Step 1).
##
## Columns: contig, start, end, mid; rank of the locus along its contig; number of
## loci on the contig; the span the contig's loci cover; and the distance from the
## locus midpoint to the outermost locus on either side. Contig lengths are not
## known, so "end" means the outermost locus, not the contig end.
##
## Usage: Rscript scripts/stage0_positions.R   -> results/stage0/locus_position.tsv

root <- Sys.getenv("VGP_PHYLO_ROOT", unset = normalizePath(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "..")))
source(file.path(root, "config", "project.R"))

m <- utils::read.delim(PATHS$manifest, stringsAsFactors = FALSE)
p <- regmatches(m$locus_id, regexec("refChr([0-9]+)_([0-9]+)_([0-9]+)$", m$locus_id))
if (any(lengths(p) != 4L)) stop("unparseable locus names, e.g. ", m$locus_id[lengths(p) != 4L][1])

pos <- data.frame(locus_id = m$locus_id, clade = m$clade,
                  contig = paste0("refChr", vapply(p, `[`, "", 2)),
                  start = as.numeric(vapply(p, `[`, "", 3)), end = as.numeric(vapply(p, `[`, "", 4)),
                  stringsAsFactors = FALSE)
pos$mid <- (pos$start + pos$end) / 2
key <- paste(pos$clade, pos$contig)
pos$n_loci_contig <- as.integer(stats::ave(pos$mid, key, FUN = length))
pos$rank <- as.integer(stats::ave(pos$mid, key, FUN = function(x) rank(x, ties.method = "first")))
pos$span_start <- stats::ave(pos$start, key, FUN = min)
pos$span_end <- stats::ave(pos$end, key, FUN = max)
pos$dist_to_outer <- pmin(pos$mid - pos$span_start, pos$span_end - pos$mid)

out <- file.path(PATHS$root, "results", "stage0", "locus_position.tsv")
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
utils::write.table(pos, out, sep = "\t", quote = FALSE, row.names = FALSE)
cat(nrow(pos), "loci on", length(unique(key)), "ancestral contigs ->", out, "\n")
