#!/usr/bin/env Rscript
## Stage 0 (minimal): where each locus sits.
##
##   1. Along the ANCESTRAL contig, from the locus name <Anc>_<Anc>refChr<N>_<start>_<end>:
##      contig, start, end, mid, rank along the contig, loci on the contig, the span
##      the contig's loci cover, and distance from the locus midpoint to the outermost
##      locus on either side. The contigs' real lengths are unknown, so "end" means the
##      outermost locus, not the contig end.
##   2. On a REFERENCE SPECIES' chromosome, from the chromosome accession in that
##      species' sequence label (results/stage0/locus_labels.tsv.gz) and its NCBI
##      sequence report (results/stage0/assembly_reports/<assembly>.tsv):
##        ref_*      the clade's reference, chosen automatically: the candidate assembly
##                   with the most loci on assembled chromosomes (ties: most loci).
##        chicken_*  GCF_016700215.2, with a size class (macro 1-5, Z, W; intermediate
##                   6-10; micro otherwise -- the usual chicken karyotype split).
##        human_*    GRCh38 (GCA_000001405.15).
##      Labels give the chromosome, not the position on it.
##
## Usage:
##   scripts/stage0_positions.R --candidates [n]   # print the top-n assemblies per clade
##                                                 # by locus coverage (fetch their reports)
##   scripts/stage0_positions.R                    # write results/stage0/locus_position.tsv
##                                                 # and results/stage0/references.tsv

args <- commandArgs(trailingOnly = TRUE)
root <- Sys.getenv("VGP_PHYLO_ROOT", unset = normalizePath(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "..")))
setwd(root)

CHICKEN <- "GCF_016700215.2"
HUMAN   <- "GCA_000001405.15"
REPORTS <- "results/stage0/assembly_reports"

m  <- utils::read.delim("data/manifest.tsv", stringsAsFactors = FALSE)
lb <- utils::read.delim(gzfile("results/stage0/locus_labels.tsv.gz"), stringsAsFactors = FALSE,
                        colClasses = "character")

## Coverage: in how many of the clade's loci each assembly occurs.
cov <- stats::aggregate(locus_id ~ clade + assembly, data = unique(lb[, c("locus_id", "clade", "assembly")]),
                        FUN = length)
names(cov)[3] <- "n_loci"
cov <- cov[order(cov$clade, -cov$n_loci), ]

if (length(args) && args[1] == "--candidates") {
  n <- if (length(args) >= 2) as.integer(args[2]) else 5L
  top <- do.call(rbind, lapply(split(cov, cov$clade), utils::head, n))
  cat(unique(c(top$assembly, CHICKEN, HUMAN)), sep = "\n")
  quit(status = 0)
}

read_report <- function(a) {
  f <- file.path(REPORTS, paste0(a, ".tsv"))
  if (!file.exists(f)) return(NULL)
  r <- utils::read.delim(f, stringsAsFactors = FALSE, colClasses = "character")
  r <- r[r$role == "assembled-molecule", ]
  ## A label may carry either the GenBank or the RefSeq accession.
  rbind(data.frame(seq_accession = r$genbank_accession, chr = r$chr_name, chr_len = as.numeric(r$length)),
        data.frame(seq_accession = r$refseq_accession, chr = r$chr_name, chr_len = as.numeric(r$length)))
}

## For each locus and a given assembly: that species' chromosome, if on an assembled one.
on_chr <- function(a) {
  rep <- read_report(a)
  if (is.null(rep)) stop("no sequence report for ", a, ": run the --candidates step and fetch it")
  x <- lb[lb$assembly == a, c("locus_id", "seq_accession")]
  x <- x[!duplicated(x$locus_id), ]
  k <- match(x$seq_accession, rep$seq_accession)
  data.frame(locus_id = x$locus_id, chr = rep$chr[k], chr_len = rep$chr_len[k],
             stringsAsFactors = FALSE)
}

## ---- reference per clade
cands <- split(cov, cov$clade)
refs <- do.call(rbind, lapply(names(cands), function(cl) {
  cc <- cands[[cl]]
  cc <- cc[file.exists(file.path(REPORTS, paste0(cc$assembly, ".tsv"))), ]
  if (!nrow(cc)) return(data.frame(clade = cl, assembly = NA, n_loci = 0L, n_on_chr = 0L, n_clade = NA))
  cc$n_on_chr <- vapply(cc$assembly, function(a) {
    o <- on_chr(a); sum(!is.na(o$chr) & o$locus_id %in% m$locus_id[m$clade == cl])
  }, integer(1))
  best <- cc[order(-cc$n_on_chr, -cc$n_loci), ][1, ]
  data.frame(clade = cl, assembly = best$assembly, n_loci = best$n_loci, n_on_chr = best$n_on_chr,
             n_clade = sum(m$clade == cl))
}))
utils::write.table(refs, "results/stage0/references.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(refs, row.names = FALSE)

## ---- ancestral contig coordinates
p <- regmatches(m$locus_id, regexec("refChr([0-9]+)_([0-9]+)_([0-9]+)$", m$locus_id))
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

## ---- reference chromosomes
add <- function(pos, o, prefix) {
  k <- match(pos$locus_id, o$locus_id)
  pos[[paste0(prefix, "_chr")]] <- o$chr[k]
  pos[[paste0(prefix, "_chr_len")]] <- o$chr_len[k]
  pos
}
pos$ref_assembly <- refs$assembly[match(pos$clade, refs$clade)]
ref_o <- do.call(rbind, lapply(unique(stats::na.omit(refs$assembly)), function(a) {
  o <- on_chr(a); o[o$locus_id %in% pos$locus_id[pos$ref_assembly %in% a], ]
}))
pos <- add(pos, ref_o, "ref")
pos <- add(pos, on_chr(CHICKEN), "chicken")
pos$chicken_class <- ifelse(is.na(pos$chicken_chr), NA,
  ifelse(pos$chicken_chr %in% c(as.character(1:5), "Z", "W"), "macro",
  ifelse(pos$chicken_chr %in% as.character(6:10), "intermediate", "micro")))
pos <- add(pos, on_chr(HUMAN), "human")

utils::write.table(pos, "results/stage0/locus_position.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nlocus_position.tsv:", nrow(pos), "loci;",
    sum(!is.na(pos$ref_chr)), "on a reference chromosome;",
    sum(!is.na(pos$chicken_chr)), "on a chicken chromosome;",
    sum(!is.na(pos$human_chr)), "on a human chromosome\n")
