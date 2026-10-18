#!/usr/bin/env Rscript
## Collapse identical and near-identical sequences before ModelFinder, for the loci
## of one shard, and prune each fixed gene tree to the sequences kept.
##
## WHY: a near-duplicate adds a near-zero terminal branch and no information, but
## still costs 2 branch-length parameters (and k of them under +H). They inflate the
## BIC penalty of the tip-heavy models and slow every fit.
##
## Rule (deterministic): sequences are visited from most to fewest unambiguous bases
## (ties by label); one is dropped if its uncorrected p-distance to an already-kept
## sequence is <= MF_NEARIDENT (default 0.001, i.e. identical or ~2 differences in
## 2500 sites), computed over the sites both have unambiguous (pairwise deletion).
## Pairs with no shared sites never count as identical. At least 4 sequences are
## always kept. Measured on the 187-locus panel: drops 1,108 of 35,494 tips (3.1%)
## in 121 loci; p <= 0 alone drops 789.
##
## Writes <scratch>/<locus>.fa (kept records, verbatim) and <scratch>/<locus>.tre,
## appends one row per locus to <outdir>/dedup/shard_NNN.tsv, and prints one job line
## per locus for modelfinder_shard.sh: id clade ntips_kept nsites fasta tree seed.
## Loci already finished (<outdir>/<clade>/<locus>.iqtree.gz) are skipped.
##
## Usage: dedup_loci.R <jobs.tsv> <shard> <nshards> <scratch> <outdir>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop("usage: dedup_loci.R <jobs.tsv> <shard> <nshards> <scratch> <outdir>")
jobs <- utils::read.delim(args[1], stringsAsFactors = FALSE)
i <- as.integer(args[2]); n <- as.integer(args[3]); scratch <- args[4]; outdir <- args[5]
thr <- as.numeric(Sys.getenv("MF_NEARIDENT", unset = "0.001"))

jobs <- jobs[(seq_len(nrow(jobs)) - 1L) %% n == i - 1L, ]
done <- file.exists(file.path(outdir, jobs$clade, paste0(jobs$locus_id, ".iqtree.gz")))
jobs <- jobs[!done, ]

tab <- file.path(outdir, "dedup", sprintf("shard_%03d.tsv", i))
dir.create(dirname(tab), recursive = TRUE, showWarnings = FALSE)
if (!file.exists(tab))
  cat("locus_id\tntips\tntips_kept\tthreshold\tdropped\n", file = tab)

read_fasta_text <- function(f) {   # named character vector of verbatim records
  ln <- readLines(f, warn = FALSE)
  h <- grep("^>", ln)
  ends <- c(h[-1] - 1L, length(ln))
  rec <- vapply(seq_along(h), function(k) paste(ln[h[k]:ends[k]], collapse = "\n"), "")
  names(rec) <- sub("\\s.*", "", sub("^>", "", ln[h]))
  rec
}

for (r in seq_len(nrow(jobs))) {
  j <- jobs[r, ]
  a <- ape::read.dna(j$fasta, format = "fasta")
  lab <- rownames(a)
  ch <- as.character(a)
  info <- rowSums(matrix(ch %in% c("a", "c", "g", "t"), nrow = nrow(ch)))
  d <- as.matrix(ape::dist.dna(a, model = "raw", pairwise.deletion = TRUE))
  ord <- order(-info, lab)
  keep <- character(); repr <- setNames(rep(NA_character_, length(lab)), lab)
  for (s in lab[ord]) {
    hit <- if (length(keep)) keep[which(d[s, keep] <= thr)] else character()
    if (length(hit)) repr[s] <- hit[1] else keep <- c(keep, s)
  }
  if (length(keep) < 4L) keep <- lab[ord][seq_len(min(4L, length(lab)))]
  drop <- setdiff(lab, keep)

  fa <- file.path(scratch, paste0(j$locus_id, ".fa"))
  tre <- file.path(scratch, paste0(j$locus_id, ".tre"))
  rec <- read_fasta_text(j$fasta)
  stopifnot(setequal(names(rec), lab))
  writeLines(rec[keep], fa)
  tr <- ape::read.tree(j$treefile)
  stopifnot(setequal(tr$tip.label, lab))
  if (length(drop)) tr <- ape::drop.tip(tr, drop)
  ape::write.tree(tr, tre)

  cat(sprintf("%s\t%d\t%d\t%g\t%s\n", j$locus_id, length(lab), length(keep), thr,
              if (length(drop)) paste0(drop, ">", repr[drop], collapse = ",") else ""),
      file = tab, append = TRUE)
  ## No trailing blank: xargs -L treats one as a line continuation.
  cat(paste(j$locus_id, j$clade, length(keep), j$nsites, fa, tre, j$seed), "\n", sep = "")
  flush(stdout())   # streamed: IQ-TREE starts on this locus while the next is prepared
}
