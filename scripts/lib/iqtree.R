## Running IQ-TREE per locus for site concordance factors, as cheaply as possible.
##
## DESIGN CONSTRAINT: we do not re-infer branch lengths and we do not select a
## substitution model. Both would be prohibitive over 112k loci, and neither is
## needed -- the gene trees here were already inferred with IQ-TREE, so their
## topology AND branch lengths are taken as given.
##
##   * -m JC is fixed for every locus. JC has NO free parameters, so IQ-TREE
##     estimates nothing. This is also a comparability requirement: if the model
##     varied per locus, sCF would partly reflect which model won. Measured cost
##     and effect on this data (largest locus, 557 tips): GTR+G4 5.7s vs JC 1.4s,
##     with mean sCF differing by 0.16 percentage points (68.06 vs 67.90). sCF is
##     a count of sites favouring each quartet resolution, so it is almost
##     insensitive to the model -- the cheap model costs essentially no accuracy.
##   * -blfix is used whenever the reference tree already carries usable branch
##     lengths, i.e. for the locus's own gene tree. Nothing is optimised at all.
##   * -blfix is NOT used for the species tree: the ASTRAL trees' branch lengths
##     are in COALESCENT units, not substitutions/site, so they cannot be held
##     fixed in a likelihood calculation. Only branch lengths are optimised there
##     (still under JC), which measured 1.7s on the largest locus.
##
## Other hard-won details:
##   * --seqtype DNA is required; auto-detection fails on these alignments and
##     IQ-TREE exits with "Unknown sequence type".
##   * Branches with sN == 0 (no informative sites in any sampled quartet) are
##     reported as sCF = 0. That is missing data, not concordance of zero, so those
##     rows are dropped -- averaging them in would bias sCF downwards precisely for
##     the shortest branches.
##   * A fixed --scfl quartet count keeps per-branch sampling effort constant
##     across loci of very different size.

#' Run IQ-TREE --scfl on one locus against a fixed reference topology.
#'
#' @param tree an ape phylo, or a path to a Newick file.
#' @param blfix TRUE to hold the reference tree's branch lengths fixed (use for
#'   gene trees); FALSE to optimise them (use for the ASTRAL species trees).
#' @return list(cf = data.frame, tree = phylo or NULL)
run_iqtree_scf <- function(fasta, tree, blfix, nquartets = 100L, threads = 1L,
                           model = IQTREE_MODEL, iqtree = IQTREE) {
  stopifnot(file.exists(fasta), is.logical(blfix), length(blfix) == 1L)
  if (!file.exists(iqtree)) stop("IQ-TREE binary not found: ", iqtree)

  dir <- tempfile("iqt_")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)

  tpath <- if (inherits(tree, "phylo")) {
    p <- file.path(dir, "ref.tre"); ape::write.tree(tree, p); p
  } else tree

  prefix <- file.path(dir, "run")
  res <- system2(iqtree,
    c("-s", shQuote(fasta), "-te", shQuote(tpath), "--scfl", nquartets,
      "--seqtype", "DNA", "-m", model, if (blfix) "-blfix",
      "--prefix", shQuote(prefix), "-T", threads, "-quiet"),
    stdout = TRUE, stderr = TRUE)
  status <- attr(res, "status")

  stat <- paste0(prefix, ".cf.stat")
  if (!file.exists(stat))
    stop("IQ-TREE produced no cf.stat (exit ", if (is.null(status)) 0 else status, "): ",
         paste(utils::tail(res, 3), collapse = " | "))

  tab <- utils::read.table(stat, header = TRUE, sep = "\t", comment.char = "#",
                           stringsAsFactors = FALSE)
  need <- c("ID", "sCF", "sDF1", "sDF2", "sN")
  miss <- setdiff(need, colnames(tab))
  if (length(miss)) stop("cf.stat missing column(s): ", paste(miss, collapse = ", "))

  tab <- tab[is.finite(tab$sN) & tab$sN > 0, , drop = FALSE]
  for (cc in c("sCF", "sDF1", "sDF2")) tab[[cc]] <- tab[[cc]] / 100
  keep <- c("ID", "sCF", "sCF_N", "sDF1", "sDF1_N", "sDF2", "sDF2_N", "sN", "Length")
  tab <- tab[, intersect(keep, colnames(tab)), drop = FALSE]

  tfile <- paste0(prefix, ".treefile")
  otree <- if (file.exists(tfile)) try(ape::read.tree(tfile), silent = TRUE) else NULL
  if (inherits(otree, "try-error")) otree <- NULL

  ## Fitted model parameters, scraped from the .iqtree report. These are estimated
  ## anyway as part of the --scfl run, so capturing them is free.
  ##
  ## Under HKY the "Rate parameter R" block fixes all transversions at 1.0 and
  ## reports the shared transition rate on the A-G (== C-T) line: that value IS
  ## kappa. "Gamma shape alpha" is the +G4 rate-heterogeneity parameter.
  kappa <- NA_real_; alpha <- NA_real_
  rfile <- paste0(prefix, ".iqtree")
  if (file.exists(rfile)) {
    ln <- readLines(rfile, warn = FALSE)
    ag <- grep("^\\s*A-G:", ln, value = TRUE)
    if (length(ag)) kappa <- suppressWarnings(as.numeric(sub(".*A-G:\\s*", "", ag[1])))
    gs <- grep("Gamma shape alpha:", ln, value = TRUE)
    if (length(gs)) alpha <- suppressWarnings(as.numeric(sub(".*alpha:\\s*", "", gs[1])))
  }

  list(cf = tab, tree = otree, kappa = kappa, alpha = alpha)
}

#' Cached form of run_iqtree_scf, keyed on (tag, clade, locus_id).
#'
#' IQ-TREE is the dominant cost of the analysis and several metrics are functions
#' of the same per-branch table, so it is computed once per (locus, reference).
#'
#' @param tag "own" or "species" -- which reference topology, so the caches for
#'   the two never collide.
iqtree_locus <- function(loc, tree, tag, blfix, nquartets = 100L, threads = 1L) {
  ## The model and quartet count are part of the cache KEY, not just the payload.
  ## Without this, changing IQTREE_MODEL would silently serve results computed
  ## under the previous model -- a stale-cache bug that is invisible in the output.
  key <- sprintf("%s_m%s_q%d", tag, gsub("[^A-Za-z0-9]", "", IQTREE_MODEL), nquartets)
  f <- file.path(PATHS$scf_cache, key, loc$clade, paste0(loc$locus_id, ".rds"))
  if (file.exists(f)) {
    v <- try(readRDS(f), silent = TRUE)
    ## A truncated file from an interrupted run must not poison every metric that
    ## reads it; fall through and recompute instead.
    if (!inherits(v, "try-error") && is.list(v) && !is.null(v$cf)) return(v)
  }
  v <- run_iqtree_scf(loc$fasta, tree, blfix = blfix,
                      nquartets = nquartets, threads = threads)
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(f, ".tmp", Sys.getpid())
  saveRDS(v, tmp, compress = TRUE)
  file.rename(tmp, f)   # atomic: concurrent workers never see a partial cache entry
  v
}

#' Prune a species tree to the taxa present in a gene tree / alignment and
#' relabel its tips to the locus's own sequence labels, so IQ-TREE can use it as
#' a fixed topology for that alignment.
#'
#' Gene-tree labels are "<accession>.<chr>|copy<N>" while species-tree tips are
#' bare accessions; these loci are single-copy (verified in the manifest) so the
#' mapping is 1:1.
species_tree_for_locus <- function(sptree, labels) {
  acc <- accession_from_label(labels)
  if (anyDuplicated(acc)) stop("locus has >1 sequence for some species; expected single-copy")
  shared <- intersect(sptree$tip.label, acc)
  if (length(shared) < 4L)
    stop("only ", length(shared), " taxa shared with the species tree; need >= 4")
  tr <- ape::keep.tip(sptree, shared)
  tr$tip.label <- labels[match(tr$tip.label, acc)]
  tr$node.label <- NULL   # ASTRAL support values are not branch labels IQ-TREE should read
  ## Coalescent-unit branch lengths would be meaningless in a likelihood
  ## calculation; drop them so IQ-TREE assigns its own (see run_iqtree_scf).
  tr$edge.length <- NULL
  tr
}
