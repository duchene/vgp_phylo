## Tree-geometry helpers.

suppressPackageStartupMessages(library(phangorn))

#' Mean tip depth below every internal node of a rooted tree.
#'
#' Returned vector is indexed by internal-node number minus Ntip, i.e.
#' depths[node - Ntip(tr)] for node in (Ntip+1):(Ntip+Nnode).
#'
#' Used as the "time from tip of tree to focal internode" that
#' PhyInformR::Approximator wants. The mean over subtended tips (rather than the
#' max) is the right summary here because these gene trees are not ultrametric:
#' tip depths below one node vary, and QIRP is a property of the typical quartet
#' spanning the internode, not the deepest one.
node_tip_depths <- function(tr) {
  nt <- Ntip(tr)
  d <- ape::dist.nodes(tr)
  tips <- seq_len(nt)
  vapply(seq_len(tr$Nnode) + nt, function(nd) {
    desc <- phangorn::Descendants(tr, nd, "tips")[[1]]
    if (!length(desc)) return(NA_real_)
    mean(d[nd, desc], na.rm = TRUE)
  }, numeric(1))
}

#' Support-weighted Robinson-Foulds discordance: the share of the gene tree's
#' total branch support that sits on splits ABSENT from the reference tree.
#'
#' Plain RF treats a weakly supported conflicting branch exactly like a strongly
#' supported one, which for single-locus gene trees -- where much of the topology
#' is barely resolved -- is the dominant source of noise. Weighting each split by
#' its own support and expressing the result as a proportion of total support
#' gives a value in [0,1] that does not scale with taxon count.
#'
#' @param gtree gene tree with node labels still in IQ-TREE "/support" form (so
#'   read_gene_tree()'s $support is available); tips may carry sequence labels.
#' @param sptree reference tree with bare accession tips.
support_weighted_rf <- function(gtree, sptree) {
  sup <- gtree$support
  if (is.null(sup) || !any(is.finite(sup))) return(NA_real_)

  ## Work on the ORIGINAL gene tree and only translate labels when building split
  ## keys. Any tip-dropping or re-rooting would break the support-to-node mapping,
  ## which is the easiest thing to get silently wrong here.
  nt <- Ntip(gtree)
  if (nt < 4L) return(NA_real_)
  if (length(sup) != gtree$Nnode) return(NA_real_)
  acc <- accession_from_label(gtree$tip.label)
  if (anyDuplicated(acc)) return(NA_real_)
  if (!all(acc %in% sptree$tip.label)) return(NA_real_)

  s <- ape::keep.tip(sptree, acc)

  ## A split is {descendant tips} vs {the rest}. That is orientation-independent,
  ## so keying on whichever side is lexicographically smaller makes gene-tree and
  ## species-tree splits directly comparable without rooting either one.
  canon <- function(side, all_labs) {
    a <- sort(side); b <- sort(setdiff(all_labs, side))
    if (length(a) < 2L || length(b) < 2L) return(NA_character_)   # trivial split
    ka <- paste(a, collapse = ","); kb <- paste(b, collapse = ",")
    if (ka < kb) ka else kb
  }

  ## Species-tree splits, in accession space.
  ss <- phangorn::Descendants(s, (Ntip(s) + 1L):(Ntip(s) + s$Nnode), "tips")
  skeys <- stats::na.omit(vapply(ss, function(ix) canon(s$tip.label[ix], s$tip.label),
                                 character(1)))

  ## Gene-tree splits, node by node, with each node's own support attached.
  gd <- phangorn::Descendants(gtree, (nt + 1L):(nt + gtree$Nnode), "tips")
  gkeys <- vapply(gd, function(ix) canon(acc[ix], acc), character(1))

  ok <- !is.na(gkeys) & is.finite(sup)
  if (!any(ok)) return(NA_real_)
  w <- sup[ok]
  conflict <- !(gkeys[ok] %in% skeys)
  denom <- sum(w)
  if (!is.finite(denom) || denom <= 0) return(NA_real_)
  sum(w[conflict]) / denom
}

#' QIRP/QIPP/QIHP summarised over every internode of a tree, given a site-rate
#' vector. Shared by the qirp_own and qirp_species metrics.
#'
#' All three outputs are probabilities in [0,1], so they are inherently comparable
#' across loci. IQ-TREE's +G4 site rates are scaled to mean ~1, which isolates
#' rate heterogeneity from the locus's overall rate and removes any dependence on
#' alignment length.
## Approximator's cost is O(internodes x sites) and it dominates everything else:
## profiled on a 509-tip / 2500-site locus, Approximator took 21.8s of 22.5s total
## (507 internodes x 2500 rates = 1.27M evaluations), against 0.06s for the site
## rates and 0.18s for the node depths.
##
## Both dimensions are therefore capped by random subsample. This is not only a
## cost fix -- it is the same comparability argument as the fixed quartet sample:
## an uncapped mean would be estimated from 5 internodes at n=8 and 554 at n=557,
## so its sampling precision would scale with tree size. Fixed caps hold effort
## constant. The rate vector only has to represent the rate DISTRIBUTION, so a
## subsample of sites is sufficient for that too.
## Set either to Inf to use everything.
QIRP_MAX_INTERNODES <- 100L
QIRP_MAX_RATES      <- 500L

qirp_summary <- function(tr, rates) {
  if (Ntip(tr) < 4L) stop("need >= 4 tips")
  if (length(rates) > QIRP_MAX_RATES) rates <- sample(rates, QIRP_MAX_RATES)
  rv <- as.matrix(rates)

  ## Approximator needs, per focal internode: t0 = the internode's own length, and
  ## t = time from the tips up to that internode. Gene trees are unrooted, so
  ## midpoint-root first to get a consistent time axis.
  rt <- if (is.rooted(tr)) tr else phangorn::midpoint(tr)
  depths <- node_tip_depths(rt)
  internal <- which(rt$edge[, 2] > Ntip(rt))
  if (!length(internal)) stop("no internal branches")
  n_internal_total <- length(internal)
  if (n_internal_total > QIRP_MAX_INTERNODES)
    internal <- sample(internal, QIRP_MAX_INTERNODES)

  res <- vapply(internal, function(e) {
    child <- rt$edge[e, 2]
    t0 <- rt$edge.length[e]
    t  <- depths[child - Ntip(rt)]
    if (!is.finite(t0) || !is.finite(t) || t0 <= 0 || t <= 0)
      return(c(NA_real_, NA_real_, NA_real_))
    v <- try(suppressWarnings(Approximator(t, t0, rv, 4)), silent = TRUE)
    if (inherits(v, "try-error") || length(v) != 3L)
      return(c(NA_real_, NA_real_, NA_real_))
    as.numeric(v)
  }, numeric(3))

  qirp <- res[1, ]; qipp <- res[2, ]; qihp <- res[3, ]
  ok <- is.finite(qirp)
  if (!any(ok)) stop("Approximator returned no finite QIRP")

  ## (QIRP, QIPP, QIHP) is a COMPOSITION: the three sum to exactly 1 (verified,
  ## max deviation 1.1e-15). So there are only TWO free dimensions, and reporting
  ## all three is redundant by construction -- which is why qirp_mean correlated
  ## -0.86 with qihp_mean. QIPP is therefore dropped as the remainder; QIRP
  ## (predicted power) and QIHP (predicted actively-wrong signal) are kept, being
  ## the two directly interpretable ones.
  m <- mean(qirp[ok])
  v <- if (sum(ok) > 1L) stats::var(qirp[ok]) else NA_real_
  ## Raw variance of a [0,1] quantity is mechanically bounded by m(1-m), so it
  ## tracks the mean for reasons that carry no information: rho(mean, var) = -0.83
  ## raw, falling to -0.54 once divided by that bound. The normalised form is
  ## reported instead -- "how unevenly is predicted power spread across internodes,
  ## relative to what is possible at this mean".
  vmax <- m * (1 - m)
  list(
    qirp_mean           = m,
    qihp_mean           = mean(qihp[ok]),
    qirp_var_rel        = if (is.finite(v) && is.finite(vmax) && vmax > 1e-12) v / vmax else NA_real_,
    n_internodes        = sum(ok),
    n_internodes_total  = n_internal_total,
    ## rates are scaled to mean 1, so their variance IS the rate-heterogeneity of
    ## the locus and is comparable across loci; the mean carries no information.
    rate_var            = stats::var(rates),
    prop_invariant_rate = mean(rates <= 1e-8)
  )
}

#' Relative per-site substitution rates, from per-site parsimony steps on a fixed
#' tree, scaled to mean 1.
#'
#' Why not IQ-TREE's +G4 site rates: that would mean fitting a rate model per
#' locus, which is exactly the cost we are avoiding. Per-site parsimony steps on
#' the already-inferred tree are effectively free (0.01s for 557 tips x 2190
#' sites) and give the quantity QIRP actually needs -- the SHAPE of the rate
#' distribution across sites.
#'
#' Scaling to mean 1 is what makes this comparable across loci: it isolates rate
#' heterogeneity from the locus's overall rate (which enters QIRP separately,
#' through the branch lengths t and t0), and removes any dependence on alignment
#' length.
#'
#' Caveat worth knowing: parsimony undercounts substitutions at fast sites, so the
#' rate distribution is mildly compressed relative to an ML estimate. It is a
#' proxy, chosen for cost.
parsimony_site_rates <- function(tr, al) {
  pd <- phangorn::as.phyDat(al)
  ## as.phyDat collapses the alignment to unique site PATTERNS. Scores come back
  ## per pattern, so expand through the index to recover one value per site --
  ## otherwise the rate distribution would be reweighted by pattern uniqueness
  ## rather than reflecting the alignment's actual site composition.
  sc <- as.numeric(phangorn::parsimony(tr, pd, site = "site"))
  idx <- attr(pd, "index")
  r <- if (!is.null(idx) && length(idx)) sc[idx] else sc
  r <- r[is.finite(r)]
  if (!length(r)) stop("no per-site parsimony scores")
  mu <- mean(r)
  if (!is.finite(mu) || mu <= 0) stop("all sites invariant; no rate variation")
  ## Invariant sites legitimately have rate 0 and are kept: a locus dominated by
  ## conserved sites genuinely carries less resolving power, and QIRP should see
  ## that. A floor is applied only because Approximator works on log scales.
  pmax(r / mu, 1e-8)
}
