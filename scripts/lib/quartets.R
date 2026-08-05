## Quartet sampling and compositional-distance helpers.
##
## Comparability note that drives the whole file: the number of quartets in a
## tree is choose(n, 4), which runs from 70 at n=8 to ~4e9 at n=557. Any metric
## that averages over "all quartets" would therefore be estimated with wildly
## different precision across loci, and would be uncomputable for the big ones.
## Every quartet-based metric here instead draws a FIXED number of quartets
## (N_QUARTETS), so sampling effort -- and hence the sampling variance of the
## mean -- is constant across loci regardless of size.

N_QUARTETS <- 2000L

#' Draw n random 4-subsets of `tips` (indices), without replacement within a
#' quartet. Returns an n x 4 integer matrix, or NULL if fewer than 4 tips.
#' For small trees the number of distinct quartets may be < n; we sample with
#' replacement across quartets, which keeps the estimator unbiased and the
#' effort constant.
sample_quartets <- function(ntips, n = N_QUARTETS) {
  if (ntips < 4L) return(NULL)
  t(vapply(seq_len(n), function(i) sample.int(ntips, 4L), integer(4)))
}

# ------------------------------------------------------- Holland's delta ----

#' Holland et al. (2002) delta for one quartet, from a patristic/pairwise
#' distance matrix. For tips i,j,k,l the three pairings give sums
#'   s1 = d(i,j)+d(k,l),  s2 = d(i,k)+d(j,l),  s3 = d(i,l)+d(j,k)
#' Sorted descending as m1 >= m2 >= m3, delta = (m1 - m2) / (m1 - m3).
#' delta = 0 for a perfectly treelike quartet, 1 for a fully "boxy" (reticulate)
#' one. Already scale-free -- it is a ratio of distance differences, so it is
#' invariant to overall rate and to alignment length.
delta_quartet <- function(d, q) {
  i <- q[1]; j <- q[2]; k <- q[3]; l <- q[4]
  s <- c(d[i, j] + d[k, l], d[i, k] + d[j, l], d[i, l] + d[j, k])
  s <- sort(s, decreasing = TRUE)
  den <- s[1] - s[3]
  if (!is.finite(den) || den <= 0) return(NA_real_)   # degenerate: all pairings equal
  (s[1] - s[2]) / den
}

#' delta over a fixed sample of quartets. Returns the raw vector so callers can
#' take whatever summary they need.
delta_plot_values <- function(d, nq = N_QUARTETS) {
  n <- nrow(d)
  qs <- sample_quartets(n, nq)
  if (is.null(qs)) return(numeric(0))
  vapply(seq_len(nrow(qs)), function(r) delta_quartet(d, qs[r, ]), numeric(1))
}

# ------------------------------------------------- compositional distance ----

#' Mean number of distinct unambiguous bases per site, rarefied to a fixed
#' subsample of m sequences.
#'
#' The unrarefied version (PhyloMAd's get.biodivstat) rises with taxon sampling
#' for purely combinatorial reasons: with 557 sequences a site is far more likely
#' to show all four bases than with 8, whatever the underlying process. Rarefying
#' to the smallest taxon count in the dataset (m = 8) makes the statistic a
#' property of the sequences rather than of how many were sampled.
#'
#' @param mat lowercase character matrix (taxa x sites)
#' @param m   subsample size; loci with fewer than m sequences return NA
#' @param reps number of independent subsamples to average over
rarefied_state_richness <- function(mat, m = 8L, reps = 20L) {
  n <- nrow(mat)
  if (n < m) return(NA_real_)
  vals <- vapply(seq_len(reps), function(i) {
    sub <- mat[sample.int(n, m), , drop = FALSE]
    k <- apply(sub, 2, function(col) length(unique(col[col %in% c("a", "c", "g", "t")])))
    mean(k[k > 0])
  }, numeric(1))
  mean(vals, na.rm = TRUE)
}

#' Turn observed counts into a composition with exact zeros smoothed away.
#'
#' Aitchison geometry is undefined on zeros, and concordance-factor compositions
#' contain them constantly (a branch with zero discordant sites). Multiplicative
#' replacement based on "a fraction of the smallest positive part" is the usual
#' recipe but misbehaves badly here: for the composition (1, 0, 0) the smallest
#' positive part is 1, so the replacement value would be enormous.
#'
#' Working from counts instead gives a principled amount of smoothing: add a
#' Jeffreys pseudocount of 1/2 to each part. A branch backed by 500 informative
#' sites is then smoothed far less than one backed by 3, which is exactly the
#' desired behaviour and keeps the transform stable across loci of any depth.
#' (Counts here are averages over sampled quartets, so they need not be integers.)
comp_from_counts <- function(counts, prior = 0.5) {
  x <- as.numeric(counts) + prior
  x / sum(x)
}

#' Centred log-ratio transform.
clr <- function(x) { lx <- log(x); lx - mean(lx) }

#' Aitchison distance between two compositions: Euclidean distance in clr space.
#' Both arguments must already be strictly positive (see comp_from_counts).
#' Scale-invariant -- it depends only on ratios within a composition -- so it is
#' directly comparable across loci of different depth, length, or taxon count.
aitchison_dist <- function(x, y) {
  x <- as.numeric(x) / sum(x)
  y <- as.numeric(y) / sum(y)
  if (any(x <= 0) || any(y <= 0)) return(NA_real_)
  sqrt(sum((clr(x) - clr(y))^2))
}
