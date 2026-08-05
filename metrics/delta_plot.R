## Holland et al. (2002) delta: quartet-level departure from treelikeness,
## i.e. a signal of reticulation / conflicting distance structure.
##
## Holland BR, Huber KT, Dress A, Moulton V (2002) "Delta plots: a tool for
## analyzing phylogenetic distance data." Mol Biol Evol 19:2051-2059.

register_metric(
  name = "delta_plot",
  description = "Holland's delta (reticulation signal) over sampled quartets: mean/var",
  needs = "alignment",
  columns = c("delta_mean", "delta_var", "delta_q90", "q_used", "ntips_used"),
  compute = function(loc) {
    al <- loc$alignment
    if (nrow(al) < 4L) stop("need >= 4 sequences")

    ## Delta is defined on a distance matrix. Using a corrected distance (TN93)
    ## rather than raw p-distance matters here: saturation compresses raw
    ## distances non-linearly, which by itself manufactures four-point violations
    ## and would confound reticulation signal with divergence.
    d <- as.matrix(ape::dist.dna(al, model = "TN93", pairwise.deletion = TRUE))
    if (!any(is.finite(d) & d > 0)) stop("no finite positive distances")
    ## TN93 is undefined for some highly divergent pairs; fall back to raw for
    ## those cells so a handful of NAs does not void the whole locus.
    if (anyNA(d)) {
      draw <- as.matrix(ape::dist.dna(al, model = "raw", pairwise.deletion = TRUE))
      d[is.na(d)] <- draw[is.na(d)]
    }

    vals <- delta_plot_values(d, N_QUARTETS)
    vals <- vals[is.finite(vals)]
    if (!length(vals)) stop("no resolvable quartets")

    list(
      delta_mean = mean(vals),
      delta_var  = if (length(vals) > 1L) stats::var(vals) else NA_real_,
      delta_q90  = as.numeric(stats::quantile(vals, 0.9)),
      q_used     = length(vals),
      ntips_used = nrow(al)
    )
  }
)
