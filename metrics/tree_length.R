## Overall evolutionary rate of the locus, taken from the existing IQ-TREE gene
## tree rather than re-estimated.
##
## COMPARABILITY WARNING, and the reason this metric reports four columns instead
## of one: TOTAL tree length is NOT comparable across loci. An unrooted binary tree
## on n tips has 2n-3 branches, so total length grows roughly linearly with taxon
## count -- a 557-tip tree will out-total an 8-tip tree at any rate. Use
## `mean_branch_length` (total / number of branches) as the rate proxy, or
## `mean_root_to_tip` if a per-lineage depth is wanted. `total_length` is reported
## only so the normalisation can be audited.

source_shared("pathnode")

register_metric(
  name = "tree_length",
  description = "Overall rate proxy from the existing gene tree (per-branch and depth)",
  needs = "tree",
  columns = c("mean_branch_length", "median_branch_length", "mean_root_to_tip",
              "total_length", "n_branches", "prop_zero_length"),
  compute = function(loc) {
    tr <- loc$tree
    bl <- tr$edge.length
    if (is.null(bl) || !any(is.finite(bl))) stop("tree has no branch lengths")
    bl <- bl[is.finite(bl)]

    ## Mean tip depth on the midpoint-rooted tree: a per-lineage measure of
    ## accumulated substitutions, so it does not grow with taxon count the way a
    ## sum over branches does.
    rooted <- if (is.rooted(tr)) tr else phangorn::midpoint(tr)
    rtt <- pathnode(rooted)[[1]]

    list(
      mean_branch_length   = mean(bl),
      median_branch_length = stats::median(bl),
      mean_root_to_tip     = mean(rtt, na.rm = TRUE),
      total_length         = sum(bl),
      n_branches           = length(bl),
      ## IQ-TREE collapses unresolvable branches to ~0; the share of these says how
      ## much of the topology is effectively unresolved, and is a proportion.
      prop_zero_length     = mean(bl <= 1e-8)
    )
  }
)
