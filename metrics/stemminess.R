## Stemminess: share of total tree length carried by internal branches.
## Uses stemmy(), vendored in scripts/lib/vendored.R from ~/scripts/stemmy.R.


register_metric(
  name = "stemminess",
  description = "Proportion of tree length on internal branches (+ size-corrected form)",
  needs = "tree",
  columns = c("stemminess", "stemminess_per_branch", "n_internal", "n_branches"),
  compute = function(loc) {
    tr <- loc$tree
    if (Ntip(tr) < 4L) stop("need >= 4 tips")
    internal <- tr$edge[, 2] > Ntip(tr)
    s <- stemmy(tr)
    ## Raw stemminess is a proportion, so it is already free of rate and length
    ## effects -- but it is NOT free of taxon-count effects: the share of edges
    ## that are internal rises from 5/13 at n=8 towards 1/2 as n grows, so bigger
    ## trees drift upwards for purely combinatorial reasons. Dividing by the
    ## proportion of branches that are internal removes that trend, giving the
    ## mean internal branch length relative to the mean branch length.
    prop_internal_edges <- sum(internal) / length(internal)
    list(
      stemminess            = s,
      stemminess_per_branch = s / prop_internal_edges,
      n_internal            = sum(internal),
      n_branches            = length(internal)
    )
  }
)
