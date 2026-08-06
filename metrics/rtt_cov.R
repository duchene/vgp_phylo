## Coefficient of variation of root-to-tip distances in the locus's own tree.
## Uses pathnode(), vendored in scripts/lib/vendored.R from ~/scripts/pathnode.R.

suppressPackageStartupMessages(library(phangorn))

register_metric(
  name = "rtt_cov",
  description = "CoV of root-to-tip distances (midpoint-rooted own tree)",
  needs = "tree",
  columns = c("rtt_cov", "rtt_mean", "ntips_used"),
  compute = function(loc) {
    tr <- loc$tree
    ## Gene trees are unrooted; rtt.cov midpoint-roots internally. Midpoint
    ## rooting is a property of the tree alone, so it stays comparable across
    ## loci with different taxon sets (no outgroup is required or assumed).
    if (Ntip(tr) < 3L) stop("need >= 3 tips")
    rooted <- if (is.rooted(tr)) tr else phangorn::midpoint(tr)
    rtts <- pathnode(rooted)[[1]]
    list(
      rtt_cov    = stats::sd(rtts, na.rm = TRUE) / mean(rtts, na.rm = TRUE),
      rtt_mean   = mean(rtts, na.rm = TRUE),
      ntips_used = Ntip(tr)
    )
  }
)
