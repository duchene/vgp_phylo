## Gene-tree vs species-tree topological distance, using generalised
## (information-theoretic) split comparisons rather than Robinson-Foulds.
##
## Traditional RF is deliberately EXCLUDED. RF counts splits present in one tree
## and absent from the other, scoring each as a hard match/mismatch. Two properties
## make it a poor choice here:
##   * No partial credit. A single rogue taxon moved across the tree destroys every
##     split it participates in, so RF reports near-maximal distance for trees that
##     are otherwise identical. With single-locus gene trees over up to 557 taxa,
##     that saturation is the normal case rather than the exception.
##   * Coarse and size-bound. RF takes only even integer values up to 2(n-3), so its
##     resolution and its range both depend on taxon count.
## The three measures below all give partial credit for splits that ALMOST match,
## and are normalised against the information content of the trees being compared,
## which is what makes them comparable across loci of 8 to 557 tips.

suppressPackageStartupMessages({
  library(phangorn)
  library(TreeDist)
})

## Of the generalised measures, the Clustering Information Distance is used.
## Rationale: it is a DISTANCE, so 0 means identical trees (Nye similarity, where
## 1 meant identical, is dropped for that reason). Against MSID the empirical
## choice barely matters -- on a 187-locus sample spanning all clades the two
## correlate at Spearman +0.987, and both hit exactly 0 for the same 5 loci -- but
## CID uses more of the available range here (0-0.305 vs 0-0.188), giving ~1.6x the
## dynamic range for downstream modelling.

register_metric(
  name = "treedist_species",
  description = "Gene-vs-species tree distance: clustering information distance + support-weighted",
  needs = c("tree", "species_tree"),
  species_tree = "astralpro",
  columns = c("cid", "discord_support_weighted", "ntips_shared"),
  compute = function(loc) {
    gt <- to_species_labels(loc$tree)
    sp <- loc$species_tree

    shared <- intersect(gt$tip.label, sp$tip.label)
    if (length(shared) < 4L) stop("only ", length(shared), " shared taxa; need >= 4")
    g <- ape::unroot(ape::keep.tip(gt, shared))
    s <- ape::unroot(ape::keep.tip(sp, shared))

    ## normalize = TRUE scales by the information content of the trees compared, so
    ## 0 = identical and 1 = as different as trees of this size can be. That is what
    ## makes the value valid across 8-557 tips.
    cid <- as.numeric(TreeDist::ClusteringInfoDistance(g, s, normalize = TRUE))

    ## Share of the gene tree's total branch support sitting on splits absent from
    ## the species tree. Weighting by support rather than counting splits means a
    ## poorly supported conflict is penalised less than a well supported one -- the
    ## dominant noise source in single-locus gene trees.
    ##
    ## No branch-length weighted distance is reported: the gene trees are in
    ## substitutions/site while the ASTRAL species tree is in coalescent units, so
    ## any length-weighted comparison between them mixes incommensurable units.
    list(
      cid                      = cid,
      discord_support_weighted = support_weighted_rf(loc$tree, sp),
      ntips_shared             = length(shared)
    )
  }
)
