## Mean transition/transversion ratio across sequence pairs.
## Wraps ~/scripts/meanTiTv.R.

register_metric(
  name = "titv",
  description = "Model-free Ti/Tv ratio over sequence pairs (see subst_model for ML kappa)",
  needs = "alignment",
  columns = c("titv", "ts_mean", "tv_mean", "pdist_mean"),
  compute = function(loc) {
    al <- loc$alignment
    if (nrow(al) < 2L) stop("need >= 2 sequences")
    ## Computed as mean(TS distance) / mean(TV distance) rather than
    ## mean(TS/TV per pair): the per-pair ratio is undefined whenever a pair has
    ## zero transversions, which is common for closely related tips in short
    ## alignments. Ratio-of-means is defined as long as any pair has a
    ## transversion, and is what ~/scripts/meanTiTv.R does.
    ##
    ## RELATIONSHIP TO HKY kappa (measured on 187 loci): rho = 0.834, so this tracks
    ## the ML rate ratio but is NOT a substitute for it. kappa/titv averages 3.00
    ## (sd 0.30) and the gap widens with divergence (2.81 -> 3.21 across divergence
    ## quartiles), because multiple hits erase observed transitions -- a bias titv
    ## cannot correct. titv also correlates -0.76 with mean p-distance vs kappa's
    ## -0.55. Prefer subst_model$kappa; this metric is retained as a model-free
    ## cross-check for loci where the model fit is doubtful.
    ts <- ape::dist.dna(al, model = "TS", pairwise.deletion = TRUE)
    tv <- ape::dist.dna(al, model = "TV", pairwise.deletion = TRUE)
    pd <- ape::dist.dna(al, model = "raw", pairwise.deletion = TRUE)
    tvm <- mean(tv, na.rm = TRUE)
    list(
      ## A ratio of two per-site rates: dimensionless, and independent of both
      ## alignment length and taxon count.
      titv       = if (is.finite(tvm) && tvm > 0) mean(ts, na.rm = TRUE) / tvm else NA_real_,
      ts_mean    = mean(ts, na.rm = TRUE),
      tv_mean    = tvm,
      pdist_mean = mean(pd, na.rm = TRUE)
    )
  }
)
