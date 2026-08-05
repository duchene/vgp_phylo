## Fitted substitution-model parameters, read off the IQ-TREE run that scf_own
## already pays for. Costs no additional IQ-TREE invocations.
##
## WHY THIS EXISTS: it supersedes the `titv` metric. `titv` is a model-free ratio of
## observed transition to transversion distances; `kappa` is the maximum-likelihood
## transition/transversion rate ratio under HKY. Measured on 187 loci they agree at
## Spearman 0.834 -- clearly the same signal -- but they are not interchangeable:
##
##   * Scale. kappa/titv averages 3.00 (sd 0.30, range 2.46-4.20), so titv
##     understates the rate ratio roughly threefold, and by a varying amount. It is
##     not a clean rescaling of kappa (lm R^2 = 0.70).
##   * Saturation. titv correlates -0.76 with mean p-distance against kappa's -0.55,
##     and the kappa/titv gap itself rises monotonically with divergence
##     (2.81 -> 2.94 -> 3.04 -> 3.21 across divergence quartiles, rho = +0.52).
##     That is multiple hits erasing observed transitions -- exactly the bias a
##     likelihood model corrects and a raw distance ratio cannot.
##
## `alpha` is the +G4 gamma shape parameter: low alpha means strong rate
## heterogeneity across sites. It is the ML counterpart of qirp_own's parsimony
## `rate_var`, so the two can be cross-checked.

register_metric(
  name = "subst_model",
  description = "ML substitution-model parameters from the scf_own run: HKY kappa, gamma alpha",
  needs = "tree",
  columns = c("kappa", "gamma_alpha"),
  compute = function(loc) {
    v <- iqtree_locus(loc, loc$treefile, tag = "own", blfix = TRUE,
                      nquartets = 100L, threads = 1L)
    ## Both are dimensionless: kappa is a rate ratio, alpha a shape parameter.
    ## Neither depends on alignment length or taxon count by construction.
    list(
      kappa       = if (is.null(v$kappa)) NA_real_ else v$kappa,
      gamma_alpha = if (is.null(v$alpha)) NA_real_ else v$alpha
    )
  }
)
