## Site concordance factors of the alignment against the locus's OWN gene tree.
## Measures how much of the site-level signal supports each branch the gene tree
## itself asserts -- i.e. internal consistency of the locus.

SCF_QUARTETS <- 100L

register_metric(
  name = "scf_own",
  description = "sCF vs own gene tree (IQ-TREE --scfl): mean/var over branches",
  needs = "tree",
  columns = c("scf_mean", "scf_var", "sdf1_mean", "sdf2_mean",
              "sn_mean", "n_branches", "scf_prop_sig"),
  compute = function(loc) {
    ## The gene tree is used exactly as inferred: blfix = TRUE holds its branch
    ## lengths, and the model has no free parameters, so IQ-TREE estimates nothing.
    tab <- iqtree_locus(loc, loc$treefile, tag = "own", blfix = TRUE,
                        nquartets = SCF_QUARTETS, threads = 1L)$cf
    if (!nrow(tab)) stop("no branch had informative sites")

    ## Binomial test per branch against the null that all three quartet
    ## resolutions are equally supported (p = 1/3), following
    ## ~/scripts/binom.scf.test.R. Reported as a PROPORTION of branches, not a
    ## count, so it does not scale with tree size.
    sig <- mapply(function(scf, sn) {
      if (!is.finite(sn) || sn < 1) return(NA)
      stats::binom.test(round(scf * sn), round(sn), p = 1/3,
                        alternative = "greater")$p.value < 0.001
    }, tab$sCF, tab$sN)

    list(
      ## sCF is a proportion of informative sites, and --scfl samples a fixed
      ## number of quartets per branch, so both the quantity and its sampling
      ## effort are comparable across loci of any size.
      scf_mean     = mean(tab$sCF),
      scf_var      = if (nrow(tab) > 1L) stats::var(tab$sCF) else NA_real_,
      sdf1_mean    = mean(tab$sDF1),
      sdf2_mean    = mean(tab$sDF2),
      sn_mean      = mean(tab$sN),
      n_branches   = nrow(tab),
      scf_prop_sig = mean(sig, na.rm = TRUE)
    )
  }
)
