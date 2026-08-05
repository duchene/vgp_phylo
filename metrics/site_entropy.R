## Site-pattern information: per-site nucleotide entropy, plus the unconstrained
## (multinomial) log-likelihood per site.

register_metric(
  name = "site_entropy",
  description = "Per-site entropy (bias-corrected) + unconstrained lnL per site",
  needs = "alignment",
  columns = c("entropy_mean", "entropy_var", "entropy_mean_mm", "entropy_var_mm",
              "entropy_mean_norm", "prop_invariant",
              "unc_lnl_per_site", "n_patterns", "nsites_used", "ntips_used"),
  compute = function(loc) {
    mat <- tolower(as.character(loc$alignment))
    n_all <- nrow(mat)
    if (n_all < 2L) stop("need >= 2 sequences")

    ## Count only unambiguous ACGT per site; gaps/N/IUPAC codes are missing data.
    counts <- apply(mat, 2, function(col) {
      tabulate(match(col, c("a", "c", "g", "t")), nbins = 4L)
    })                                   # 4 x nsites
    n_obs <- colSums(counts)
    keep <- n_obs >= 2L                  # entropy is undefined/0 with < 2 observations
    if (!any(keep)) stop("no site has >= 2 unambiguous bases")
    counts <- counts[, keep, drop = FALSE]
    n_obs <- n_obs[keep]

    p <- sweep(counts, 2, n_obs, "/")
    plogp <- ifelse(p > 0, p * log(p), 0)
    H <- -colSums(plogp)                 # natural-log entropy, in [0, log 4]

    ## Plug-in entropy is biased DOWNWARDS by roughly (k-1)/(2n), where k is the
    ## number of observed states and n the number of sequences at that site. Since
    ## n varies from 8 to 557 across these loci, that bias alone would create a
    ## spurious correlation between entropy and locus size. The Miller-Madow
    ## correction removes the leading term and is what makes this comparable.
    k_obs <- colSums(counts > 0)
    H_mm <- H + (k_obs - 1) / (2 * n_obs)

    ## Unconstrained (multinomial) log-likelihood: sum over distinct site patterns
    ## of n_i*log(n_i), minus nsites*log(nsites). This is the likelihood ceiling no
    ## substitution model can exceed, and equivalently a measure of how repetitive
    ## the site patterns are.
    ##
    ## This is the one statistic retained from PhyloMAd. It is the algebraic
    ## equivalent of ~/phylomad/codeFolder/testStatistics/get.unconstrained.R --
    ## verified identical (all.equal TRUE) on a 509x2500 locus -- reimplemented here
    ## so this metric carries no dependency on the phylomad checkout, and running
    ## ~1.6x faster. The other PhyloMAd statistics were dropped as redundant: their
    ## rescaled forms correlated -1.00 with stemminess and +0.99 with
    ## aitchison_basefreq, whereas this one correlates only +0.22 with entropy.
    ##
    ## Computed on the FULL alignment (all columns, gaps and N included as states),
    ## matching phylomad, not on the ACGT-filtered subset used for entropy above.
    pat <- do.call(paste0, asplit(mat, 1))
    tb <- tabulate(match(pat, unique(pat)))
    unc <- sum(tb * log(tb)) - length(pat) * log(length(pat))

    list(
      entropy_mean      = mean(H),
      entropy_var       = if (length(H) > 1L) stats::var(H) else NA_real_,
      entropy_mean_mm   = mean(H_mm),
      entropy_var_mm    = if (length(H_mm) > 1L) stats::var(H_mm) else NA_real_,
      ## Per-site, because the raw sum scales with alignment length.
      unc_lnl_per_site  = unc / length(pat),
      ## Raw diagnostic: distinct site patterns. NOT comparable across loci -- it
      ## rises with both taxon count and alignment length.
      n_patterns        = length(tb),
      ## Scaled to [0,1] by the maximum possible entropy for 4 states, which makes
      ## the value interpretable as "fraction of maximal site disorder".
      entropy_mean_norm = mean(H) / log(4),
      prop_invariant    = mean(k_obs == 1L),
      nsites_used       = ncol(counts),
      ntips_used        = n_all
    )
  }
)
