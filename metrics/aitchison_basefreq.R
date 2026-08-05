## Aitchison (compositional) distances between the BASE COMPOSITIONS of the
## sequences in a locus.
##
## Each sequence contributes a 4-part composition (A, C, G, T proportions). These
## are proportions constrained to sum to 1, so they live on the simplex and
## ordinary Euclidean distance between them is not meaningful -- differences in one
## base are necessarily offset by the others. The correct geometry is Aitchison's:
## Euclidean distance between centred log-ratio transforms, which depends only on
## the RATIOS between bases.
##
## What this measures: compositional heterogeneity among taxa. It is the same
## signal the chi-squared composition test targets (phylomad's get.chisqstat, also
## reported by phylomad_adequacy as Cramer's V) but in proper compositional
## geometry, and summarised per locus as a mean and variance over sequence pairs.
##
## Comparability: Aitchison distance is scale-invariant and dimensionless, so it
## does not depend on alignment length or on overall divergence. Two second-order
## dependencies are worth knowing:
##   * Sampling noise in each composition falls as 1/sqrt(n_sites), so very short
##     alignments inflate distances slightly. With 1200-2500 sites here the effect
##     is small, and `nsites_used` is reported so it can be checked or regressed out.
##   * The number of pairs is choose(n,2), so the PRECISION of the mean varies with
##     taxon count even though its expectation does not. `ntips_used` is reported
##     for the same reason.

register_metric(
  name = "aitchison_basefreq",
  description = "Aitchison distance between per-sequence base compositions: mean/var",
  needs = "alignment",
  columns = c("ait_mean", "ait_var", "ait_max",
              "ait_totvar", "ait_to_centroid_mean",
              "gc_mean", "gc_sd", "nsites_used", "ntips_used"),
  compute = function(loc) {
    mat <- tolower(as.character(loc$alignment))
    n <- nrow(mat)
    if (n < 2L) stop("need >= 2 sequences")

    ## Per-sequence ACGT counts; gaps, N and IUPAC ambiguity codes are missing data
    ## and simply do not contribute to the composition.
    cnt <- t(apply(mat, 1, function(r) tabulate(match(r, c("a", "c", "g", "t")), nbins = 4L)))
    tot <- rowSums(cnt)
    keep <- tot >= 20L                 # a composition from <20 bases is noise
    if (sum(keep) < 2L) stop("fewer than 2 sequences with >= 20 unambiguous bases")
    cnt <- cnt[keep, , drop = FALSE]

    ## Jeffreys pseudocount, so a base absent from one sequence does not make the
    ## log-ratio infinite. Smoothing scales with the evidence available: a sequence
    ## with 2500 bases is barely perturbed, a short one is smoothed more.
    comp <- t(apply(cnt, 1, comp_from_counts))

    Z <- t(apply(comp, 1, clr))        # rows = sequences, cols = clr coordinates
    d <- stats::dist(Z)                # Euclidean in clr space == Aitchison
    dv <- as.numeric(d)
    dv <- dv[is.finite(dv)]
    if (!length(dv)) stop("no finite Aitchison distances")

    ## Compositional centre (clr-space centroid) and total variance -- the standard
    ## one-number summary of dispersion for compositional data, equal to the mean
    ## squared Aitchison distance to the centre. Independent of taxon count.
    ctr <- colMeans(Z)
    d2c <- sqrt(rowSums(sweep(Z, 2, ctr)^2))

    gc <- (cnt[, 2] + cnt[, 3]) / rowSums(cnt)

    list(
      ait_mean             = mean(dv),
      ait_var              = if (length(dv) > 1L) stats::var(dv) else NA_real_,
      ait_max              = max(dv),
      ait_totvar           = mean(d2c^2),
      ait_to_centroid_mean = mean(d2c),
      gc_mean              = mean(gc),
      gc_sd                = if (length(gc) > 1L) stats::sd(gc) else NA_real_,
      nsites_used          = ncol(mat),
      ntips_used           = nrow(cnt)
    )
  }
)
