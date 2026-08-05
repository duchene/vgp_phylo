## Quartet Internode Resolution Probability (Townsend, Su & Tekle 2012), per
## internode of the locus's OWN gene tree, via PhyInformR::Approximator.
##
## QIRP asks: given this locus's distribution of site rates, what is the
## probability that a quartet spanning this internode is resolved CORRECTLY?
## Approximator() returns the triple (QIRP, QIPP, QIHP) = probability correct,
## polytomous, and homoplasious (incorrectly resolved).
##
## Nothing is re-inferred: the internode lengths (t0) and tip-to-node depths (t)
## are read straight off the existing IQ-TREE gene tree, and the site-rate vector
## comes from per-site parsimony on that same fixed tree.

suppressWarnings(suppressPackageStartupMessages(library(PhyInformR)))

register_metric(
  name = "qirp_own",
  description = "PhyInformR QIRP + QIHP per internode of own tree (compositional: 2 free dims)",
  needs = c("tree", "alignment"),
  columns = c("qirp_mean", "qihp_mean", "qirp_var_rel",
              "n_internodes", "n_internodes_total", "rate_var", "prop_invariant_rate"),
  compute = function(loc) {
    tr <- loc$tree
    if (Ntip(tr) < 4L) stop("need >= 4 tips")
    rates <- parsimony_site_rates(tr, loc$alignment)
    qirp_summary(tr, rates)
  }
)
