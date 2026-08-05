## Central paths and defaults. Sourced by every script; edit here, not in scripts.

VGP_ROOT <- Sys.getenv("VGP_PHYLO_ROOT", unset = "/home/davidd/vgp_phylo")

PATHS <- list(
  root          = VGP_ROOT,
  loci          = file.path(VGP_ROOT, "data", "loci"),
  gene_trees    = file.path(VGP_ROOT, "data", "gene_trees"),
  species_trees = file.path(VGP_ROOT, "data", "species_trees"),
  metadata      = file.path(VGP_ROOT, "data", "metadata"),
  manifest      = file.path(VGP_ROOT, "data", "manifest.tsv"),
  results       = file.path(VGP_ROOT, "results", "metrics"),
  ## Cached IQ-TREE per-branch concordance tables. IQ-TREE dominates the runtime
  ## of this analysis, so several metrics read one cached table per locus rather
  ## than each re-running it.
  scf_cache     = file.path(VGP_ROOT, "results", "scf_tables"),
  merged        = file.path(VGP_ROOT, "results", "locus_metrics.tsv"),
  logs          = file.path(VGP_ROOT, "logs"),
  figures       = file.path(VGP_ROOT, "figures")
)

## Where metric definitions are looked up, in order. The project directory holds
## metrics specific to this dataset; the shared directory holds ones reused across
## projects. First match wins, so a project metric can shadow a shared one.
METRIC_PATH <- c(
  file.path(VGP_ROOT, "metrics"),
  Sys.getenv("PHYLO_METRICS_SHARED", unset = "/home/davidd/scripts/metrics")
)

## Reusable single-function scripts (~/scripts), sourced by name via
## source_shared(). Kept separate from METRIC_PATH: these are helpers, not metrics.
SHARED_SCRIPTS <- Sys.getenv("PHYLO_SCRIPTS_SHARED", unset = "/home/davidd/scripts")

## External binaries.
IQTREE <- Sys.getenv("IQTREE", unset = "/home/davidd/iqtree3")

## --scfl is likelihood-based, so it REQUIRES a substitution model. The model is
## held fixed across all loci rather than selected per locus: otherwise sCF would
## partly reflect which model happened to win, and nothing would be comparable.
##
## HKY+G4 is the chosen compromise, measured on 187 loci spanning all clades:
##
##   model    wall (187 loci, 12 cores)   vs GTR+G4: rho / mean diff / max diff
##   GTR+G4   49 s                        --
##   HKY+G4   26 s                        0.9973 / +0.0003 / 0.0139
##   JC       12 s                        0.9896 / +0.0019 / 0.0371
##
## HKY+G4 is statistically indistinguishable from the richest model for this
## purpose at half the cost, while still modelling the two things that demonstrably
## vary in this dataset: transition/transversion bias (the titv metric spans
## 0.69-1.45) and across-site rate heterogeneity. JC is cheaper again and preserves
## locus RANKING well (rho 0.99 on the mean, 0.95 on the variance), so it is
## defensible if cost becomes binding -- but it assumes equal base frequencies,
## which aitchison_basefreq shows is false here, and pairs poorly with branch
## lengths inferred under a richer model.
##
## Only model parameters are estimated; branch lengths stay fixed via -blfix.
## The model is part of the IQ-TREE cache key, so changing it invalidates cleanly.
IQTREE_MODEL <- Sys.getenv("IQTREE_MODEL", unset = "HKY+G4")

## Reference species trees, by short name. Metrics that compare a gene tree to a
## species tree take one of these keys.
SPECIES_TREES <- c(
  astralpro = file.path(PATHS$species_trees, "astralpro.l6.tre"),
  wastral   = file.path(PATHS$species_trees, "wastral.merged.nwk")
)

## Tip-label aliases applied when a species tree is loaded. Empty: the one
## historical mismatch (human in astralpro.l6.tre) was fixed in the file itself
## -- see README.md, "Species tree tip label fix". Kept as a hook in case another
## label discrepancy turns up.
TIP_ALIASES <- character(0)

## This machine is shared -- other jobs run here. Keep defaults well below nproc.
DEFAULTS <- list(
  cores            = 8L,    # cores used *within* one metric
  metric_jobs      = 4L,    # metrics run concurrently (total load ~ cores*jobs)
  block            = 2000L  # loci per output flush; also the resume granularity
)

PHYLOMAD_DIR <- Sys.getenv("PHYLOMAD_DIR", unset = "/home/davidd/phylomad")
