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

## Where metric definitions are looked up. Only the in-repo directory by default,
## so a checkout behaves identically on any machine. Set PHYLO_METRICS_SHARED to
## prepend an external directory if you want to develop a metric outside the repo;
## anything relied on for published results belongs in metrics/.
METRIC_PATH <- c(
  file.path(VGP_ROOT, "metrics"),
  Sys.getenv("PHYLO_METRICS_SHARED", unset = NA_character_)
)
METRIC_PATH <- METRIC_PATH[!is.na(METRIC_PATH) & nzchar(METRIC_PATH)]

## NOTE: this project is self-contained. Helper functions that used to be sourced
## from ~/scripts at run time are now vendored verbatim in scripts/lib/vendored.R,
## with their upstream md5s recorded, so results never depend on an unversioned
## directory outside the repository.

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

## Core budget. Two things make the naive reading of `nproc` wrong here:
##
## 1. nproc reports 64, but this is a 32-PHYSICAL-core Threadripper 5975WX with
##    2-way SMT. Compute-bound work gets little from an SMT sibling, so the real
##    ceiling is 32, not 64.
## 2. IQ-TREE likelihood is memory-bandwidth-bound. A 156-tip x 2500-site HKY+G4
##    partial-likelihood array is ~48 MB, which already exceeds one 32 MiB L3
##    instance (128 MiB total, 4 instances, 8 physical cores each). Concurrent jobs
##    therefore thrash L3 and stream from DRAM.
##
## Measured consequence: scaling 12 -> 32 workers multiplied throughput by only
## 1.16x (size-normalised), i.e. per-core efficiency fell 57%. The same locus took
## 1.0s at 12-core load and 5.3s at 32-core load. Do not assume more cores helps.
DEFAULTS <- list(
  cores            = 8L,    # cores used *within* one metric
  metric_jobs      = 2L,    # metrics run concurrently (total load ~ cores*jobs)
  block            = 2000L  # loci per output flush; also the resume granularity
)

## (PHYLOMAD_DIR removed: the phylomad_adequacy metric was dropped as redundant, so
## nothing here reads ~/phylomad any more. See metrics/README.md.)
