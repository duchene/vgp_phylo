## Shared data layer for the metric <-> model link analyses
## (scripts/make_link_figures.R, make_hub_landscape.R, make_path_analysis.R).
##
## load_link_data() returns one row per locus: the 15 locus metrics, the best-model
## outputs, ancestral positions, and the biological covariates. Every analysis
## variable is already transformed (logs where the scale is multiplicative) and
## named by a short key; VARS gives its label and block.
##
## Four locus metrics are the same quantity as a ModelFinder output and are NOT in
## VARS: kappa (HKY) ~ model ti/tv (rho 0.85 after clade + size), gamma alpha ~
## FreeRate CV (-0.95), mean branch length ~ ModelFinder's (0.99), alignment GC ~
## model GC (0.98). Keeping both would make any cross-block analysis rediscover
## them. The ModelFinder versions are kept; they come from each locus's best model.

VARS <- rbind(
  ## locus metrics (block "metric")
  data.frame(key = "scf_sig",  label = "sCF signif. branches",     block = "metric"),
  data.frame(key = "scf",      label = "sCF mean",                 block = "metric"),
  data.frame(key = "sn",       label = "informative sites/branch", block = "metric"),
  data.frame(key = "qirp",     label = "QIRP (P correct)",         block = "metric"),
  data.frame(key = "qihp",     label = "QIHP (P incorrect)",       block = "metric"),
  data.frame(key = "delta",    label = "Holland delta",            block = "metric"),
  data.frame(key = "cid",      label = "CID vs species tree",      block = "metric"),
  data.frame(key = "comphet",  label = "compositional het.",       block = "metric"),
  data.frame(key = "rttcov",   label = "root-to-tip CoV",          block = "metric"),
  data.frame(key = "entropy",  label = "site entropy",             block = "metric"),
  data.frame(key = "stem",     label = "stemminess",               block = "metric"),
  data.frame(key = "unclnl",   label = "uncon. lnL/site",          block = "metric"),
  ## best-model outputs (block "model")
  data.frame(key = "gc",       label = "GC",                       block = "model"),
  data.frame(key = "titv",     label = "ti/tv",                    block = "model"),
  data.frame(key = "rate_cv",  label = "FreeRate CV",              block = "model"),
  data.frame(key = "het_gain", label = "rate-het. dBIC/site",      block = "model"),
  data.frame(key = "fast_rate", label = "fastest-class rate",      block = "model"),
  data.frame(key = "k",        label = "rate classes k",           block = "model"),
  data.frame(key = "rate",     label = "rate per branch",          block = "model"),
  stringsAsFactors = FALSE)
VLAB <- setNames(VARS$label, VARS$key)

## Six lineage groups, as in scripts/make_figures.R and pca_scores.R.
lineage_of <- function(cl) {
  ifelse(cl == "BirdsAnc0", "Birds",
  ifelse(cl == "SquamataAnc0", "Squamates",
  ifelse(cl == "CartilaginousFishesAnc0", "Cartilaginous fishes",
  ifelse(grepl("^RayFinnedFishes", cl), "Ray-finned fishes",
  ifelse(cl %in% c("AmphibiansAnc0", "AnuraAnc0", "AnuraAnc1", "GymnophionaAnc0"), "Amphibians",
  ifelse(grepl("^Vertebrates", cl), "Deep vertebrates", "Jawless fishes"))))))
}
LINEAGES <- c("Birds", "Squamates", "Cartilaginous fishes", "Amphibians", "Deep vertebrates",
              "Ray-finned fishes")

load_link_data <- function(root = PATHS$root) {
  rd <- function(...) utils::read.delim(file.path(root, ...), stringsAsFactors = FALSE)
  src <- file.path(root, "results", "published", "locus_metrics.tsv.gz")
  lm_ <- utils::read.delim(gzfile(src), stringsAsFactors = FALSE)
  mp  <- rd("results", "modelfinder", "model_params.tsv")
  pos <- rd("results", "stage0", "locus_position.tsv")
  chr <- rd("results", "stage0", "locus_chromosome.tsv")
  ## heterogeneity gain: BIC of GTR+F with one rate class minus the best model's
  ## BIC, per site. Read only those rows from the 141 MB all-models table.
  am <- utils::read.delim(pipe(sprintf("zcat '%s' | awk -F'\\t' 'NR == 1 || $2 == \"GTR+F\"'",
                          file.path(root, "results", "modelfinder", "all_models.tsv.gz"))),
                          stringsAsFactors = FALSE)
  i <- match(lm_$locus_id, mp$locus_id); stopifnot(!anyNA(i))
  mp <- mp[i, ]
  d <- data.frame(locus_id = lm_$locus_id, clade = lm_$clade, lineage = lineage_of(lm_$clade),
                  ntips = mp$ntips_kept, nsites = lm_$nsites, subst = mp$subst,
                  stringsAsFactors = FALSE)
  d$scf_sig  <- lm_$scf_own.scf_prop_sig
  d$scf      <- lm_$scf_own.scf_mean
  d$sn       <- lm_$scf_own.sn_mean
  d$qirp     <- lm_$qirp_own.qirp_mean
  d$qihp     <- lm_$qirp_own.qihp_mean
  d$delta    <- lm_$delta_plot.delta_mean
  d$cid      <- lm_$treedist_species.cid
  d$comphet  <- log10(lm_$aitchison_basefreq.ait_totvar + 1e-4)
  d$rttcov   <- lm_$rtt_cov.rtt_cov
  d$entropy  <- lm_$site_entropy.entropy_mean_mm
  d$stem     <- lm_$stemminess.stemminess_per_branch
  d$unclnl   <- lm_$site_entropy.unc_lnl_per_site
  d$gc       <- mp$gc
  d$titv     <- log10(mp$titv)
  d$rate_cv  <- mp$rate_cv
  bic_e      <- am$BIC[match(d$locus_id, am$locus_id)]
  d$het_gain <- log10(pmax((bic_e - mp$BIC) / mp$nsites, 0) + 1e-3)
  d$fast_rate <- log10(mp$r_fast_rate)
  d$k        <- mp$k
  d$rate     <- log10(mp$mean_branch_length)
  ## positions and biological covariates
  ip <- match(d$locus_id, pos$locus_id)
  d$contig <- pos$contig[ip]; d$mid <- pos$mid[ip]
  d$n_loci_contig <- pos$n_loci_contig[ip]; d$dist_to_outer <- pos$dist_to_outer[ip]
  d$log_dist_end <- log10(pos$dist_to_outer[ip] + 1e3)   # 1 kb floor: outermost loci sit at ~1 kb
  d$log_contig_span <- log10(pos$span_end[ip] - pos$span_start[ip] + 1e3)   # contig size proxy
  d$gc_align <- lm_$aitchison_basefreq.gc_mean
  ic <- match(d$locus_id, chr$locus_id)
  d$chicken_chr <- chr$chicken_chr[ic]; d$chicken_class <- chr$chicken_class[ic]
  d$log_chicken_chr_len <- log10(chr$chicken_chr_len[ic])
  d$human_chr <- chr$human_chr[ic]
  d$log_human_chr_len <- log10(chr$human_chr_len[ic])
  sm <- file.path(root, "results", "stage0", "locus_softmask.tsv")
  if (file.exists(sm)) {
    s <- utils::read.delim(sm, stringsAsFactors = FALSE)
    d$softmask <- s$softmask[match(d$locus_id, s$locus_id)]
  } else d$softmask <- NA_real_
  d
}

## Within-clade, size-free version of a variable: residual of its rank on clade,
## log tips and log sites (as the partial panel of figure A2), rescaled to unit SD.
## Ranks keep a few extreme loci from steering the linear models built on top.
## With a single clade in `d` the clade term drops out.
resid_clade_size <- function(x, d, rank = TRUE) {
  r <- rep(NA_real_, length(x)); ok <- is.finite(x)
  y <- if (rank) base::rank(x[ok]) else x[ok]
  lt <- log(d$ntips[ok]); ls <- log(d$nsites[ok]); cl <- d$clade[ok]
  r[ok] <- stats::resid(if (length(unique(cl)) > 1) stats::lm(y ~ cl + lt + ls) else stats::lm(y ~ lt + ls))
  r / stats::sd(r, na.rm = TRUE)
}
