#!/usr/bin/env Rscript
## Produce the committed result artifacts from the local working outputs.
##
## The working files (results/locus_metrics.tsv, results/locus_quality.tsv) are
## gitignored: together they are 120 MB, and R writes 15 significant digits, which
## for concordance factors and probabilities is noise past about the seventh.
##
## This writes results/published/ instead, which IS committed:
##   locus_metrics.tsv.gz   all 112,199 loci x 64 columns, 6 significant digits
##   locus_quality.tsv.gz   composite scores and within-clade ranks
##   data_dictionary.tsv    every column: source metric, meaning, direction
##   summary_stats.tsv      n / NAs / quantiles per column, so the tables can be
##                          sanity-checked without downloading them
##
## Rounding to 6 significant digits halves the gzipped size (39 -> 21 MB) for a
## maximum relative change of 5e-06. Separating "working" from "published" also
## keeps regenerated full-precision files from showing up as noisy diffs.
##
## Usage: Rscript scripts/export_results.R

root <- Sys.getenv("VGP_PHYLO_ROOT", unset = normalizePath(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "..")))
source(file.path(root, "config", "project.R"))

OUT <- file.path(root, "results", "published")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
SIGDIG <- 6L

write_gz <- function(d, name) {
  num <- vapply(d, is.numeric, logical(1))
  d[num] <- lapply(d[num], signif, digits = SIGDIG)
  p <- file.path(OUT, paste0(name, ".tsv.gz"))
  con <- gzfile(p, "wt", compression = 6)
  utils::write.table(d, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  close(con)
  message(sprintf("  %-24s %6.1f MB  (%d x %d)", basename(p),
                  file.size(p) / 1048576, nrow(d), ncol(d)))
  invisible(d)
}

# ------------------------------------------------------------------- tables ----
message("exporting tables ...")
met <- read.delim(PATHS$merged, stringsAsFactors = FALSE)
qual_path <- file.path(root, "results", "locus_quality.tsv")
qual <- if (file.exists(qual_path)) read.delim(qual_path, stringsAsFactors = FALSE) else NULL

write_gz(met, "locus_metrics")
if (!is.null(qual)) write_gz(qual, "locus_quality") else
  message("  (no locus_quality.tsv -- run scripts/make_figures.R first)")

# --------------------------------------------------------- data dictionary ----
## Curated, because an auto-generated dictionary of column names tells a reader
## nothing they could not already guess. `direction` records which way is
## DESIRABLE for phylogenomics; NA means there is no defensible direction, and
## those columns are described but never scored.
D <- function(col, meaning, direction = NA, note = "")
  data.frame(column = col, meaning = meaning, better = direction, note = note,
             stringsAsFactors = FALSE)

dict <- rbind(
  D("locus_id", "Locus identifier: <Anc>_<Anc>refChr<N>_<start>_<end>"),
  D("clade",    "Ancestral reference genome the locus was called against (22 levels)"),
  D("ntips",    "Sequences in the alignment; equals nspecies (loci are single-copy)"),
  D("nsites",   "Alignment columns", NA, "2500 is a pipeline cap; 51.7% of loci sit at it"),

  D("aitchison_basefreq.ait_mean",  "Mean Aitchison distance between per-sequence base compositions", "lower", "compositional heterogeneity among taxa"),
  D("aitchison_basefreq.ait_var",   "Variance of those pairwise distances", NA),
  D("aitchison_basefreq.ait_max",   "Maximum pairwise Aitchison distance", NA),
  D("aitchison_basefreq.ait_totvar","Compositional total variance = mean squared distance to the compositional centre", "lower", "standard one-number dispersion summary; preferred over ait_mean"),
  D("aitchison_basefreq.ait_to_centroid_mean", "Mean Aitchison distance to the compositional centre", "lower"),
  D("aitchison_basefreq.gc_mean",   "Mean GC content across sequences", NA),
  D("aitchison_basefreq.gc_sd",     "SD of GC content across sequences", "lower"),
  D("aitchison_basefreq.nsites_used","Alignment columns used", NA),
  D("aitchison_basefreq.ntips_used","Sequences with >= 20 unambiguous bases", NA),

  D("delta_plot.delta_mean", "Holland's delta over 2000 sampled quartets: departure from treelikeness", "lower", "0 = perfectly treelike, 1 = fully boxy/reticulate"),
  D("delta_plot.delta_var",  "Variance of delta across quartets", NA),
  D("delta_plot.delta_q90",  "90th percentile of delta", "lower", "catches a reticulate minority a mean would hide"),
  D("delta_plot.q_used",     "Quartets that resolved", NA),
  D("delta_plot.ntips_used", "Sequences used", NA),

  D("qirp_own.qirp_mean", "Mean quartet internode resolution probability (P correct) over internodes", "higher", "Townsend, Su & Tekle 2012, via PhyInformR::Approximator"),
  D("qirp_own.qihp_mean", "Mean probability an internode is resolved INCORRECTLY (homoplasy)", "lower", "no other metric captures actively misleading signal"),
  D("qirp_own.qirp_var_rel", "Variance of QIRP across internodes / m(1-m)", NA, "normalised because raw variance of a [0,1] quantity tracks its mean mechanically"),
  D("qirp_own.n_internodes", "Internodes evaluated", NA, "capped at 100 per locus so effort is size-independent"),
  D("qirp_own.n_internodes_total", "Internodes available before capping", NA),
  D("qirp_own.rate_var", "Variance of per-site rates, scaled to mean 1", NA, "IS the locus's rate heterogeneity; see subst_model.gamma_alpha for the ML counterpart"),
  D("qirp_own.prop_invariant_rate", "Share of zero-rate (invariant) sites", NA),

  D("rtt_cov.rtt_cov",   "CoV of root-to-tip distances, midpoint-rooted", "lower", "among-lineage rate variation; scale-free"),
  D("rtt_cov.rtt_mean",  "Mean root-to-tip distance", NA),
  D("rtt_cov.ntips_used","Tips used", NA),

  D("scf_own.scf_mean",  "Mean site concordance factor against the locus's own gene tree", NA, "no defensible direction: high sCF also tracks shallow, poorly resolved trees"),
  D("scf_own.scf_var",   "Variance of sCF across branches", NA),
  D("scf_own.sdf1_mean", "Mean site discordance factor, alternative quartet 1", NA),
  D("scf_own.sdf2_mean", "Mean site discordance factor, alternative quartet 2", NA),
  D("scf_own.sn_mean",   "Mean informative sites per branch", "higher", "raw information content"),
  D("scf_own.n_branches","Branches with informative sites (sN > 0)", NA),
  D("scf_own.scf_prop_sig", "Proportion of branches whose sCF beats p=1/3 (binomial, alpha 0.001)", "higher", "a proportion, so it does not scale with tree size"),

  D("site_entropy.entropy_mean",    "Mean per-site Shannon entropy (plug-in, natural log)", NA, "biased downward by ~(k-1)/2n; prefer the _mm columns"),
  D("site_entropy.entropy_var",     "Variance of per-site entropy (plug-in)", NA),
  D("site_entropy.entropy_mean_mm", "Mean per-site entropy, Miller-Madow bias-corrected", NA, "USE THIS; the correction removes an artefactual correlation with taxon count"),
  D("site_entropy.entropy_var_mm",  "Variance of bias-corrected per-site entropy", NA),
  D("site_entropy.entropy_mean_norm","Mean entropy / log 4, i.e. fraction of maximal site disorder", NA),
  D("site_entropy.prop_invariant",  "Share of sites with a single observed state", NA),
  D("site_entropy.unc_lnl_per_site","Unconstrained (multinomial) log-likelihood per site", NA, "likelihood ceiling no model can exceed; the one retained PhyloMAd statistic"),
  D("site_entropy.n_patterns",      "Distinct site patterns", NA, "NOT comparable across loci: rises with both taxon count and length"),
  D("site_entropy.nsites_used",     "Sites with >= 2 unambiguous bases", NA),
  D("site_entropy.ntips_used",      "Sequences used", NA),

  D("stemminess.stemminess",           "Share of total tree length on internal branches", NA, "rises with taxon count for combinatorial reasons"),
  D("stemminess.stemminess_per_branch","Stemminess / proportion of branches that are internal", NA, "USE THIS: the size-corrected form"),
  D("stemminess.n_internal",           "Internal branches", NA),
  D("stemminess.n_branches",           "All branches (2n-3 when fully resolved)", NA),

  D("subst_model.kappa",       "ML transition/transversion rate ratio under HKY", NA, "prefer this over titv: less saturation-biased (rho -0.55 vs -0.76 with divergence)"),
  D("subst_model.gamma_alpha", "ML gamma shape parameter (+G4)", NA, "low alpha = strong across-site rate heterogeneity; IQ-TREE ceiling 999 for 2 loci"),

  D("titv.titv",       "mean(TS distance) / mean(TV distance) over sequence pairs", NA, "model-free, NOT kappa; values < 1 are expected"),
  D("titv.ts_mean",    "Mean transition distance", NA),
  D("titv.tv_mean",    "Mean transversion distance", NA),
  D("titv.pdist_mean", "Mean raw p-distance", NA, "divergence context for the other columns"),

  D("tree_length.mean_branch_length",  "Total tree length / number of branches", NA, "USE THIS as the overall rate proxy"),
  D("tree_length.median_branch_length","Median branch length", NA),
  D("tree_length.mean_root_to_tip",    "Mean tip depth, midpoint-rooted", NA, "per-lineage accumulated substitutions"),
  D("tree_length.total_length",        "Sum of branch lengths", NA, "NOT comparable across loci: grows with taxon count"),
  D("tree_length.n_branches",          "Branches", NA),
  D("tree_length.prop_zero_length",    "Share of ~zero-length branches", "lower", "how much of the topology is effectively unresolved"),

  D("treedist_species.cid", "Clustering Information Distance to the ASTRAL species tree, normalised", "lower", "0 = identical. Circular if the goal is to infer that tree"),
  D("treedist_species.discord_support_weighted", "Share of the gene tree's total branch support on splits absent from the species tree", "lower"),
  D("treedist_species.ntips_shared", "Taxa shared with the species tree", NA)
)

miss <- setdiff(names(met), dict$column)
extra <- setdiff(dict$column, names(met))
if (length(miss))  warning("columns absent from the dictionary: ", paste(miss, collapse = ", "))
if (length(extra)) warning("dictionary lists non-existent columns: ", paste(extra, collapse = ", "))
dict <- dict[match(intersect(names(met), dict$column), dict$column), ]

utils::write.table(dict, file.path(OUT, "data_dictionary.tsv"), sep = "\t",
                   quote = FALSE, row.names = FALSE, na = "")
message(sprintf("  %-24s %d columns documented (%d undocumented)",
                "data_dictionary.tsv", nrow(dict), length(miss)))

# ----------------------------------------------------------- summary stats ----
num <- vapply(met, is.numeric, logical(1))
ss <- do.call(rbind, lapply(names(met)[num], function(k) {
  x <- met[[k]]
  q <- stats::quantile(x, c(0, .05, .25, .5, .75, .95, 1), na.rm = TRUE)
  data.frame(column = k, n = sum(!is.na(x)), n_na = sum(is.na(x)),
             mean = signif(mean(x, na.rm = TRUE), SIGDIG),
             sd = signif(stats::sd(x, na.rm = TRUE), SIGDIG),
             min = signif(q[1], SIGDIG), p05 = signif(q[2], SIGDIG),
             p25 = signif(q[3], SIGDIG), median = signif(q[4], SIGDIG),
             p75 = signif(q[5], SIGDIG), p95 = signif(q[6], SIGDIG),
             max = signif(q[7], SIGDIG), row.names = NULL, stringsAsFactors = FALSE)
}))
utils::write.table(ss, file.path(OUT, "summary_stats.tsv"), sep = "\t",
                   quote = FALSE, row.names = FALSE)
message(sprintf("  %-24s %d numeric columns", "summary_stats.tsv", nrow(ss)))
message("done -- results/published/ is the committed artifact set")
