#!/usr/bin/env Rscript
## Collect the per-locus ModelFinder outputs (scripts/submit_modelfinder.sh) into
## two tables in <outdir>:
##
##   best_model.tsv      one row per locus: the BIC-best model, its fit statistics
##                       and every parameter IQ-TREE reports for it -- the six GTR
##                       exchangeabilities, base frequencies, FreeRate rates and
##                       weights (R1..R15), GHOST class weights and tree lengths (H1..H15),
##                       and the tips kept/dropped by dedup_loci.R.
##   all_models.tsv.gz   one row per (locus, model tested): lnL, df, tree length,
##                       BIC and dBIC to the best -- the whole ModelFinder table.
##
## Per-branch lengths of the best model (one set per class for +H) are not tabulated;
## they stay in <outdir>/<clade>/<locus>.treefile.gz.
##
## Usage: scripts/parse_modelfinder.R <outdir> [cores]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("usage: parse_modelfinder.R <outdir> [cores]")
outdir <- normalizePath(args[1])
cores  <- if (length(args) >= 2L) as.integer(args[2]) else 1L

jobs <- utils::read.delim(file.path(outdir, "jobs.tsv"), stringsAsFactors = FALSE)

num <- function(x) suppressWarnings(as.numeric(x))
grab <- function(ln, pat) {   # first number after `pat`, NA if the line is absent
  h <- grep(pat, ln, value = TRUE)
  if (!length(h)) return(NA_real_)
  num(sub(paste0(".*", pat, "\\s*"), "", h[1]) |> sub(pattern = "\\s.*", replacement = ""))
}
nums_after <- function(ln, pat) {   # all numbers on the line after `pat`
  h <- grep(pat, ln, value = TRUE)
  if (!length(h)) return(numeric())
  num(strsplit(trimws(sub(paste0(".*", pat), "", h[1])), "\\s+")[[1]])
}
K <- 15L   # = -cmax

parse_iqtree <- function(f) {
  ln <- readLines(gzfile(f), warn = FALSE)
  out <- list(
    best_model = sub(".*:\\s*", "", grep("^Best-fit model according to BIC:", ln, value = TRUE)[1]),
    lnL    = grab(ln, "Log-likelihood of the tree:"),
    df     = grab(ln, "Number of free parameters \\(#branches \\+ #model parameters\\):"),
    AIC    = grab(ln, "Akaike information criterion \\(AIC\\) score:"),
    AICc   = grab(ln, "Corrected Akaike information criterion \\(AICc\\) score:"),
    BIC    = grab(ln, "Bayesian information criterion \\(BIC\\) score:"),
    tree_length     = grab(ln, "Total tree length \\(sum of branch lengths\\):"),
    internal_length = grab(ln, "Sum of internal branch lengths:")
  )
  ## Exchangeabilities (G-T = 1 reference). Equal for JC/F81, kappa-tied for K2P/HKY.
  for (p in c("A-C", "A-G", "A-T", "C-G", "C-T", "G-T"))
    out[[paste0("r_", sub("-", "", p))]] <- grab(ln, paste0("^\\s*", p, ":"))
  ## Base frequencies: estimated counts for +F, 0.25 each for JC/K2P.
  eq <- any(grepl("State frequencies: \\(equal frequencies\\)", ln))
  for (b in c("A", "C", "G", "T"))
    out[[paste0("pi_", b)]] <- if (eq) 0.25 else grab(ln, paste0("pi\\(", b, "\\) ="))
  out$freq_type <- if (eq) "equal" else
    sub(".*State frequencies:\\s*", "", grep("State frequencies:", ln, value = TRUE)[1])
  ## FreeRate: "Site proportion and rates:  (w1,r1) (w2,r2) ..."
  rw <- grep("^Site proportion and rates:", ln, value = TRUE)
  pr <- if (length(rw)) regmatches(rw[1], gregexpr("\\(([^)]+)\\)", rw[1]))[[1]] else character()
  pr <- do.call(rbind, lapply(strsplit(gsub("[()]", "", pr), ","), num))
  ## GHOST: class weights and per-class total tree lengths.
  hw <- nums_after(ln, "Heterotachy weights:")
  ht <- nums_after(ln, "Heterotachy tree lengths:")
  out$ncat <- if (length(hw)) length(hw) else if (!is.null(pr)) nrow(pr) else 1L
  for (k in seq_len(K)) {
    out[[paste0("R", k, "_rate")]]   <- if (!is.null(pr) && k <= nrow(pr)) pr[k, 2] else NA_real_
    out[[paste0("R", k, "_weight")]] <- if (!is.null(pr) && k <= nrow(pr)) pr[k, 1] else NA_real_
  }
  for (k in seq_len(K)) {
    out[[paste0("H", k, "_weight")]]  <- if (k <= length(hw)) hw[k] else NA_real_
    out[[paste0("H", k, "_treelen")]] <- if (k <= length(ht)) ht[k] else NA_real_
  }
  out
}

parse_models <- function(f, nsites) {   # .model.gz: "MODEL: lnL df treelen"
  ## The file is a checkpoint: skip its non-model keys (initTree, finishedFastMLTree, ...).
  ln <- grep("^(JC|F81|K2P|HKY|TN|GTR)[A-Za-z0-9+*]*: ", readLines(gzfile(f), warn = FALSE), value = TRUE)
  p  <- strsplit(sub(":", "", ln), "\\s+")
  d  <- data.frame(model = vapply(p, `[`, "", 1),
                   lnL = num(vapply(p, `[`, "", 2)), df = num(vapply(p, `[`, "", 3)),
                   tree_length = num(vapply(p, `[`, "", 4)), stringsAsFactors = FALSE)
  d$BIC  <- -2 * d$lnL + d$df * log(nsites)
  d$dBIC <- d$BIC - min(d$BIC)
  d
}

one <- function(i) {
  j <- jobs[i, ]
  base <- file.path(outdir, j$clade, j$locus_id)
  if (!file.exists(paste0(base, ".iqtree.gz"))) return(NULL)
  b <- parse_iqtree(paste0(base, ".iqtree.gz"))
  a <- parse_models(paste0(base, ".model.gz"), j$nsites)
  ## Cross-check: the tabulated BIC-best must be the model IQ-TREE reported.
  b$check_best <- identical(a$model[which.min(a$BIC)], b$best_model)
  b$n_models <- nrow(a)
  list(best = c(list(locus_id = j$locus_id, clade = j$clade, ntips = j$ntips, nsites = j$nsites), b),
       all  = c(list(locus_id = rep(j$locus_id, nrow(a))), as.list(a)))
}

## Bind by column, not with do.call(rbind, <112k data frames>), which takes hours.
bind_cols <- function(rows) {
  nm <- names(rows[[1]])
  as.data.frame(setNames(lapply(nm, function(n) unlist(lapply(rows, `[[`, n), use.names = FALSE)), nm),
                stringsAsFactors = FALSE, check.names = FALSE)
}

res <- parallel::mclapply(seq_len(nrow(jobs)), function(i)
  tryCatch(one(i), error = function(e) structure(conditionMessage(e), class = "parse_error")),
  mc.cores = cores)
bad <- vapply(res, inherits, logical(1), "parse_error")
if (any(bad)) {
  cat(sum(bad), "loci could not be parsed, e.g.", jobs$locus_id[which(bad)[1]], ":",
      unclass(res[[which(bad)[1]]]), "\n")
  res <- res[!bad]
}
res <- res[!vapply(res, is.null, logical(1))]
best <- bind_cols(lapply(res, `[[`, "best"))
all  <- bind_cols(lapply(res, `[[`, "all"))

## Sequences collapsed before fitting (dedup_loci.R); the last row wins if a locus
## was prepared more than once (resubmitted task).
dd <- do.call(rbind, lapply(Sys.glob(file.path(outdir, "dedup", "shard_*.tsv")), utils::read.delim,
                            stringsAsFactors = FALSE, colClasses = c(dropped = "character")))
if (!is.null(dd)) {
  dd <- dd[!duplicated(dd$locus_id, fromLast = TRUE), ]
  k <- match(best$locus_id, dd$locus_id)
  best <- data.frame(best[, 1:3], ntips_kept = dd$ntips_kept[k], best[, -(1:3)],
                     dedup_threshold = dd$threshold[k], dropped = dd$dropped[k],
                     check.names = FALSE, stringsAsFactors = FALSE)
}

utils::write.table(best, file.path(outdir, "best_model.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
gz <- gzfile(file.path(outdir, "all_models.tsv.gz"), "w")
utils::write.table(all, gz, sep = "\t", quote = FALSE, row.names = FALSE); close(gz)
cat(nrow(best), "of", nrow(jobs), "loci parsed;", sum(!best$check_best),
    "where the recomputed BIC-best differs from IQ-TREE's\n")
## Non-zero exit if anything is missing, so a dependent step (archiving) stops.
if (nrow(best) < nrow(jobs) || any(bad)) quit(status = 1)
