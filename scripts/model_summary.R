#!/usr/bin/env Rscript
## Per-locus model parameters for plotting along the genome, derived from each
## locus's OWN best model (best_model.tsv from parse_modelfinder.R).
##
## DESIGN: deliberately not made comparable across models. The point is to show
## model selection as it came out: a locus whose best model is JC has ti/tv = 1 and
## equal base frequencies BECAUSE JC was selected, and that is what the figures show.
## The best model's name, family and number of classes are kept beside every value.
##
##   titv        transition/transversion RATE ratio of the best model:
##               mean(A-G, C-T) / mean(A-C, A-T, C-G, G-T) exchangeability.
##               JC/F81 = 1, K2P/HKY = kappa, TN/GTR = the ratio of means.
##   gc, at_skew, gc_skew   from the model's base frequencies (0.25 each for JC/K2P).
##   rate_cv     coefficient of variation of the FreeRate site-rate distribution,
##               sqrt(sum w_i (r_i - rbar)^2) / rbar; 0 when the best model has one
##               rate class (no +R).
##   r_slow_rate, r_slow_weight / r_fast_rate, r_fast_weight
##               slowest and fastest FreeRate class; NA without +R.
##   tree_length total branch length of the best model on the fixed gene tree
##               (subs/site, summed over ntips_kept tips); mean_branch_length divides
##               by the 2n - 3 branches.
##
## Usage: scripts/model_summary.R <best_model.tsv> <out.tsv>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("usage: model_summary.R <best_model.tsv> <out.tsv>")
b <- utils::read.delim(args[1], stringsAsFactors = FALSE, check.names = FALSE)

m <- b$best_model
subst <- sub("\\+.*", "", m)
rate_family <- ifelse(grepl("\\+R[0-9]+$", m), "R", ifelse(grepl("\\+H[0-9]+$", m), "H", "E"))
k <- rep(1L, length(m))
k[rate_family != "E"] <- as.integer(sub(".*\\+[RH]", "", m[rate_family != "E"]))

## JC/F81 have no rate block to read on some outputs: their exchangeabilities are 1.
r <- b[, paste0("r_", c("AC", "AG", "AT", "CG", "CT", "GT"))]
eqr <- subst %in% c("JC", "F81")
r[eqr, ] <- lapply(r[eqr, ], function(x) ifelse(is.na(x), 1, x))
titv <- rowMeans(r[, c("r_AG", "r_CT")]) / rowMeans(r[, c("r_AC", "r_AT", "r_CG", "r_GT")])

p <- b[, paste0("pi_", c("A", "C", "G", "T"))]

## FreeRate classes: columns R<i>_rate / R<i>_weight, NA beyond k.
K <- sum(grepl("^R[0-9]+_rate$", names(b)))
rr <- as.matrix(b[, paste0("R", seq_len(K), "_rate")])
ww <- as.matrix(b[, paste0("R", seq_len(K), "_weight")])
rbar <- rowSums(rr * ww, na.rm = TRUE)
cv <- sqrt(rowSums(ww * (rr - rbar)^2, na.rm = TRUE)) / rbar
cv[rate_family != "R"] <- 0
first_last <- function(which_fun) {
  i <- apply(rr, 1, function(x) if (all(is.na(x))) NA_integer_ else which_fun(x))
  list(rate = rr[cbind(seq_len(nrow(rr)), i)], weight = ww[cbind(seq_len(nrow(ww)), i)])
}
slow <- first_last(which.min); fast <- first_last(which.max)

nt <- if ("ntips_kept" %in% names(b)) b$ntips_kept else b$ntips
out <- data.frame(
  locus_id = b$locus_id, clade = b$clade, ntips = b$ntips, ntips_kept = nt, nsites = b$nsites,
  best_model = m, subst = subst, freq = ifelse(b$freq_type == "equal", "equal", "F"),
  rate_family = rate_family, k = k,
  titv = titv, gc = p$pi_G + p$pi_C,
  at_skew = (p$pi_A - p$pi_T) / (p$pi_A + p$pi_T), gc_skew = (p$pi_G - p$pi_C) / (p$pi_G + p$pi_C),
  rate_cv = cv, r_slow_rate = slow$rate, r_slow_weight = slow$weight,
  r_fast_rate = fast$rate, r_fast_weight = fast$weight,
  tree_length = b$tree_length, internal_length = b$internal_length,
  mean_branch_length = b$tree_length / (2 * nt - 3),
  lnL = b$lnL, BIC = b$BIC, df = b$df,
  stringsAsFactors = FALSE
)
out <- cbind(out, r, p)
utils::write.table(out, args[2], sep = "\t", quote = FALSE, row.names = FALSE)
cat(nrow(out), "loci ->", args[2], "\n")
print(table(subst = out$subst, rate = paste0(out$rate_family, ifelse(out$rate_family == "E", "", out$k))))
