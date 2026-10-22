#!/usr/bin/env Rscript
## Summary figures for assessing which loci are best for phylogenomics.
##
## Reads results/locus_metrics.tsv and writes figures/ as PDF (vector, for
## publication) plus PNG (for quick viewing). Also writes
## results/locus_quality.tsv -- the per-locus composite scores behind figure 5.
##
## Usage: Rscript scripts/make_figures.R

root <- Sys.getenv("VGP_PHYLO_ROOT", unset = normalizePath(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "..")))
source(file.path(root, "config", "project.R"))
FIG <- PATHS$figures
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------ palette ----
## Reference data-viz palette, used unchanged (its validation results are
## documented: first three categorical slots clear all-pairs CVD/normal-vision
## floors; blue is the sequential hue; blue<->red with a GRAY midpoint is the
## diverging pair). Never a rainbow, never a hue at the diverging midpoint.
SURF   <- "#fcfcfb"; INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#898781"
GRID   <- "#e1e0d9"; AXIS <- "#c3c2b7"
CAT    <- c("#2a78d6", "#eb6834", "#1baf7a")            # slots 1-3 only
SEQ    <- c("#cde2fb","#b7d3f6","#9ec5f4","#86b6ef","#6da7ec","#5598e7",
            "#3987e5","#2a78d6","#256abf","#1c5cab","#184f95","#104281","#0d366b")
DIV_LO <- "#2a78d6"; DIV_MID <- "#f0efec"; DIV_HI <- "#d03b3b"
GOOD   <- "#0ca30c"

## Empty regions must recede to the chart surface, not sit on the lightest
## blue step -- otherwise "no data" reads as "a little data".
seq_ramp <- colorRampPalette(c(SURF, SEQ))
div_ramp <- colorRampPalette(c(DIV_LO, DIV_MID, DIV_HI))

## Draw one figure twice: PDF (vector, for publication) and PNG (for viewing).
## The draw function is called once per device so both get identical output.
render <- function(name, w, h, draw) {
  pdf(file.path(FIG, paste0(name, ".pdf")), width = w, height = h,
      bg = SURF, useDingbats = FALSE)
  draw(); invisible(dev.off())
  png(file.path(FIG, paste0(name, ".png")), width = w * 150, height = h * 150,
      res = 150, bg = SURF, type = "cairo")
  draw(); invisible(dev.off())
  message("  -> ", name, ".pdf / .png")
}
base_par <- function(...) {
  par(bg = SURF, fg = AXIS, col.axis = MUTED, col.lab = INK2, col.main = INK,
      family = "sans", las = 1, cex.axis = 0.75, cex.lab = 0.85, cex.main = 0.95,
      tcl = -0.2, mgp = c(2.1, 0.45, 0), ...)
}
grid_h <- function(at) abline(h = at, col = GRID, lwd = 0.6)
grid_v <- function(at) abline(v = at, col = GRID, lwd = 0.6)

# --------------------------------------------------------------------- data ----
message("reading merged metrics ...")
d <- read.delim(PATHS$merged, stringsAsFactors = FALSE)

## The metric set carried through the figures, with display names and the
## direction that is DESIRABLE for phylogenomics. dir = 0 means "no defensible
## direction" -- such metrics are described but never scored (see FIG 5).
M <- rbind(
  data.frame(col="scf_own.scf_prop_sig",                    nm="sCF signif. branches", dir= 1, grp="signal"),
  data.frame(col="scf_own.scf_mean",                        nm="sCF mean",             dir= 0, grp="signal"),
  data.frame(col="scf_own.sn_mean",                         nm="informative sites/branch", dir= 1, grp="signal"),
  data.frame(col="qirp_own.qirp_mean",                      nm="QIRP (P correct)",     dir= 1, grp="signal"),
  data.frame(col="qirp_own.qihp_mean",                      nm="QIHP (P incorrect)",   dir=-1, grp="conflict"),
  data.frame(col="delta_plot.delta_mean",                   nm="Holland delta",        dir=-1, grp="conflict"),
  data.frame(col="treedist_species.cid",                    nm="CID vs species tree",  dir=-1, grp="conflict"),
  data.frame(col="aitchison_basefreq.ait_totvar",           nm="compositional het.",   dir=-1, grp="bias"),
  data.frame(col="rtt_cov.rtt_cov",                         nm="root-to-tip CoV",      dir=-1, grp="bias"),
  data.frame(col="site_entropy.entropy_mean_mm",            nm="site entropy",         dir= 0, grp="bias"),
  data.frame(col="subst_model.kappa",                       nm="kappa (Ti/Tv)",        dir= 0, grp="model"),
  data.frame(col="subst_model.gamma_alpha",                 nm="gamma alpha",          dir= 0, grp="model"),
  data.frame(col="stemminess.stemminess_per_branch",        nm="stemminess",           dir= 0, grp="tree"),
  data.frame(col="tree_length.mean_branch_length",          nm="mean branch length",   dir= 0, grp="tree"),
  data.frame(col="site_entropy.unc_lnl_per_site",           nm="uncon. lnL/site",      dir= 0, grp="tree"),
  stringsAsFactors = FALSE)

X <- d[, M$col]; colnames(X) <- M$nm
## Two columns need a transform before any distance, correlation or PCA:
## compositional heterogeneity spans 4 orders of magnitude, and gamma alpha has a
## long tail plus IQ-TREE's 999 ceiling (2 loci).
X[["compositional het."]] <- log10(X[["compositional het."]] + 1e-4)
X[["gamma alpha"]]        <- log10(pmin(X[["gamma alpha"]], 100))

## Lineage groups: 22 clades is far more than any categorical palette supports, so
## clade is shown by FACET, never by hue. Groups are the natural vertebrate
## lineages behind the ancestral references.
grp_of <- function(cl) {
  ifelse(cl == "BirdsAnc0", "Birds",
  ifelse(cl == "SquamataAnc0", "Squamates",
  ifelse(cl == "CartilaginousFishesAnc0", "Cartilaginous fishes",
  ifelse(grepl("^RayFinnedFishes", cl), "Ray-finned fishes",
  ifelse(cl %in% c("AmphibiansAnc0","AnuraAnc0","AnuraAnc1","GymnophionaAnc0"), "Amphibians",
  ifelse(grepl("^Vertebrates", cl), "Deep vertebrates", "Jawless fishes"))))))
}
d$lineage <- grp_of(d$clade)
LIN <- c("Birds","Squamates","Cartilaginous fishes","Ray-finned fishes",
         "Amphibians","Deep vertebrates","Jawless fishes")

ok <- complete.cases(X)
Xc <- X[ok, ]; clade <- d$clade[ok]; lin <- d$lineage[ok]
message(sprintf("  %d loci, %d metrics, %d complete", nrow(d), ncol(X), sum(ok)))

# ============================================================== FIG 1 ==========
## How much of each metric is just "which clade is this?" -- the finding that
## governs every downstream ranking. Form: magnitude across categories -> bars.
## One series, so one hue and no legend.
message("fig 1: clade confound")
r2 <- vapply(colnames(Xc), function(k)
  summary(lm(Xc[[k]] ~ clade))$r.squared, numeric(1))
r2 <- sort(r2)
pc <- prcomp(Xc, scale. = TRUE)
r2_pc1 <- summary(lm(pc$x[,1] ~ clade))$r.squared

render("fig1_clade_confound", 7.5, 5, function() {
  base_par(mar = c(3.6, 11.5, 2.6, 1.2))
  bp <- barplot(r2, horiz = TRUE, names.arg = names(r2), col = CAT[1], border = NA,
                xlim = c(0, 0.8), xlab = expression("variance explained by clade  (" * R^2 * ")"),
                main = "Most of each metric is explained by which clade the locus came from",
                cex.names = 0.75, xaxs = "i")
  grid_v(seq(0, 0.8, 0.2)); box(col = AXIS, lwd = 0.6)
  barplot(r2, horiz = TRUE, col = CAT[1], border = NA, add = TRUE, axes = FALSE, names.arg = rep("", length(r2)))
  ## Direct labels on the bars rather than a value on every axis tick.
  text(r2 + 0.012, bp, sprintf("%.2f", r2), adj = 0, cex = 0.62, col = INK2, xpd = NA)
  abline(v = 0.5, col = DIV_HI, lwd = 1, lty = 2)
  text(0.51, 1.2, "half the variance", col = DIV_HI, cex = 0.62, adj = 0, srt = 90)
  mtext(sprintf("PC1 of all metrics: R² = %.2f. Rank loci WITHIN clade, not globally.", r2_pc1),
        side = 1, line = 2.5, cex = 0.68, col = INK2, adj = 0)
})

# ============================================================== FIG 2 ==========
## Which metrics are redundant. Correlation has polarity, so a diverging ramp with
## a gray midpoint; ordered by hierarchical clustering so blocks are visible.
message("fig 2: redundancy")
C <- cor(Xc, method = "spearman")
ord <- hclust(as.dist(1 - abs(C)), method = "average")$order
Co <- C[ord, ord]; n <- ncol(Co)

render("fig2_metric_redundancy", 7.6, 6.8, function() {
  base_par(mar = c(9.5, 9.5, 2.8, 3.4))
  ## Show one triangle only: the matrix is symmetric, so drawing both halves
  ## doubles the ink and makes the sparse value labels look arbitrary.
  Cm <- Co; Cm[lower.tri(Cm, diag = TRUE)] <- NA
  image(1:n, 1:n, Cm, col = div_ramp(101), zlim = c(-1, 1), axes = FALSE,
        xlab = "", ylab = "",
        main = "Spearman correlation between metrics (clustered)")
  axis(1, at = 1:n, labels = colnames(Co), las = 2, cex.axis = 0.68, col = NA, col.ticks = NA)
  axis(2, at = 1:n, labels = colnames(Co), las = 1, cex.axis = 0.68, col = NA, col.ticks = NA)
  ## 2px surface gap between cells, per the mark spec.
  abline(h = seq(0.5, n + 0.5), v = seq(0.5, n + 0.5), col = SURF, lwd = 1.6)
  ## Label only the strong pairs -- never a number in every cell.
  for (a in 1:n) for (b in 1:n) if (a < b && abs(Co[a,b]) >= 0.55)
    text(a, b, sprintf("%.2f", Co[a,b]), cex = 0.52,
         col = ifelse(abs(Co[a,b]) > 0.8, SURF, INK))
  box(col = AXIS, lwd = 0.6)
  ## Colour key
  xk <- grconvertX(0.965, "ndc", "user"); yk <- seq(1, n)
  ys <- seq(min(yk), max(yk), length.out = 101)
  rect(xk, ys[-101], xk + (n * 0.022), ys[-1], col = div_ramp(101), border = NA, xpd = NA)
  text(xk + n * 0.028, c(min(ys), mean(ys), max(ys)), c("-1", "0", "+1"),
       adj = 0, cex = 0.6, col = INK2, xpd = NA)
  mtext(paste("Near-duplicate pairs (|rho| >= 0.55) are labelled: informative sites/sCF-significance,",
              "branch length/sCF, branch length/kappa, QIRP/QIHP, composition/entropy."),
        side = 1, line = 8.2, cex = 0.62, col = INK2, adj = 0)
})

# ============================================================== FIG 3 ==========
## Dimensionality and what the axes mean. Two panels, one series each.
message("fig 3: PCA")
ve <- summary(pc)$importance[2, 1:8]
render("fig3_pca_axes", 9.5, 4.6, function() {
  base_par(mfrow = c(1, 2), mar = c(3.8, 3.8, 2.8, 1.2))
  bp <- barplot(ve * 100, names.arg = paste0("PC", 1:8), col = CAT[1], border = NA,
                ylim = c(0, 45), ylab = "variance explained (%)",
                main = "Four axes carry three quarters of the variation")
  grid_h(seq(0, 40, 10)); box(col = AXIS, lwd = 0.6)
  barplot(ve * 100, col = CAT[1], border = NA, add = TRUE, axes = FALSE, names.arg = rep("", 8))
  text(bp, ve * 100 + 1.6, sprintf("%.0f", ve * 100), cex = 0.62, col = INK2)
  text(bp[4], 38, sprintf("cumulative to PC4: %.0f%%", 100 * sum(ve[1:4])),
       cex = 0.68, col = INK2, adj = 0.5)

  ## A biplot was tried here and abandoned: six metrics load near the origin with
  ## near-identical angles, so their labels overprint into an unreadable blob
  ## whatever radial offset is applied. A loadings heatmap has no collisions, and
  ## shows four components instead of two. Loadings have sign, so: diverging ramp,
  ## gray midpoint.
  np <- 4
  L <- pc$rotation[, 1:np, drop = FALSE]
  L <- L[order(L[, 1]), , drop = FALSE]
  nm <- nrow(L)
  par(mar = c(3.2, 11.5, 2.8, 3.0))
  image(1:np, 1:nm, t(L), col = div_ramp(101), zlim = c(-0.6, 0.6), axes = FALSE,
        xlab = "", ylab = "", main = "Metric loadings on the first four axes")
  axis(1, at = 1:np, labels = sprintf("PC%d\n%.0f%%", 1:np, ve[1:np] * 100),
       cex.axis = 0.68, col = NA, col.ticks = NA, padj = 0.5)
  axis(2, at = 1:nm, labels = rownames(L), las = 1, cex.axis = 0.66,
       col = NA, col.ticks = NA)
  abline(h = seq(0.5, nm + 0.5), v = seq(0.5, np + 0.5), col = SURF, lwd = 1.6)
  for (a in 1:np) for (b in 1:nm) if (abs(L[b, a]) >= 0.25)
    text(a, b, sprintf("%+.2f", L[b, a]), cex = 0.5,
         col = ifelse(abs(L[b, a]) > 0.45, SURF, INK))
  box(col = AXIS, lwd = 0.6)
  xk <- np + 0.75; ys <- seq(1, nm, length.out = 101)
  rect(xk, ys[-101], xk + 0.22, ys[-1], col = div_ramp(101), border = NA, xpd = NA)
  text(xk + 0.30, c(1, (1 + nm) / 2, nm), c("-0.6", "0", "+0.6"),
       adj = 0, cex = 0.58, col = INK2, xpd = NA)
})

# ============================================================== FIG 4 ==========
## The trade-off that makes "best" non-obvious: per-site concordance and
## informativeness pull in opposite directions along divergence. Shown WITHIN
## clade (figure 1's finding): each metric is a robust z on its own clade, so the
## density is the trade-off among comparable loci, not the gap between clades.
## Row 2 asks how much of QIRP is anything beyond branch length.
## 112k points, so density (sequential ramp), not a scatter.
message("fig 4: trade-off, within clade")
wz <- function(u) {                 # robust z within clade: (x - median) / MAD
  ave(u, clade, FUN = function(v) {
    sc <- mad(v); if (!is.finite(sc) || sc == 0) sc <- sd(v)
    if (!is.finite(sc) || sc == 0) return(v * 0)
    (v - median(v)) / sc
  })
}
## rank residual on clade + log tips + log sites, for the size-free rho
nt <- log(d$ntips[ok]); ns <- log(d$nsites[ok])
sz <- function(u) resid(lm(rank(u) ~ clade + nt + ns))
W  <- as.data.frame(lapply(Xc, wz), check.names = FALSE)
## CID is discrete on small trees and some clades have a tiny MAD, so its z bands
## into stripes; use the within-clade percentile instead.
PCT <- "CID vs species tree"
W[[PCT]] <- 100 * ave(Xc[[PCT]], clade, FUN = function(v) (rank(v) - 0.5) / length(v))
S  <- as.data.frame(lapply(Xc, sz), check.names = FALSE)
lim <- function(u) quantile(u, c(0.005, 0.995))
## raw-scale density panel, also used by figure 5
den_panel <- function(x, y, xlab, ylab, main, sub = NULL) {
  fin <- is.finite(x) & is.finite(y)
  smoothScatter(x[fin], y[fin], colramp = seq_ramp, nrpoints = 0,
                xlab = xlab, ylab = ylab, main = main)
  box(col = AXIS, lwd = 0.6)
  if (!is.null(sub)) mtext(sub, side = 3, line = 0.1, cex = 0.6, col = INK2)
  rho <- suppressWarnings(cor(x[fin], y[fin], method = "spearman"))
  legend("topright", legend = sprintf("rho = %+.2f", rho), bty = "n",
         cex = 0.68, text.col = INK2)
}
## the same on within-clade z, with the size-free rho beside it
zden_panel <- function(kx, ky, xlab, ylab, main) {
  x <- W[[kx]]; y <- W[[ky]]
  smoothScatter(x, y, colramp = seq_ramp, nrpoints = 0, xlim = lim(x), ylim = lim(y),
                xlab = xlab, ylab = ylab, main = main)
  abline(h = 0, v = 0, col = GRID, lwd = 0.6)
  box(col = AXIS, lwd = 0.6)
  r_w <- cor(x, y, method = "spearman")
  r_s <- cor(S[[kx]], S[[ky]])
  legend("topright", legend = c(sprintf("within clade  rho = %+.2f", r_w),
                                sprintf("+ size-free   rho = %+.2f", r_s)),
         bty = "o", bg = adjustcolor(SURF, 0.85), box.col = NA, cex = 0.62, text.col = INK2)
}
zl <- function(k) paste0(k, if (k == PCT) "  (percentile in clade)" else "  (z in clade)")
tradeoff_pairs <- list(
  c("mean branch length",       "sCF mean",            "Divergence vs per-site concordance"),
  c("informative sites/branch", "sCF mean",            "Information vs concordance"),
  c("QIRP (P correct)",         "CID vs species tree", "Predicted power vs reference agreement"),
  c("Holland delta",            "CID vs species tree", "Treelikeness vs reference agreement"),
  c("QIRP (P correct)",         "sCF mean",            "Predicted power vs per-site concordance"),
  c("QIRP (P correct)",         "compositional het.", "Predicted power vs composition"),
  c("QIRP (P correct)",         "Holland delta",       "Predicted power vs reticulation"),
  c("QIRP (P correct)",         "mean branch length",  "Predicted power vs divergence"))
## per-lineage within-clade rho, to the log, so the pooled panels can be checked
for (pr in tradeoff_pairs) {
  r <- vapply(split(seq_along(lin), lin), function(i)
    cor(W[[pr[1]]][i], W[[pr[2]]][i], method = "spearman"), numeric(1))
  message(sprintf("  %-26s ~ %-20s %s", pr[1], pr[2],
                  paste(sprintf("%s %+.2f", substr(names(r), 1, 5), r), collapse = "  ")))
}
render("fig4_tradeoff", 13, 6.8, function() {
  base_par(mfrow = c(2, 4), mar = c(3.8, 3.9, 3.0, 1.0), oma = c(2.2, 0, 2.0, 0))
  for (pr in tradeoff_pairs) zden_panel(pr[1], pr[2], zl(pr[1]), zl(pr[2]), pr[3])
  mtext("Metrics disagree about what a good locus looks like, within clade", outer = TRUE,
        line = 0.4, cex = 0.95, col = INK, font = 2)
  mtext(paste("Axes: robust z within clade ((x - median) / MAD), CID as percentile within clade;",
              "shading = locus density (n = 112,199),",
              "central 99% shown. Size-free rho: rank residual on clade + log tips + log sites.",
              "Bottom row: QIRP against four other axes."),
        outer = TRUE, side = 1, line = 0.6, cex = 0.62, col = INK2)
})

# ============================================================== FIG 5 ==========
## A transparent composite, standardised WITHIN clade (figure 1's finding), plus
## the Pareto alternative that needs no weights.
message("fig 5: composite quality")
scored <- M[M$dir != 0, ]
## Robust within-clade z: (x - median) / MAD, so a heavy tail in one clade does not
## dominate, then signed so higher is always better.
##
## NOTE: use ave(), not `tapply(...)[clade]`. The latter looks like a per-locus
## lookup but returns the WHOLE clade vector for every locus, so unlisting it tries
## to build a ~10^9-element vector -- it hangs rather than erroring.
robust_z <- function(u) {
  sc <- mad(u, na.rm = TRUE)
  if (!is.finite(sc) || sc == 0) sc <- sd(u, na.rm = TRUE)
  if (!is.finite(sc) || sc == 0) return(rep(0, length(u)))
  (u - median(u, na.rm = TRUE)) / sc
}
Z <- sapply(scored$nm, function(k) ave(Xc[[k]], clade, FUN = robust_z))
Z <- Z * matrix(scored$dir, nrow(Z), ncol(Z), byrow = TRUE)
Z[!is.finite(Z)] <- 0
score <- rowMeans(Z)
q <- quantile(score, c(0.9, 0.99), na.rm = TRUE)

## Pareto front on the four least-redundant scored axes.
## Non-dominated within clade, on a 4000-locus subsample per clade for tractability.
## Weight-free alternative to the composite: a locus must beat its own clade's
## median on EVERY scored axis at once. This replaced a Pareto front, which is
## O(n^2) per clade and took >11 min on 112k loci -- far too slow for what it added.
## This rule is O(n) and answers the same question: good on all counts, no weights.
allgood <- rowSums(Z > 0) == ncol(Z)
message(sprintf("  %d loci beat their clade median on all %d axes", sum(allgood), ncol(Z)))

render("fig5_locus_quality", 10, 6.4, function() {
  base_par(mfrow = c(2, 3), mar = c(3.8, 3.9, 3.0, 1.2), oma = c(2.4, 0, 2.0, 0))
  h <- hist(score, breaks = 90, plot = FALSE)
  plot(h, col = CAT[1], border = NA, xlab = "composite quality (within-clade z)",
       ylab = "loci", main = "Composite score", freq = TRUE)
  grid_h(pretty(h$counts)); box(col = AXIS, lwd = 0.6)
  plot(h, col = CAT[1], border = NA, add = TRUE, freq = TRUE)
  abline(v = q, col = GOOD, lwd = 1.3, lty = c(2, 1))
  text(q, max(h$counts) * c(0.95, 0.8), c("top 10%", "top 1%"), col = GOOD,
       cex = 0.62, adj = -0.1)

  bxp <- split(score, factor(lin, levels = LIN))
  par(mar = c(7.6, 3.9, 3.0, 1.2))
  boxplot(bxp, col = CAT[1], border = INK2, outline = FALSE, las = 2,
          ylab = "composite quality", main = "By lineage (centred within clade)",
          cex.axis = 0.62, lwd = 0.7, staplewex = 0)
  grid_h(pretty(score)); abline(h = 0, col = MUTED, lwd = 0.7, lty = 3)
  box(col = AXIS, lwd = 0.6)

  par(mar = c(3.8, 3.9, 3.0, 1.2))
  den_panel(Xc[["QIRP (P correct)"]], Xc[["sCF signif. branches"]],
            "QIRP (P correct)", "sCF signif. branches", "Two signal axes")
  pts <- which(score >= q[2])
  points(Xc[["QIRP (P correct)"]][pts], Xc[["sCF signif. branches"]][pts],
         pch = 16, cex = 0.16, col = GOOD)
  legend("bottomright", legend = "top 1%", pch = 16, col = GOOD, pt.cex = 0.7,
         bty = "n", cex = 0.66, text.col = INK2)

  ## How many top loci each clade contributes, as a share of its own size --
  ## within-clade standardisation should make this roughly flat by construction;
  ## departures show clades whose distributions are too tight to standardise well.
  ## Both tables must share levels: the smallest clades (2 loci) contribute nothing
  ## to the global top 10%, so a bare table() of the subset drops them and the
  ## division is non-conformable.
  cf <- factor(clade, levels = sort(unique(clade)))
  tab <- as.numeric(table(cf[score >= q[1]])) / as.numeric(table(cf))
  names(tab) <- levels(cf)
  tab <- sort(tab[is.finite(tab)])
  par(mar = c(3.8, 9.8, 3.0, 1.2))
  barplot(as.numeric(tab), horiz = TRUE, names.arg = names(tab), col = CAT[1],
          border = NA, cex.names = 0.5, xlab = "share of clade in global top 10%",
          main = "Standardisation check")
  abline(v = 0.10, col = DIV_HI, lty = 2, lwd = 1); box(col = AXIS, lwd = 0.6)

  par(mar = c(3.8, 3.9, 3.0, 1.2))
  den_panel(Xc[["compositional het."]], Xc[["QIHP (P incorrect)"]],
            expression(log[10] * " compositional heterogeneity"), "QIHP (P incorrect)",
            "Two risk axes")
  points(Xc[["compositional het."]][pts], Xc[["QIHP (P incorrect)"]][pts],
         pch = 16, cex = 0.16, col = GOOD)

  h2 <- hist(score[allgood], breaks = 60, plot = FALSE)
  plot(h2, col = CAT[3], border = NA, xlab = "composite quality",
       ylab = "loci", main = sprintf("Good on ALL %d axes (n = %s)", ncol(Z),
                                     format(sum(allgood), big.mark = ",")))
  grid_h(pretty(h2$counts)); box(col = AXIS, lwd = 0.6)
  plot(h2, col = CAT[3], border = NA, add = TRUE)

  mtext("Ranking loci for phylogenomics", outer = TRUE, line = 0.4, cex = 0.95,
        col = INK, font = 2)
  mtext(paste("Composite = mean of 8 robust within-clade z-scores, signed so higher is better.",
              "Weights are equal by default -- see scripts/make_figures.R to change them."),
        outer = TRUE, side = 1, line = 0.9, cex = 0.68, col = INK2)
})

# ------------------------------------------------------------------- output ----
out <- data.frame(locus_id = d$locus_id[ok], clade = clade, lineage = lin,
                  quality = score, good_on_all_axes = allgood,
                  rank_within_clade = ave(-score, clade, FUN = function(u) rank(u, ties.method = "min")),
                  Z, check.names = FALSE)
out <- out[order(-out$quality), ]
utils::write.table(out, file.path(root, "results", "locus_quality.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("wrote results/locus_quality.tsv (%s loci); figures in %s",
                format(nrow(out), big.mark = ","), FIG))
