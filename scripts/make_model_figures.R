#!/usr/bin/env Rscript
## Figure set A of the best-model landscape (docs/PLAN, Step 1): what model selection
## chose, and how the chosen models' parameters relate to one another and to the
## existing locus metrics. No genomic positions needed.
##
## Values are each locus's OWN best model (results/modelfinder/model_params.tsv),
## deliberately not made comparable: e.g. ti/tv is exactly 1 under JC/F81, and
## rate_cv is 0 when the best model has a single rate class.
##
## Writes figures/model/A1..A3 as PDF + PNG.
## Usage: Rscript scripts/make_model_figures.R

root <- Sys.getenv("VGP_PHYLO_ROOT", unset = normalizePath(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "..")))
source(file.path(root, "config", "project.R"))
FIG <- file.path(PATHS$figures, "model")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------ palette ----
## Same reference palette and conventions as scripts/make_figures.R. Categorical
## slots 1-4 are used only where the four substitution models sit ADJACENT (stacked
## bars), the pairing the palette is validated for; never all-pairs (scatter).
SURF   <- "#fcfcfb"; INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#898781"
GRID   <- "#e1e0d9"; AXIS <- "#c3c2b7"
CAT    <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100")
SEQ    <- c("#cde2fb","#b7d3f6","#9ec5f4","#86b6ef","#6da7ec","#5598e7",
            "#3987e5","#2a78d6","#256abf","#1c5cab","#184f95","#104281","#0d366b")
DIV_LO <- "#2a78d6"; DIV_MID <- "#f0efec"; DIV_HI <- "#d03b3b"
seq_ramp <- colorRampPalette(c(SURF, SEQ))
div_ramp <- colorRampPalette(c(DIV_LO, DIV_MID, DIV_HI))

render <- function(name, w, h, draw) {
  pdf(file.path(FIG, paste0(name, ".pdf")), width = w, height = h, bg = SURF, useDingbats = FALSE)
  draw(); invisible(dev.off())
  png(file.path(FIG, paste0(name, ".png")), width = w * 150, height = h * 150, res = 150,
      bg = SURF, type = "cairo")
  draw(); invisible(dev.off())
  message("  -> model/", name, ".pdf / .png")
}
base_par <- function(...) {
  par(bg = SURF, fg = AXIS, col.axis = MUTED, col.lab = INK2, col.main = INK,
      family = "sans", las = 1, cex.axis = 0.75, cex.lab = 0.85, cex.main = 0.95,
      tcl = -0.2, mgp = c(2.1, 0.45, 0), ...)
}
grid_h <- function(at) abline(h = at, col = GRID, lwd = 0.6)
grid_v <- function(at) abline(v = at, col = GRID, lwd = 0.6)
legend_box <- function(...) legend(..., bty = "n", text.col = INK2, cex = 0.72)

# --------------------------------------------------------------------- data ----
p <- read.delim(file.path(PATHS$root, "results", "modelfinder", "model_params.tsv"),
                stringsAsFactors = FALSE)
SUBST <- c("GTR", "TN", "HKY", "K2P", "F81", "JC")   # most to least parameter-rich
SUBST <- SUBST[SUBST %in% p$subst]
stopifnot(length(SUBST) <= length(CAT))
clades <- names(sort(table(p$clade), decreasing = TRUE))
clades <- clades[table(p$clade)[clades] >= 50]          # smaller clades: too few to show shares

# ------------------------------------------------- A1: what was selected ----
render("A1_model_selection", 12, 6.2, function() {
  layout(matrix(1:3, 1), widths = c(1.25, 1.25, 1))
  ## (a) substitution model share per clade
  base_par(mar = c(5.2, 11.5, 2.4, 1))
  tab <- prop.table(table(factor(p$clade, levels = rev(clades)), factor(p$subst, levels = SUBST)), 1)
  b <- barplot(t(tab), horiz = TRUE, col = CAT[seq_along(SUBST)], border = SURF, lwd = 1.5,
               names.arg = sprintf("%s (%s)", rownames(tab), format(table(p$clade)[rownames(tab)], big.mark = ",")),
               axes = FALSE, space = 0.25, xlim = c(0, 1))
  axis(1, at = seq(0, 1, 0.25), labels = paste0(seq(0, 100, 25), "%"), col = AXIS)
  title(main = "a  Best substitution model", adj = 0, line = 1)
  legend_box("bottom", inset = c(0, -0.13), horiz = TRUE, xpd = NA, legend = SUBST,
             fill = CAT[seq_along(SUBST)], border = SURF)
  ## (b) number of rate classes per clade (ordinal: sequential ramp from step 250)
  kk <- sort(unique(p$k))
  kcol <- colorRampPalette(SEQ[4:13])(length(kk))
  lab <- paste0("k=", kk)
  tabk <- prop.table(table(factor(p$clade, levels = rev(clades)), factor(p$k, levels = kk)), 1)
  base_par(mar = c(5.2, 1, 2.4, 1))
  barplot(t(tabk), horiz = TRUE, col = kcol, border = SURF, lwd = 1.5, names.arg = rep("", nrow(tabk)),
          axes = FALSE, space = 0.25, xlim = c(0, 1))
  axis(1, at = seq(0, 1, 0.25), labels = paste0(seq(0, 100, 25), "%"), col = AXIS)
  title(main = "b  Rate classes (k) of the best model", adj = 0, line = 1)
  legend_box("bottom", inset = c(0, -0.13), horiz = TRUE, xpd = NA, legend = lab, fill = kcol, border = SURF)
  nh <- sum(p$rate_family == "H")
  mtext(sprintf("k=1: no +R.  +H (GHOST) best on %d loci, counted under its k", nh), side = 3, line = 0.1,
        cex = 0.55, col = MUTED, adj = 0)
  ## (c) mean k by data size: more tips and sites -> more classes detectable
  tb <- cut(p$ntips_kept, c(0, 25, 50, 100, 150, 200, 300, 400, 600), dig.lab = 4)
  sb <- cut(p$nsites, c(1199, 1500, 1800, 2100, 2400, 2500), dig.lab = 4)
  mk <- tapply(p$k, list(tb, sb), mean); nk <- table(tb, sb); mk[nk < 20] <- NA
  base_par(mar = c(5.2, 6.2, 2.4, 3.6))
  brks <- seq(1, ceiling(max(mk, na.rm = TRUE)), length.out = 9)
  image(seq_len(nrow(mk)), seq_len(ncol(mk)), mk, col = seq_ramp(10)[3:10], breaks = brks,
        axes = FALSE, xlab = "", ylab = "")
  axis(1, at = seq_len(nrow(mk)), labels = sub("\\((\\d+),(\\d+)\\]", "\\1-\\2", rownames(mk)), las = 2, col = AXIS)
  axis(2, at = seq_len(ncol(mk)), labels = sub("\\((\\d+),(\\d+)\\]", "\\1-\\2", colnames(mk)), col = AXIS)
  title(xlab = "tips (after collapsing)", line = 2.9); title(ylab = "sites", line = 4.6)
  title(main = "c  Mean k by data size", adj = 0, line = 1)
  ## colour key
  usr <- par("usr"); x0 <- usr[2] + 0.25; ys <- seq(usr[3], usr[4], length.out = length(brks))
  rect(x0, ys[-length(ys)], x0 + 0.35, ys[-1], col = seq_ramp(10)[3:10], border = SURF, xpd = NA)
  text(x0 + 0.45, ys[c(1, length(ys))], sprintf("%.1f", brks[c(1, length(brks))]), adj = 0,
       cex = 0.65, col = INK2, xpd = NA)
  mtext("blank: < 20 loci", side = 1, line = 4, cex = 0.55, col = MUTED, adj = 1)
})

# ------------------------------------------------ A2: correlation matrix ----
VARS <- c(titv = "ti/tv (rate ratio)", gc = "GC (model freqs)", at_skew = "AT skew",
          gc_skew = "GC skew", rate_cv = "FreeRate CV", r_slow_weight = "slowest class weight",
          r_slow_rate = "slowest class rate", k = "rate classes k",
          tree_length = "tree length", mean_branch_length = "tree length / branches")
EXT <- c("scf_own.scf_mean" = "sCF mean", "delta_plot.delta_mean" = "delta (Holland)",
         "site_entropy.entropy_mean" = "site entropy", "qirp_own.rate_var" = "parsimony rate var",
         "aitchison_basefreq.ait_mean" = "composition heterog.", "rtt_cov.rtt_cov" = "root-to-tip CV",
         "treedist_species.cid" = "gene-species tree dist")
lm_ <- read.delim(PATHS$merged, stringsAsFactors = FALSE, check.names = FALSE)
X <- cbind(p[, names(VARS)], lm_[match(p$locus_id, lm_$locus_id), names(EXT)])
colnames(X) <- c(VARS, EXT)
## Partial version: Spearman on ranks after removing log tips, log sites and clade.
resid_rank <- function(x) {
  r <- rep(NA_real_, length(x)); ok <- !is.na(x)
  r[ok] <- resid(lm(rank(x[ok]) ~ log(p$ntips_kept[ok]) + log(p$nsites[ok]) + p$clade[ok]))
  r
}
C1 <- cor(X, method = "spearman", use = "pairwise.complete.obs")
C2 <- cor(sapply(X, resid_rank), use = "pairwise.complete.obs")
ord <- hclust(as.dist(1 - abs(C1)), "average")$order
render("A2_parameter_correlations", 13, 6.6, function() {
  layout(matrix(1:3, 1), widths = c(1, 1, 0.12))
  for (i in 1:2) {
    C <- list(C1, C2)[[i]][ord, ord]; n <- ncol(C)
    base_par(mar = c(9.5, if (i == 1) 10 else 1, 3, 0.5))
    image(1:n, 1:n, C[, n:1], col = div_ramp(21), breaks = seq(-1, 1, length.out = 22), axes = FALSE,
          xlab = "", ylab = "")
    abline(h = 0.5 + 0:n, v = 0.5 + 0:n, col = SURF, lwd = 1.5)
    axis(1, at = 1:n, labels = colnames(C), las = 2, col = NA, col.axis = INK2, cex.axis = 0.72)
    if (i == 1) axis(2, at = n:1, labels = colnames(C), col = NA, col.axis = INK2, cex.axis = 0.72)
    big <- which(abs(C) >= 0.4 & row(C) != col(C), arr.ind = TRUE)   # label only strong pairs
    text(big[, 2], n + 1 - big[, 1], sprintf("%.2f", C[big]), cex = 0.5,
         col = ifelse(abs(C[big]) > 0.7, SURF, INK))
    title(main = c("a  Spearman correlation across loci",
                   "b  Partial: after log tips, log sites and clade")[i], adj = 0, line = 1)
  }
  base_par(mar = c(9.5, 0.5, 3, 3))
  image(1, seq(-1, 1, length.out = 21), matrix(seq(-1, 1, length.out = 21), 1), col = div_ramp(21),
        axes = FALSE, xlab = "", ylab = "")
  axis(4, at = c(-1, -0.5, 0, 0.5, 1), col = NA, col.axis = INK2)
})

# ---------------------------------------------- A3: pairwise relationships ----
## One row per clade (clades with >= 1,000 loci), one column per pair, axes SHARED
## down each column so clades compare directly. Rate per branch is strongly clade-
## structured (~0.025 subs/site in Birds and deep Vertebrates vs 0.06-0.09 in fishes
## and Anura), so pooled panels would mix groups. 2-D histograms of log counts,
## scaled within each panel; bins coarser for smaller clades.
PAIRS <- list(
  list(x = "gc", y = "titv", xl = "GC (model base frequencies)", yl = "ti/tv rate ratio",
       main = "GC vs ti/tv"),
  list(x = "mean_branch_length", y = "rate_cv", xl = "tree length / branches", yl = "FreeRate CV",
       main = "Rate vs rate heterogeneity", logx = TRUE),
  list(x = "gc", y = "mean_branch_length", xl = "GC (model base frequencies)", yl = "tree length / branches",
       main = "GC vs rate", logy = TRUE),
  list(x = "at_skew", y = "gc_skew", xl = "AT skew (A-T)/(A+T)", yl = "GC skew (G-C)/(G+C)",
       main = "Strand composition skews"))
tr <- function(v, lg) if (lg) log10(v) else v
for (k in seq_along(PAIRS)) {   # shared limits per column, from all loci
  pr <- PAIRS[[k]]; lx <- isTRUE(pr$logx); ly <- isTRUE(pr$logy)
  PAIRS[[k]]$xlim <- quantile(tr(p[[pr$x]], lx), c(0.002, 0.998), na.rm = TRUE)
  PAIRS[[k]]$ylim <- quantile(tr(p[[pr$y]], ly), c(0.002, 0.998), na.rm = TRUE)
}
dens <- function(dd, pr, first_row, last_row) {
  lx <- isTRUE(pr$logx); ly <- isTRUE(pr$logy)
  x <- tr(dd[[pr$x]], lx); y <- tr(dd[[pr$y]], ly)
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  nb <- if (length(x) >= 10000) 60 else if (length(x) >= 3000) 40 else 28
  bx <- seq(pr$xlim[1], pr$xlim[2], length.out = nb + 1); by <- seq(pr$ylim[1], pr$ylim[2], length.out = nb + 1)
  h <- table(cut(x, bx, include.lowest = TRUE), cut(y, by, include.lowest = TRUE))
  z <- log10(unclass(h)); z[!is.finite(z)] <- NA
  base_par(mar = c(if (last_row) 3.4 else 0.8, 3.4, if (first_row) 2.2 else 0.6, 0.6))
  plot(NA, xlim = pr$xlim, ylim = pr$ylim, axes = FALSE, xlab = "", ylab = "", xaxs = "i", yaxs = "i")
  if (any(!is.na(z)))
    image(head(bx, -1) + diff(bx) / 2, head(by, -1) + diff(by) / 2, z, col = seq_ramp(40)[6:40],
          add = TRUE)
  lab_axis <- function(side, lg, rng, labels = TRUE) {
    at <- pretty(rng); at <- at[at >= rng[1] & at <= rng[2]]
    axis(side, at = at, labels = if (!labels) FALSE else if (lg) signif(10^at, 2) else at, col = AXIS)
  }
  lab_axis(1, lx, pr$xlim, labels = last_row); lab_axis(2, ly, pr$ylim)
  if (last_row) title(xlab = pr$xl, line = 2.2)
  if (first_row) title(main = pr$main, adj = 0, line = 0.7)
  r <- suppressWarnings(cor(x, y, method = "spearman"))
  mtext(sprintf("rho %.2f", r), side = 3, line = -1, adj = 0.03, cex = 0.55, col = INK2)
}
A3_CLADES <- clades[table(p$clade)[clades] >= 1000]
render("A3_parameter_pairs", 12, 2.15 * length(A3_CLADES) + 1.2, function() {
  par(mfrow = c(length(A3_CLADES), length(PAIRS)), oma = c(0.5, 9.5, 2.6, 0.5))
  for (cl in A3_CLADES) {
    dd <- p[p$clade == cl, ]
    for (k in seq_along(PAIRS)) {
      dens(dd, PAIRS[[k]], first_row = cl == A3_CLADES[1], last_row = cl == tail(A3_CLADES, 1))
      if (k == 1) {
        nk <- sum(dd$freq == "equal"); n1 <- sum(dd$k == 1)
        mtext(sprintf("%s\n%s loci\n%s K2P (GC 0.50)\n%s single-class (CV 0)", cl,
                      format(nrow(dd), big.mark = ","), format(nk, big.mark = ","), format(n1, big.mark = ",")),
              side = 2, line = 3, las = 1, adj = 1, cex = 0.55, col = INK2)
      }
    }
  }
  mtext(paste("Each clade's own best-model values; axes shared down each column; shading = log count within panel.",
              "Clades with >= 1,000 loci. Model-fixed values: K2P loci sit at GC 0.50 and skews 0,0;",
              "single-class loci at CV 0."), side = 3, outer = TRUE, line = 1, adj = 0, cex = 0.55, col = MUTED)
})
message("done")
