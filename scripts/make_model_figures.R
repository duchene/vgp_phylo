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
## 2-D histograms (log counts): 112k points would overplot. One blue sequential ramp.
dens <- function(x, y, xlab, ylab, main, logx = FALSE, logy = FALSE, nb = 70, note = NULL) {
  ok <- is.finite(x) & is.finite(y) & (!logx | x > 0) & (!logy | y > 0)
  x <- x[ok]; y <- y[ok]; tx <- if (logx) log10(x) else x; ty <- if (logy) log10(y) else y
  qx <- quantile(tx, c(0.002, 0.998)); qy <- quantile(ty, c(0.002, 0.998))
  bx <- seq(qx[1], qx[2], length.out = nb + 1); by <- seq(qy[1], qy[2], length.out = nb + 1)
  h <- table(cut(tx, bx, include.lowest = TRUE), cut(ty, by, include.lowest = TRUE))
  z <- log10(unclass(h)); z[!is.finite(z)] <- NA
  base_par(mar = c(3.6, 3.8, 2.6, 1))
  image(head(bx, -1) + diff(bx) / 2, head(by, -1) + diff(by) / 2, z, col = seq_ramp(40)[6:40],
        xlab = "", ylab = "", axes = FALSE)
  lab_axis <- function(side, lg, rng) {
    if (lg) { at <- pretty(rng); axis(side, at = at, labels = signif(10^at, 2), col = AXIS) }
    else axis(side, col = AXIS)
  }
  lab_axis(1, logx, qx); lab_axis(2, logy, qy)
  title(xlab = xlab, line = 2.3); title(ylab = ylab, line = 2.6)
  title(main = main, adj = 0, line = 0.8)
  r <- cor(x, y, method = "spearman")
  mtext(sprintf("Spearman %.2f   n = %s%s", r, format(length(x), big.mark = ","),
                if (is.null(note)) "" else paste0("   ", note)),
        side = 3, line = -1.1, adj = 0.02, cex = 0.58, col = INK2)
}
render("A3_parameter_pairs", 10, 9, function() {
  par(mfrow = c(2, 2))
  dens(p$gc, p$titv, "GC (model base frequencies)", "ti/tv rate ratio", "a  GC vs ti/tv",
       note = sprintf("stripe at 0.50: %s K2P loci (equal freqs)", format(sum(p$freq == "equal"), big.mark = ",")))
  dens(p$mean_branch_length, p$rate_cv, "tree length / branches (subs/site)", "FreeRate CV of site rates",
       "b  Overall rate vs rate heterogeneity", logx = TRUE,
       note = sprintf("CV = 0 for %s single-class loci", format(sum(p$k == 1), big.mark = ",")))
  dens(p$gc, p$mean_branch_length, "GC (model base frequencies)", "tree length / branches (subs/site)",
       "c  GC vs overall rate", logy = TRUE,
       note = sprintf("stripe at 0.50: %s K2P loci", format(sum(p$freq == "equal"), big.mark = ",")))
  dens(p$at_skew, p$gc_skew, "AT skew (A-T)/(A+T)", "GC skew (G-C)/(G+C)", "d  Strand composition skews",
       note = sprintf("point at 0,0: %s K2P loci", format(sum(p$freq == "equal"), big.mark = ",")))
})
message("done")
