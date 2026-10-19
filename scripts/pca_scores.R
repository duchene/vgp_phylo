#!/usr/bin/env Rscript
## Per-locus positions on the first two principal components, computed two ways,
## plus the two-panel figure comparing them.
##
## Writes results/published/pc_scores.tsv.gz, results/published/pc_loadings.tsv,
## and figures/fig6_pc_space.{pdf,png}.
##
## Reads results/published/locus_metrics.tsv.gz if present, otherwise the
## uncompressed merged table, so it runs on a bare clone with no data/ directory.
##
## Usage: Rscript scripts/pca_scores.R

root <- Sys.getenv("VGP_PHYLO_ROOT", unset = normalizePath(file.path(dirname(sub("^--file=", "",
          grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "..")))
FIG <- file.path(root, "figures")
PUB <- file.path(root, "results", "published")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(PUB, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------ palette ----
## Identical to scripts/make_figures.R, so figure 6 sits with figures 1 to 5.
SURF <- "#fcfcfb"; INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#898781"
GRID <- "#e1e0d9"; AXIS <- "#c3c2b7"
CAT  <- c("#2a78d6", "#eb6834", "#1baf7a")
SEQ  <- c("#cde2fb","#b7d3f6","#9ec5f4","#86b6ef","#6da7ec","#5598e7",
          "#3987e5","#2a78d6","#256abf","#1c5cab","#184f95","#104281","#0d366b")
DIV_HI <- "#d03b3b"
seq_ramp <- colorRampPalette(c(SURF, SEQ))

## cairo on the Linux workstation. capabilities("cairo") reports TRUE on a Mac
## without XQuartz even though the library then fails to load, so test it by
## opening a throwaway device rather than trusting the flag.
PNG_TYPE <- local({
  f <- tempfile(fileext = ".png")
  ok <- tryCatch({
    suppressWarnings(png(f, type = "cairo")); dev.off()
    file.exists(f) && file.size(f) > 0
  }, error = function(e) FALSE)
  unlink(f)
  if (ok) "cairo" else if (isTRUE(capabilities("aqua"))) "quartz" else "Xlib"
})

render <- function(name, w, h, draw) {
  pdf(file.path(FIG, paste0(name, ".pdf")), width = w, height = h,
      bg = SURF, useDingbats = FALSE)
  draw(); invisible(dev.off())
  png(file.path(FIG, paste0(name, ".png")), width = w * 150, height = h * 150,
      res = 150, bg = SURF, type = PNG_TYPE)
  draw(); invisible(dev.off())
  message("  -> ", name, ".pdf / .png")
}
base_par <- function(...) {
  par(bg = SURF, fg = AXIS, col.axis = MUTED, col.lab = INK2, col.main = INK,
      family = "sans", las = 1, cex.axis = 0.75, cex.lab = 0.85, cex.main = 0.95,
      tcl = -0.2, mgp = c(2.1, 0.45, 0), ...)
}

# --------------------------------------------------------------------- data ----
src <- file.path(PUB, "locus_metrics.tsv.gz")
if (!file.exists(src)) src <- file.path(root, "results", "locus_metrics.tsv")
message("reading ", basename(src), " ...")
d <- read.delim(if (grepl("[.]gz$", src)) gzfile(src) else src, stringsAsFactors = FALSE)

## The same 15 metrics, display names and transforms as scripts/make_figures.R.
M <- rbind(
  data.frame(col="scf_own.scf_prop_sig",             nm="sCF signif. branches"),
  data.frame(col="scf_own.scf_mean",                 nm="sCF mean"),
  data.frame(col="scf_own.sn_mean",                  nm="informative sites/branch"),
  data.frame(col="qirp_own.qirp_mean",               nm="QIRP (P correct)"),
  data.frame(col="qirp_own.qihp_mean",               nm="QIHP (P incorrect)"),
  data.frame(col="delta_plot.delta_mean",            nm="Holland delta"),
  data.frame(col="treedist_species.cid",             nm="CID vs species tree"),
  data.frame(col="aitchison_basefreq.ait_totvar",    nm="compositional het."),
  data.frame(col="rtt_cov.rtt_cov",                  nm="root-to-tip CoV"),
  data.frame(col="site_entropy.entropy_mean_mm",     nm="site entropy"),
  data.frame(col="subst_model.kappa",                nm="kappa (Ti/Tv)"),
  data.frame(col="subst_model.gamma_alpha",          nm="gamma alpha"),
  data.frame(col="stemminess.stemminess_per_branch", nm="stemminess"),
  data.frame(col="tree_length.mean_branch_length",   nm="mean branch length"),
  data.frame(col="site_entropy.unc_lnl_per_site",    nm="uncon. lnL/site"),
  stringsAsFactors = FALSE)

X <- d[, M$col]; colnames(X) <- M$nm
X[["compositional het."]] <- log10(X[["compositional het."]] + 1e-4)
X[["gamma alpha"]]        <- log10(pmin(X[["gamma alpha"]], 100))

grp_of <- function(cl) {
  ifelse(cl == "BirdsAnc0", "Birds",
  ifelse(cl == "SquamataAnc0", "Squamates",
  ifelse(cl == "CartilaginousFishesAnc0", "Cartilaginous fishes",
  ifelse(grepl("^RayFinnedFishes", cl), "Ray-finned fishes",
  ifelse(cl %in% c("AmphibiansAnc0","AnuraAnc0","AnuraAnc1","GymnophionaAnc0"), "Amphibians",
  ifelse(grepl("^Vertebrates", cl), "Deep vertebrates", "Jawless fishes"))))))
}
ok  <- complete.cases(X)
Xc  <- X[ok, ]
id  <- d$locus_id[ok]; clade <- d$clade[ok]; lin <- grp_of(clade)
message(sprintf("  %d loci, %d metrics, %d complete", nrow(d), ncol(X), sum(ok)))

# ---------------------------------------------------------------- two PCAs ----
## UNCORRECTED: metrics as measured, scaled to unit variance. Position on these
## axes reflects which clade the locus came from as much as the locus itself.
pc_raw <- prcomp(Xc, scale. = TRUE)

## CORRECTED: each metric first standardised within its own clade, as a robust z,
## exactly the transform used for the composite in scripts/make_figures.R. Clade
## differences are removed before the axes are found, so position reflects how a
## locus compares with others called against the same ancestral reference.
robust_z <- function(u) {
  sc <- mad(u, na.rm = TRUE)
  if (!is.finite(sc) || sc == 0) sc <- sd(u, na.rm = TRUE)
  if (!is.finite(sc) || sc == 0) return(rep(0, length(u)))
  (u - median(u, na.rm = TRUE)) / sc
}
Z <- sapply(colnames(Xc), function(k) ave(Xc[[k]], clade, FUN = robust_z))
Z[!is.finite(Z)] <- 0
pc_cor <- prcomp(Z, scale. = TRUE)

## Principal component signs are arbitrary, so fix them by rule rather than by
## whatever the decomposition returned. PC1 is oriented so that informative sites
## per branch loads positive, PC2 so that compositional heterogeneity loads
## negative. Both axes then point the same way in both panels and across reruns.
orient <- function(pc) {
  s1 <- sign(pc$rotation["informative sites/branch", 1])
  s2 <- -sign(pc$rotation["compositional het.", 2])
  pc$rotation[, 1] <- pc$rotation[, 1] * s1; pc$x[, 1] <- pc$x[, 1] * s1
  pc$rotation[, 2] <- pc$rotation[, 2] * s2; pc$x[, 2] <- pc$x[, 2] * s2
  pc
}
pc_raw <- orient(pc_raw); pc_cor <- orient(pc_cor)

ve_raw <- summary(pc_raw)$importance[2, 1:2]
ve_cor <- summary(pc_cor)$importance[2, 1:2]
r2 <- function(v) summary(lm(v ~ clade))$r.squared
r2_raw <- c(r2(pc_raw$x[, 1]), r2(pc_raw$x[, 2]))
r2_cor <- c(r2(pc_cor$x[, 1]), r2(pc_cor$x[, 2]))

# -------------------------------------------------------------------- write ----
out <- data.frame(
  locus_id = id, clade = clade, lineage = lin,
  pc1_uncorrected = round(pc_raw$x[, 1], 5), pc2_uncorrected = round(pc_raw$x[, 2], 5),
  pc1_withinclade = round(pc_cor$x[, 1], 5), pc2_withinclade = round(pc_cor$x[, 2], 5),
  stringsAsFactors = FALSE)
con <- gzfile(file.path(PUB, "pc_scores.tsv.gz"), "w")
write.table(out, con, sep = "\t", quote = FALSE, row.names = FALSE); close(con)
message("  -> results/published/pc_scores.tsv.gz  (", nrow(out), " loci)")

load_tab <- data.frame(
  metric = rownames(pc_raw$rotation),
  pc1_uncorrected = round(pc_raw$rotation[, 1], 4), pc2_uncorrected = round(pc_raw$rotation[, 2], 4),
  pc1_withinclade = round(pc_cor$rotation[, 1], 4), pc2_withinclade = round(pc_cor$rotation[, 2], 4),
  row.names = NULL, stringsAsFactors = FALSE)
load_tab <- rbind(load_tab, data.frame(
  metric = c("VARIANCE EXPLAINED", "R2 ON CLADE"),
  pc1_uncorrected = round(c(ve_raw[1], r2_raw[1]), 4), pc2_uncorrected = round(c(ve_raw[2], r2_raw[2]), 4),
  pc1_withinclade = round(c(ve_cor[1], r2_cor[1]), 4), pc2_withinclade = round(c(ve_cor[2], r2_cor[2]), 4),
  stringsAsFactors = FALSE))
write.table(load_tab, file.path(PUB, "pc_loadings.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
message("  -> results/published/pc_loadings.tsv")

# ======================================================================= FIG 6 =
## 112k points, so density rather than a scatter, on the same sequential ramp as
## figure 4. Lineage cannot be a hue at 7 levels in a scatter, so each lineage is
## marked by its centroid instead. Their spread IS the clade confound: wide in the
## left panel, collapsed onto the origin in the right one.
LIN <- c("Birds","Squamates","Cartilaginous fishes","Ray-finned fishes",
         "Amphibians","Deep vertebrates","Jawless fishes")
LAB <- c("Birds","Squamates","Cart. fishes","Ray-finned","Amphibians","Deep vert.","Jawless")

## Push labels apart vertically until they clear each other, then recentre on
## where they started, so the column of labels tracks the points it describes.
spread_y <- function(y, dy) {
  o <- order(y); ys <- y[o]
  for (i in seq_along(ys)[-1]) if (ys[i] - ys[i - 1] < dy) ys[i] <- ys[i - 1] + dy
  ys <- ys - (mean(ys) - mean(y))
  ys[order(o)]
}

pc_panel <- function(px, py, ve, r2v, main, sub) {
  ## Trim the extreme 0.1% per axis: a handful of outliers otherwise set the
  ## limits and squeeze every real locus into the middle sixth of the panel.
  xl <- quantile(px, c(0.0005, 0.9995)); yl <- quantile(py, c(0.0005, 0.9995))
  smoothScatter(px, py, colramp = seq_ramp, nrpoints = 0, xlim = xl, ylim = yl,
                xlab = sprintf("PC1  (%.0f%% of variance)", ve[1] * 100),
                ylab = sprintf("PC2  (%.0f%% of variance)", ve[2] * 100),
                main = main)
  abline(h = 0, v = 0, col = GRID, lwd = 0.6)
  box(col = AXIS, lwd = 0.6)
  mtext(sub, side = 3, line = 0.1, cex = 0.6, col = INK2)
  cx <- tapply(px, factor(lin, LIN), median); cy <- tapply(py, factor(lin, LIN), median)
  ## Direct labels in a column at the right edge, each joined to its centroid by a
  ## leader. Labels at the point itself overprint whenever medians sit close, which
  ## is exactly what the right-hand panel is meant to show.
  ly <- spread_y(cy, 0.075 * diff(yl))
  lx <- xl[1] + 0.995 * diff(xl)
  segments(cx, cy, lx, ly, col = AXIS, lwd = 0.5)
  points(cx, cy, pch = 21, bg = CAT[2], col = SURF, cex = 1.15, lwd = 0.8)
  text(lx, ly, LAB, pos = 2, offset = 0.25, cex = 0.6, col = INK, xpd = NA)
  legend("bottomleft", pch = 21, pt.bg = CAT[2], col = SURF, pt.cex = 1.0,
         legend = "lineage median", bty = "n", cex = 0.62, text.col = INK2)
}

render("fig6_pc_space", 10, 5.2, function() {
  base_par(mfrow = c(1, 2), mar = c(3.9, 3.9, 3.0, 1.2), oma = c(2.4, 0, 2.0, 0))
  pc_panel(pc_raw$x[, 1], pc_raw$x[, 2], ve_raw, r2_raw,
           "Uncorrected: metrics as measured",
           sprintf("clade explains %.0f%% of PC1 and %.0f%% of PC2",
                   r2_raw[1] * 100, r2_raw[2] * 100))
  pc_panel(pc_cor$x[, 1], pc_cor$x[, 2], ve_cor, r2_cor,
           "Corrected: metrics standardised within clade",
           sprintf("clade explains %.0f%% of PC1 and %.0f%% of PC2",
                   r2_cor[1] * 100, r2_cor[2] * 100))
  mtext("Where each locus sits in metric space, before and after removing clade",
        outer = TRUE, line = 0.4, cex = 0.95, col = INK, font = 2)
  mtext(paste0("Shading = locus density (n = ", format(nrow(out), big.mark = ","),
               "). PC1 rises with informative sites per branch, PC2 falls with ",
               "compositional heterogeneity. Axis limits trim the extreme 0.1%."),
        outer = TRUE, side = 1, line = 0.8, cex = 0.66, col = INK2)
})

message("done.")
