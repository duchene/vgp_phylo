#!/usr/bin/env Rscript
## Figure set C: how the locus metrics relate to the substitution-model results,
## kept to a few axes rather than every pairwise correlation.
##
##   C1  variance partitioning: for each model output, the share explained by clade,
##       by data size (log tips, log sites) and by the locus metrics, uniquely and
##       jointly; and the same for each metric with the model outputs as the third set.
##   C2  redundancy analysis (RDA, a multi-response linear model reduced to its main
##       axes) of the model outputs on the metrics, within clade and size-free, as a
##       correlation circle; small panels repeat it within each lineage.
##   C3  the composition link: GC against sCF, ti/tv, rate-heterogeneity gain and
##       compositional heterogeneity, within clade, per lineage.
##   C4  hub variables: a few variables that stand for the others (results/links/).
##
## Variables and transforms: scripts/lib/link_data.R (four metrics that duplicate a
## model output are left out there). Within clade and size-free means: rank, then the
## residual on clade + log tips + log sites (as in figure A2b).
##
## Writes figures/model/C1..C4 (PDF + PNG), results/links/hub_selection.tsv and
## results/links/rda_summary.tsv. Usage: Rscript scripts/make_link_figures.R

root <- Sys.getenv("VGP_PHYLO_ROOT", unset = normalizePath(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "..")))
source(file.path(root, "config", "project.R"))
source(file.path(PATHS$root, "scripts", "lib", "link_data.R"))
FIG <- file.path(PATHS$figures, "model")
OUT <- file.path(PATHS$root, "results", "links")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------ palette ----
SURF   <- "#fcfcfb"; INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#898781"
GRID   <- "#e1e0d9"; AXIS <- "#c3c2b7"
CAT    <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100")
SEQ    <- c("#cde2fb","#b7d3f6","#9ec5f4","#86b6ef","#6da7ec","#5598e7",
            "#3987e5","#2a78d6","#256abf","#1c5cab","#184f95","#104281","#0d366b")
DIV_LO <- "#2a78d6"; DIV_MID <- "#f0efec"; DIV_HI <- "#d03b3b"
seq_ramp <- colorRampPalette(c(SURF, SEQ))
div_ramp <- colorRampPalette(c(DIV_LO, DIV_MID, DIV_HI))
COL_METRIC <- CAT[1]; COL_MODEL <- CAT[2]   # the two blocks, everywhere in set C

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
      family = "sans", las = 1, cex.axis = 0.72, cex.lab = 0.8, cex.main = 0.85,
      tcl = -0.2, mgp = c(2, 0.4, 0), ...)
}
grid_h <- function(at) abline(h = at, col = GRID, lwd = 0.6)
grid_v <- function(at) abline(v = at, col = GRID, lwd = 0.6)
legend_box <- function(..., cex = 0.72) legend(..., bty = "n", text.col = INK2, cex = cex)

# --------------------------------------------------------------------- data ----
message("reading ...")
d <- load_link_data()
METRIC <- VARS$key[VARS$block == "metric"]; MODEL <- VARS$key[VARS$block == "model"]
d <- d[complete.cases(d[, VARS$key]), ]
message(sprintf("  %d loci with every variable", nrow(d)))
R <- sapply(VARS$key, function(v) resid_clade_size(d[[v]], d))   # within clade, size-free

# ------------------------------------------------- C1: variance partitioning ----
## Three predictor sets A (clade), B (size), C (the other block). Unique fraction of
## a set = R2(all three) - R2(the other two); shared = R2(all) - sum of uniques.
## Plain least squares on the transformed values (not ranks), so the clade and size
## fractions are on the scale the variables are reported in.
r2 <- function(y, M) { f <- lm.fit(cbind(1, M), y); 1 - sum(f$residuals^2) / sum((y - mean(y))^2) }
partition <- function(ys, third) {
  A <- model.matrix(~ clade, d)[, -1]; B <- cbind(log(d$ntips), log(d$nsites))
  C <- as.matrix(d[, third])
  t(sapply(ys, function(v) {
    y <- d[[v]]; full <- r2(y, cbind(A, B, C))
    u <- c(clade = full - r2(y, cbind(B, C)), size = full - r2(y, cbind(A, C)),
           other = full - r2(y, cbind(A, B)))
    c(u, shared = full - sum(u), full = full)
  }))
}
vp_model  <- partition(MODEL, METRIC)
vp_metric <- partition(METRIC, MODEL)
PART_COL <- c(clade = "#9ec5f4", size = CAT[4], other = NA, shared = "#c3c2b7")
render("C1_variance_partitioning", 13, 6.2, function() {
  layout(matrix(1:2, 1), widths = c(1, 1.3))
  for (i in 1:2) {
    vp <- list(vp_model, vp_metric)[[i]]
    vp <- vp[order(vp[, "other"]), , drop = FALSE]
    oc <- c(COL_METRIC, COL_MODEL)[i]
    cols <- PART_COL; cols["other"] <- oc
    base_par(mar = c(4, 10.5, 3.6, 1))
    m <- t(pmax(vp[, c("other", "shared", "clade", "size")], 0))
    barplot(m, horiz = TRUE, col = cols[rownames(m)], border = SURF, lwd = 1.2, space = 0.3,
            names.arg = VLAB[colnames(m)], axes = FALSE, xlim = c(0, 1), cex.names = 0.75)
    grid_v(seq(0, 1, 0.25))
    barplot(m, horiz = TRUE, col = cols[rownames(m)], border = SURF, lwd = 1.2, space = 0.3,
            names.arg = rep("", ncol(m)), axes = FALSE, xlim = c(0, 1), add = TRUE)
    axis(1, at = seq(0, 1, 0.25), labels = paste0(seq(0, 100, 25), "%"), col = AXIS)
    title(xlab = "share of variance across loci (R2)", line = 2.2)
    title(main = c("a  Model outputs: what explains them", "b  Locus metrics: what explains them")[i],
          adj = 0, line = 2)
    lab <- c(c("locus metrics (unique)", "model outputs (unique)")[i], "shared", "clade (unique)",
             "data size (unique)")
    legend_box("top", inset = c(0, -0.085), horiz = TRUE, xpd = NA, legend = lab,
               fill = cols[c("other", "shared", "clade", "size")], border = SURF, cex = 0.62)
  }
  mtext(paste("Least squares on the transformed variables. Unique = R2 lost when that set is dropped;",
              "shared = explained but not attributable to one set; blank = unexplained.",
              sprintf("%s loci.", format(nrow(d), big.mark = ","))),
        side = 1, outer = TRUE, line = -1.2, adj = 0.02, cex = 0.55, col = MUTED)
})

# ----------------------------------------------- C2: redundancy analysis ----
## Y (model outputs) on X (metrics), both within clade and size-free (columns of R,
## unit SD). The fitted values are what X explains of Y; their principal axes are the
## RDA axes. Arrows are correlations with the RDA locus scores (a correlation circle),
## so metric and model arrows are on one scale. Lineage panels: the same pooled axes,
## correlations recomputed within each lineage, i.e. whether the structure replicates.
Xr <- R[, METRIC]; Yr <- R[, MODEL]
fit <- lm.fit(Xr, Yr)
Fh <- fit$fitted.values
pa <- prcomp(Fh, center = FALSE)
explained <- sum(Fh^2) / sum(Yr^2)
axis_share <- pa$sdev^2 / sum(apply(Yr, 2, var))          # of total Y variance
## orient: FreeRate CV points right on RDA1, GC up on RDA2, for stable reruns
sc <- Fh %*% pa$rotation[, 1:2]
sgn <- c(sign(cor(Yr[, "rate_cv"], sc[, 1])), sign(cor(Yr[, "gc"], sc[, 2])))
sc <- sweep(sc, 2, sgn, `*`)
cc <- cancor(Xr, Yr)$cor
rda_corr <- function(ix) rbind(cor(Xr[ix, ], sc[ix, ]), cor(Yr[ix, ], sc[ix, ]))
C_all <- rda_corr(seq_len(nrow(d)))
write.table(data.frame(quantity = c("Y variance explained by X", paste0("RDA", 1:4, " share of Y variance"),
                                    paste0("canonical correlation ", 1:4)),
                       value = round(c(explained, axis_share[1:4], cc[1:4]), 4)),
            file.path(OUT, "rda_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
## Push apart labels on the same side of the circle that would overlap vertically.
spread <- function(y, side, gap) {
  for (sd_ in c(TRUE, FALSE)) {
    j <- which(side == sd_); if (length(j) < 2) next
    o <- j[order(y[j])]
    for (it in 1:50) {
      moved <- FALSE
      for (a in seq_along(o)[-1]) if (y[o[a]] - y[o[a - 1]] < gap) {
        mv <- (gap - (y[o[a]] - y[o[a - 1]])) / 2
        y[o[a]] <- y[o[a]] + mv; y[o[a - 1]] <- y[o[a - 1]] - mv; moved <- TRUE
      }
      if (!moved) break
    }
  }
  y
}
circle <- function(Cm, main, big, n = NULL) {
  base_par(mar = if (big) c(5.2, 3.4, 3.2, 1) else c(1.2, 1.2, 2, 0.6))
  lim <- if (big) c(-1.08, 1.08) else c(-1, 1)
  plot(NA, xlim = lim, ylim = lim, asp = 1, axes = FALSE, xlab = "", ylab = "")
  th <- seq(0, 2 * pi, length.out = 200)
  for (rr in c(0.5, 1)) lines(rr * cos(th), rr * sin(th), col = GRID, lwd = 0.8)
  abline(h = 0, v = 0, col = GRID, lwd = 0.8)
  isM <- rownames(Cm) %in% MODEL
  arrows(0, 0, Cm[, 1], Cm[, 2], length = if (big) 0.07 else 0.03, lwd = if (big) 1.8 else 1.2,
         col = ifelse(isM, COL_MODEL, COL_METRIC))
  ## labels: all arrows >= 0.2 long on the big panel, the 3 longest on small ones
  len <- sqrt(rowSums(Cm^2))
  show <- if (big) len >= 0.2 else rank(-len) <= 3
  x <- Cm[show, 1]; y <- Cm[show, 2]; r <- pmax(sqrt(x^2 + y^2), 1e-9)
  ly <- spread(y + 0.05 * y / r, x >= 0, if (big) 0.055 else 0.12)
  for (right in c(TRUE, FALSE)) {
    j <- (x >= 0) == right
    if (any(j)) text(x[j] + 0.05 * x[j] / r[j], ly[j], VLAB[rownames(Cm)[show][j]], cex = if (big) 0.62 else 0.5,
                     adj = c(if (right) 0 else 1, 0.5), col = ifelse(isM[show][j], COL_MODEL, COL_METRIC), xpd = NA)
  }
  if (big) {
    axis(1, at = c(-1, -0.5, 0, 0.5, 1), col = AXIS); axis(2, at = c(-1, -0.5, 0, 0.5, 1), col = AXIS)
    title(xlab = sprintf("correlation with RDA1 (%.0f%% of model-output variance)", 100 * axis_share[1]), line = 2.1)
    title(ylab = sprintf("correlation with RDA2 (%.0f%%)", 100 * axis_share[2]), line = 2.1)
  }
  title(main = main, adj = 0, line = if (big) 1.8 else 0.6, cex.main = if (big) 0.85 else 0.72)
  if (!is.null(n)) mtext(sprintf("%s loci", format(n, big.mark = ",")), side = 3, line = -0.4, adj = 0,
                         cex = 0.5, col = MUTED)
}
render("C2_rda_correlation_circle", 13, 7.4, function() {
  par(oma = c(0, 0, 1.6, 0))
  layout(cbind(1, matrix(2:7, 2, byrow = TRUE)), widths = c(2.1, 1, 1, 1))
  short <- names(which(sqrt(rowSums(C_all^2)) < 0.2))
  circle(C_all, "a  Model outputs on locus metrics (RDA), within clade and size-free", TRUE)
  legend_box("bottomleft", inset = c(-0.02, -0.01), legend = c("locus metric", "model output"),
             col = c(COL_METRIC, COL_MODEL), lwd = 2)
  mtext(sprintf("Metrics explain %.0f%% of model-output variance; RDA1+2 hold %.0f%% of that. Canonical r: %s.",
                100 * explained, 100 * sum(axis_share[1:2]) / explained, paste(sprintf("%.2f", cc[1:3]), collapse = ", ")),
        side = 1, line = 3.3, adj = 0, cex = 0.52, col = MUTED)
  if (length(short)) mtext(paste0("Unlabelled (|r| < 0.2): ", paste(VLAB[short], collapse = ", "), "."),
                           side = 1, line = 4.1, adj = 0, cex = 0.52, col = MUTED)
  for (L in LINEAGES) {
    ix <- which(d$lineage == L)
    circle(rda_corr(ix), L, FALSE, length(ix))
  }
  mtext("b  The same axes, correlations recomputed within each lineage (3 longest arrows labelled)",
        side = 3, outer = TRUE, line = 0.2, adj = 0.97, cex = 0.62, col = INK, font = 2)
})

# ------------------------------------------------ C3: the composition link ----
## Within clade and size-free GC (x) against four responses (rows), per lineage
## (columns). Shading: log count; line: median in 20 GC bins; band: 25-75%.
C3_Y <- c(scf = "sCF mean", titv = "ti/tv", het_gain = "rate-het. dBIC/site", comphet = "compositional het.")
render("C3_gc_links_by_lineage", 13, 8.4, function() {
  par(mfrow = c(length(C3_Y), length(LINEAGES)), oma = c(3, 8.5, 3, 0.5))
  xl <- c(-2.2, 2.2); yl <- c(-2.5, 2.5)
  for (v in names(C3_Y)) for (L in LINEAGES) {
    ix <- d$lineage == L; x <- R[ix, "gc"]; y <- R[ix, v]
    base_par(mar = c(0.8, 1.8, if (v == names(C3_Y)[1]) 2 else 0.6, 0.4))
    plot(NA, xlim = xl, ylim = yl, axes = FALSE, xlab = "", ylab = "", xaxs = "i", yaxs = "i")
    nb <- if (sum(ix) >= 10000) 50 else 32
    bx <- seq(xl[1], xl[2], length.out = nb + 1); by <- seq(yl[1], yl[2], length.out = nb + 1)
    h <- table(cut(x, bx), cut(y, by)); z <- log10(unclass(h)); z[!is.finite(z)] <- NA
    if (any(!is.na(z))) image(head(bx, -1) + diff(bx) / 2, head(by, -1) + diff(by) / 2, z,
                              col = seq_ramp(40)[4:22], add = TRUE)
    g <- cut(x, quantile(x, 0:20 / 20), include.lowest = TRUE)
    xm <- tapply(x, g, median); q <- sapply(c(0.25, 0.5, 0.75), function(p) tapply(y, g, quantile, p))
    polygon(c(xm, rev(xm)), c(q[, 1], rev(q[, 3])), col = adjustcolor(SEQ[10], 0.18), border = NA)
    lines(xm, q[, 2], col = SEQ[11], lwd = 2)
    abline(h = 0, col = AXIS, lwd = 0.6)
    axis(2, at = c(-2, 0, 2), col = AXIS, labels = L == LINEAGES[1])
    if (v == tail(names(C3_Y), 1)) axis(1, at = c(-2, 0, 2), col = AXIS)
    mtext(sprintf("rho %.2f", cor(x, y)), side = 3, line = -1, adj = 0.04, cex = 0.5, col = INK2)
    if (v == names(C3_Y)[1]) mtext(sprintf("%s (%s)", L, format(sum(ix), big.mark = ",")), side = 3,
                                   line = 0.4, cex = 0.55, col = INK, font = 2)
    if (L == LINEAGES[1]) mtext(C3_Y[[v]], side = 2, line = 2, las = 1, adj = 1, cex = 0.6, col = INK2)
  }
  mtext("GC (model base frequencies), within clade and size-free (SD units)", side = 1, outer = TRUE,
        line = 1.8, cex = 0.62, col = INK2)
  mtext(paste("All axes: rank residual on clade + log tips + log sites, in SD units. Line: median in 20 GC bins;",
              "band: 25-75%; rho: correlation of the residuals."), side = 3, outer = TRUE, line = 1.4,
        adj = 0, cex = 0.58, col = MUTED)
})

# ---------------------------------------------------- C4: hub variables ----
## Hubs: a few variables, each standing for a group of strongly inter-correlated
## ones, and only weakly correlated with each other. On the within-clade,
## size-free correlations of all 19 variables:
##   1. factors: principal components with eigenvalue > 1 (Kaiser), varimax-rotated;
##   2. each variable joins the factor it loads on most, if |loading| >= 0.5;
##   3. a factor's hub is the member with the highest mean r^2 with the other
##      members (ties: highest mean r^2 with all variables).
## Coverage: R^2 of each variable on the hubs alone.
Cr <- cor(R)
ev <- eigen(Cr, symmetric = TRUE)
nf <- sum(ev$values > 1)
Lv <- unclass(varimax(ev$vectors[, 1:nf] %*% diag(sqrt(ev$values[1:nf])))$loadings)
rownames(Lv) <- colnames(R)
Lv <- Lv[, order(-colSums(Lv^2)), drop = FALSE]
Lv <- sweep(Lv, 2, sign(colSums(Lv)), `*`)
prim <- apply(abs(Lv), 1, which.max); load1 <- Lv[cbind(seq_len(nrow(Lv)), prim)]
member <- setNames(ifelse(abs(load1) >= 0.5, prim, NA), rownames(Lv))
C2_ <- Cr^2; diag(C2_) <- NA
overall <- rowMeans(C2_, na.rm = TRUE)
hubs <- unlist(lapply(seq_len(nf), function(f) {
  m <- names(which(member == f)); if (!length(m)) return(NULL)
  if (length(m) == 1) return(m)
  within <- sapply(m, function(v) mean(C2_[v, setdiff(m, v)]))
  m[order(-round(within, 3), -overall[m])][1]
}))
cover <- sapply(colnames(R), function(v)
  if (v %in% hubs) 1 else summary(lm(R[, v] ~ R[, hubs]))$r.squared)
hs <- data.frame(variable = colnames(R), label = VLAB[colnames(R)], block = VARS$block[match(colnames(R), VARS$key)],
                 factor = paste0("F", prim), loading = round(load1, 3), member = !is.na(member),
                 hub = colnames(R) %in% hubs, mean_r2_all = round(overall, 4),
                 coverage_r2_by_hubs = round(cover, 3), stringsAsFactors = FALSE)
hs <- hs[order(prim, -abs(load1)), ]
write.table(hs, file.path(OUT, "hub_selection.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
message("  hubs: ", paste(hubs, collapse = ", "), sprintf("  (%d factors)", nf))
render("C4_hub_variables", 13, 7.8, function() {
  par(oma = c(1, 0, 0, 0))
  layout(matrix(1:3, 1), widths = c(1, 0.1, 0.42))
  o <- hs$variable; n <- length(o); C <- Cr[o, o]
  base_par(mar = c(10.5, 10, 3, 0.5))
  image(1:n, 1:n, C[, n:1], col = div_ramp(21), breaks = seq(-1, 1, length.out = 22), axes = FALSE,
        xlab = "", ylab = "")
  abline(h = 0.5 + 0:n, v = 0.5 + 0:n, col = SURF, lwd = 1.5)
  ## factor blocks
  f <- hs$factor; st <- which(!duplicated(f)); en <- c(st[-1] - 1, n)
  rect(st - 0.5, n + 1 - en - 0.5, en + 0.5, n + 1 - st + 0.5, border = INK, lwd = 1.4)
  text((st + en) / 2, n + 1 - st + 0.5, unique(f), pos = 3, cex = 0.6, col = INK2, xpd = NA, offset = 0.2)
  lab <- paste0(hs$label, ifelse(hs$hub, "  ◆", ""))
  labcol <- ifelse(hs$block == "model", COL_MODEL, COL_METRIC)
  axis(1, at = 1:n, labels = FALSE, col = NA)
  text(1:n, -0.2, lab, srt = 90, adj = 1, xpd = NA, cex = 0.66, col = labcol, font = ifelse(hs$hub, 2, 1))
  text(0.2, n:1, lab, adj = 1, xpd = NA, cex = 0.66, col = labcol, font = ifelse(hs$hub, 2, 1))
  big <- which(abs(C) >= 0.4 & row(C) != col(C), arr.ind = TRUE)
  text(big[, 2], n + 1 - big[, 1], sprintf("%.2f", C[big]), cex = 0.48, col = ifelse(abs(C[big]) > 0.7, SURF, INK))
  title(main = "a  Within-clade, size-free correlations, grouped by rotated factor (◆ = hub)", adj = 0, line = 1.6)
  ## colour key
  base_par(mar = c(10.5, 0.5, 3, 3))
  image(1, seq(-1, 1, length.out = 21), matrix(seq(-1, 1, length.out = 21), 1), col = div_ramp(21),
        axes = FALSE, xlab = "", ylab = "")
  axis(4, at = c(-1, -0.5, 0, 0.5, 1), col = NA, col.axis = INK2)
  ## coverage, one bar per matrix row
  base_par(mar = c(10.5, 1, 3, 1))
  plot(NA, xlim = c(0, 1), ylim = c(0.5, n + 0.5), yaxs = "i", axes = FALSE, xlab = "", ylab = "")
  grid_v(seq(0, 1, 0.25))
  yy <- n:1
  rect(0, yy - 0.38, hs$coverage_r2_by_hubs, yy + 0.38, border = NA,
       col = ifelse(hs$hub, "#c3c2b7", labcol))
  axis(1, at = seq(0, 1, 0.25), col = AXIS)
  title(main = "b  R2 explained by the hubs", adj = 0, line = 1.6)
  title(xlab = "R2 (hubs: 1 by definition)", line = 2)
  mtext(sprintf("%d factors (eigenvalue > 1) after varimax hold %.0f%% of the variance. Blue labels: locus metrics; orange: model outputs.",
                nf, 100 * sum(ev$values[1:nf]) / ncol(R)),
        side = 1, outer = TRUE, line = -0.8, adj = 0.02, cex = 0.55, col = MUTED)
})
message("done")
