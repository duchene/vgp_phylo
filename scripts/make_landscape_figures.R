#!/usr/bin/env Rscript
## Figure sets B (along ancestral contigs) and D (deep vertebrate clades) of the
## best-model landscape (docs/PLAN, Step 1). Values are each locus's OWN best model
## (results/modelfinder/model_params.tsv); positions are the ancestral-genome
## coordinates in the locus names (results/stage0/locus_position.tsv).
##
##   B1  correlogram: how similar loci are as a function of distance along the same
##       ancestral contig, against a within-contig permutation null.
##   B2  tracks along the most locus-dense contigs of each clade.
##   B3  parameters against distance to the outermost locus of the contig.
##   D1  parameter distributions across the deep vertebrate clades.
##
## Only clades whose contigs carry many loci are used along the genome: Birds,
## Squamates, Cartilaginous fishes and VertebratesAnc5 (see docs/PLAN).
## Usage: Rscript scripts/make_landscape_figures.R

root <- Sys.getenv("VGP_PHYLO_ROOT", unset = normalizePath(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "..")))
source(file.path(root, "config", "project.R"))
FIG <- file.path(PATHS$figures, "model")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
set.seed(1)

# ------------------------------------------------------------------ palette ----
SURF   <- "#fcfcfb"; INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#898781"
GRID   <- "#e1e0d9"; AXIS <- "#c3c2b7"
CAT    <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100")   # adjacent use only
SEQ    <- c("#cde2fb","#b7d3f6","#9ec5f4","#86b6ef","#6da7ec","#5598e7",
            "#3987e5","#2a78d6","#256abf","#1c5cab","#184f95","#104281","#0d366b")
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

# --------------------------------------------------------------------- data ----
p   <- read.delim(file.path(PATHS$root, "results", "modelfinder", "model_params.tsv"), stringsAsFactors = FALSE)
pos <- read.delim(file.path(PATHS$root, "results", "stage0", "locus_position.tsv"), stringsAsFactors = FALSE)
d <- merge(p, pos[, c("locus_id", "contig", "mid", "n_loci_contig", "rank", "dist_to_outer")], by = "locus_id")
d$log_rate <- log10(d$mean_branch_length)

PARAM <- c(gc = "GC (model freqs)", titv = "ti/tv rate ratio", rate_cv = "FreeRate CV",
           k = "rate classes k", log_rate = "log10 tree length / branches")
ALONG <- c("BirdsAnc0", "SquamataAnc0", "CartilaginousFishesAnc0", "VertebratesAnc5")
SHORT <- c(BirdsAnc0 = "Birds", SquamataAnc0 = "Squamates", CartilaginousFishesAnc0 = "Cart. fishes",
           VertebratesAnc5 = "Vert. Anc5")
## Clade identity in line charts (all pairs visible): categorical slots 1-3, plus
## neutral ink for the 4th (the palette validates only 3 slots all-pairs).
CLADE_COL <- setNames(c("#2a78d6", "#eb6834", "#1baf7a", "#52514e"), ALONG)
MINLOCI <- 10L   # contigs with fewer loci carry no along-contig information

# ---------------------------------------------------------- B1: correlogram ----
## For loci on the same contig, the mean product of within-clade z-scored ranks of a
## parameter, by distance bin (a Moran's-I-like correlogram on ranks). Null: shuffle
## values among the loci of each contig (keeps every contig's mean and the positions,
## destroys order along the contig), 50 times; the band is its 2.5-97.5% range.
## Positive values above the band at short distance = local similarity beyond what
## the contig as a whole shares.
BINS <- c(0, 2e4, 5e4, 1e5, 2.5e5, 5e5, 1e6, 2e6, 6e6)
correlogram <- function(dd, v, nperm = 50) {
  z <- as.numeric(scale(rank(dd[[v]])))
  ctg <- split(seq_len(nrow(dd)), dd$contig)
  ctg <- ctg[lengths(ctg) >= MINLOCI]
  pairs <- do.call(rbind, lapply(ctg, function(ix) {
    cb <- utils::combn(ix, 2); cbind(cb[1, ], cb[2, ])
  }))
  bin <- cut(abs(dd$mid[pairs[, 1]] - dd$mid[pairs[, 2]]), BINS, include.lowest = TRUE)
  stat <- function(zz) tapply(zz[pairs[, 1]] * zz[pairs[, 2]], bin, mean)
  obs <- stat(z)
  null <- replicate(nperm, {
    zp <- z
    for (ix in ctg) zp[ix] <- z[ix][sample.int(length(ix))]
    stat(zp)
  })
  q <- function(x, pr) if (all(is.na(x))) NA_real_ else quantile(x, pr, na.rm = TRUE)   # empty bins
  np <- table(bin); thin <- as.vector(np) < 200   # too few pairs to show
  obs[thin] <- NA
  list(obs = obs, lo = ifelse(thin, NA, apply(null, 1, q, 0.025)),
       hi = ifelse(thin, NA, apply(null, 1, q, 0.975)), npairs = np)
}
cg <- lapply(setNames(ALONG, ALONG), function(cl) {
  message("  correlogram ", cl)
  dd <- d[d$clade == cl, ]
  lapply(setNames(names(PARAM), names(PARAM)), function(v) correlogram(dd, v))
})
binlab <- c("0-20k", "20-50k", "50-100k", "0.1-0.25M", "0.25-0.5M", "0.5-1M", "1-2M", "2-6M")
render("B1_correlogram_along_contigs", 13, 8, function() {
  par(mfrow = c(length(ALONG), length(PARAM)), oma = c(3.5, 9, 2.5, 0.5))
  for (cl in ALONG) for (v in names(PARAM)) {
    r <- cg[[cl]][[v]]
    base_par(mar = c(1.2, 2.6, 1.2, 0.6))
    yl <- range(c(r$obs, r$lo, r$hi, 0), na.rm = TRUE)
    plot(NA, xlim = c(0.5, length(r$obs) + 0.5), ylim = yl, axes = FALSE, xlab = "", ylab = "")
    grid_h(pretty(yl)); abline(h = 0, col = AXIS, lwd = 0.8)
    x <- seq_along(r$obs); ok <- !is.na(r$lo)
    polygon(c(x[ok], rev(x[ok])), c(r$lo[ok], rev(r$hi[ok])), col = "#e1e0d9", border = NA)
    lines(x, r$obs, col = SEQ[8], lwd = 2); points(x, r$obs, pch = 21, bg = SEQ[8], col = SURF, cex = 1.1)
    axis(2, at = pretty(yl), col = AXIS)
    if (cl == tail(ALONG, 1)) axis(1, at = x, labels = binlab, las = 2, col = AXIS, cex.axis = 0.62)
    if (cl == ALONG[1]) mtext(PARAM[[v]], side = 3, line = 0.6, cex = 0.62, col = INK, font = 2)
    if (v == names(PARAM)[1])
      mtext(sprintf("%s\n%s loci", SHORT[[cl]],
                    format(sum(d$clade == cl & d$n_loci_contig >= MINLOCI), big.mark = ",")),
            side = 2, line = 3.2, las = 1, adj = 1, cex = 0.62, col = INK2)
  }
  mtext("distance between loci on the same ancestral contig", side = 1, outer = TRUE, line = 2.4,
        cex = 0.65, col = INK2)
  mtext(paste("Correlogram on ranks (blue) vs within-contig permutation (grey band, 95%).",
              "Contigs with >= 10 loci; bins with < 200 locus pairs omitted."), side = 3, outer = TRUE, line = 1, adj = 0, cex = 0.62, col = MUTED)
})

# --------------------------------------------------------------- B2: tracks ----
## Per clade, the 3 contigs with the most loci; one row per parameter. Points are
## loci; the line is a running median over 7 neighbours. Best substitution model
## shown as a coloured rug (adjacent categories only, so slots 1-4 are safe).
SUBST <- c("GTR", "TN", "HKY", "K2P")
runmed_k <- function(y, k = 7) if (length(y) >= k) stats::runmed(y, k, endrule = "median") else y
render("B2_tracks_along_contigs", 14, 11, function() {
  ctgs <- unlist(lapply(ALONG, function(cl) {
    tb <- sort(table(d$contig[d$clade == cl]), decreasing = TRUE)
    paste(cl, names(tb)[1:3])
  }))
  par(mfrow = c(length(PARAM) + 1, length(ctgs)), oma = c(3, 10, 3.5, 0.5))
  rows <- c(names(PARAM), "subst")
  for (rw in rows) for (cc in ctgs) {
    cl <- sub(" .*", "", cc); ct <- sub(".* ", "", cc)
    dd <- d[d$clade == cl & d$contig == ct, ]; dd <- dd[order(dd$mid), ]
    x <- dd$mid / 1e6
    base_par(mar = c(0.6, 1.8, 0.6, 0.4))
    if (rw == "subst") {
      plot(NA, xlim = range(x), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
      rect(x - diff(range(x)) / 400, 0.1, x + diff(range(x)) / 400, 0.9,
           col = CAT[match(dd$subst, SUBST)], border = NA)
      axis(1, col = AXIS, cex.axis = 0.6)
    } else {
      yl <- range(d[[rw]][d$clade %in% ALONG], finite = TRUE)
      yl <- quantile(d[[rw]][d$clade %in% ALONG], c(0.005, 0.995), na.rm = TRUE)
      plot(x, dd[[rw]], ylim = yl, pch = 16, cex = 0.45, col = "#86b6ef", axes = FALSE, xlab = "", ylab = "")
      grid_h(pretty(yl))
      lines(x, runmed_k(dd[[rw]]), col = SEQ[10], lwd = 1.6)
      axis(2, at = pretty(yl), col = AXIS, cex.axis = 0.58)
    }
    if (rw == rows[1]) mtext(sprintf("%s\n%s (%d loci)", SHORT[[cl]], ct, nrow(dd)),
                             side = 3, line = 0.3, cex = 0.55, col = INK, font = 2)
    if (cc == ctgs[1]) mtext(if (rw == "subst") "best subst. model" else PARAM[[rw]],
                             side = 2, line = 2.6, las = 1, adj = 1, cex = 0.55, col = INK2)
  }
  mtext("position on the ancestral contig (Mb)", side = 1, outer = TRUE, line = 1.6, cex = 0.62, col = INK2)
  legend(grconvertX(0.72, "ndc", "user"), grconvertY(0.993, "ndc", "user"), xpd = NA, horiz = TRUE,
         bty = "n", legend = SUBST, fill = CAT, border = SURF, cex = 0.7, text.col = INK2, yjust = 0.5)
  mtext("Points: loci; line: running median of 7. Each locus shows its own best model's values.",
        side = 3, outer = TRUE, line = 2.4, adj = 0, cex = 0.62, col = MUTED)
})

# ------------------------------------------- B3: distance to outermost locus ----
render("B3_distance_to_contig_ends", 13, 3.6, function() {
  par(mfrow = c(1, length(PARAM)), oma = c(0, 0, 2, 7.5))
  brk <- c(0, 5e4, 1e5, 2.5e5, 5e5, 1e6, 6e6)
  bl <- c("<50k", "50-100k", "0.1-0.25M", "0.25-0.5M", "0.5-1M", ">1M")
  lc <- CLADE_COL[ALONG]
  for (v in names(PARAM)) {
    base_par(mar = c(4.8, 3, 2, 0.6))
    med <- sapply(ALONG, function(cl) {
      dd <- d[d$clade == cl & d$n_loci_contig >= MINLOCI, ]
      g <- cut(dd$dist_to_outer, brk, include.lowest = TRUE)
      md <- tapply(dd[[v]], g, median, na.rm = TRUE); md[table(g) < 50] <- NA; md   # thin bins
    })
    yl <- range(med, na.rm = TRUE)
    matplot(seq_along(bl), med, type = "o", lty = 1, lwd = 2, pch = 16, col = lc, ylim = yl,
            axes = FALSE, xlab = "", ylab = "")
    grid_h(pretty(yl))
    axis(1, at = seq_along(bl), labels = bl, las = 2, col = AXIS, cex.axis = 0.62); axis(2, col = AXIS)
    title(main = PARAM[[v]], adj = 0, line = 0.6)
  }
  legend(grconvertX(0.905, "ndc", "user"), grconvertY(0.6, "ndc", "user"), xpd = NA, bty = "n", legend = SHORT[ALONG],
         col = lc, lwd = 2, pch = 16, cex = 0.72, text.col = INK2)
  mtext("Median by distance from the locus to the outermost locus of its ancestral contig (contigs >= 10 loci; bins < 50 loci omitted)",
        side = 3, outer = TRUE, line = 0.4, adj = 0, cex = 0.62, col = MUTED)
})

# ----------------------------------------------- D1: deep vertebrate clades ----
DEEP <- sort(unique(d$clade[grepl("^VertebratesAnc", d$clade)]))
DEEP <- DEEP[table(d$clade)[DEEP] >= 20]   # VertebratesAnc0 has 2 loci
render("D1_deep_vertebrate_clades", 13, 7, function() {
  vv <- c(names(PARAM), "ntips_kept")
  lab <- c(PARAM, ntips_kept = "tips (after collapsing)")
  par(mfrow = c(2, 4), oma = c(0, 0, 2.5, 0))
  for (v in vv) {
    base_par(mar = c(5.5, 3.4, 2, 0.6))
    y <- split(d[[v]], factor(d$clade, levels = DEEP))
    yl <- quantile(unlist(y), c(0.005, 0.995), na.rm = TRUE)
    plot(NA, xlim = c(0.5, length(DEEP) + 0.5), ylim = yl, axes = FALSE, xlab = "", ylab = "")
    grid_h(pretty(yl))
    for (i in seq_along(y)) {
      q <- quantile(y[[i]], c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = TRUE)
      segments(i, q[1], i, q[5], col = SEQ[5], lwd = 1.5)
      rect(i - 0.28, q[2], i + 0.28, q[4], col = "#b7d3f6", border = SEQ[8], lwd = 1)
      segments(i - 0.28, q[3], i + 0.28, q[3], col = SEQ[12], lwd = 2.5)
    }
    axis(1, at = seq_along(DEEP), labels = sprintf("%s\n(%s)", sub("Vertebrates", "", DEEP),
         format(lengths(y), big.mark = ",")), col = AXIS, padj = 0.6, cex.axis = 0.62)
    axis(2, at = pretty(yl), col = AXIS)
    title(main = lab[[v]], adj = 0, line = 0.6)
  }
  ## last panel: substitution-model share
  base_par(mar = c(5.5, 3.4, 2, 0.6))
  tab <- prop.table(table(factor(d$clade, levels = DEEP), factor(d$subst, levels = SUBST)), 1)
  barplot(t(tab), col = CAT, border = SURF, names.arg = sub("Vertebrates", "", DEEP), axes = FALSE,
          cex.names = 0.62)
  axis(2, at = seq(0, 1, 0.25), labels = paste0(seq(0, 100, 25), "%"), col = AXIS)
  title(main = "best substitution model", adj = 0, line = 0.6)
  legend("bottom", inset = -0.32, xpd = NA, horiz = TRUE, bty = "n", legend = SUBST, fill = CAT,
         border = SURF, cex = 0.62, text.col = INK2)
  mtext(paste("Deep vertebrate clades (loci per clade in brackets; VertebratesAnc0 omitted, 2 loci).",
              "Box: 25-75%, line: median, whiskers: 10-90%."),
        side = 3, outer = TRUE, line = 0.8, adj = 0, cex = 0.62, col = MUTED)
})
message("done")
