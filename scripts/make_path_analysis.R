#!/usr/bin/env Rscript
## Path analysis linking the hub variables (results/links/hub_selection.tsv) to each
## other and to genomic predictors: a hypothesis about how they might be linked,
## written as a recursive (acyclic) set of linear regressions, not a test of it.
##
## Variables, all within clade and size-free (rank residual on clade + log tips +
## log sites, SD units, recomputed within each fitted set):
##   genomic predictors (exogenous, allowed to correlate freely):
##     log distance to the outermost locus of the ancestral contig, log span of the
##     contig's loci (a contig-size proxy), soft-masked (repeat) fraction, and, in
##     Birds, log length of the chicken chromosome the locus lies on;
##   hubs (endogenous), in the order of ORDER below: each hub is regressed on every
##     predictor and every hub before it.
## The model is saturated (every allowed path is estimated), so it fits perfectly by
## construction; what it gives is a decomposition of each correlation into direct and
## indirect paths UNDER THE ASSUMED ORDER. Reversing two adjacent hubs changes the
## direction of their arrow, not the fit.
##
## Fitted on loci of contigs with >= 10 loci (distance to the contig end is
## meaningless otherwise): pooled over the four clades of set B, and separately for
## Birds (with chicken chromosome length), Squamates and Cartilaginous fishes.
## Intervals: 2.5-97.5% of 200 bootstrap replicates resampling whole contigs within
## clade (neighbouring loci are not independent, figure B1).
##
## Writes results/links/path_coefficients.tsv and figures/model/P1, P2.
## Usage: Rscript scripts/make_path_analysis.R

root <- Sys.getenv("VGP_PHYLO_ROOT", unset = normalizePath(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "..")))
source(file.path(root, "config", "project.R"))
source(file.path(PATHS$root, "scripts", "lib", "link_data.R"))
FIG <- file.path(PATHS$figures, "model")
OUT <- file.path(PATHS$root, "results", "links")
set.seed(1)
NBOOT <- 200L

# ------------------------------------------------------------------ palette ----
SURF   <- "#fcfcfb"; INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#898781"
GRID   <- "#e1e0d9"; AXIS <- "#c3c2b7"
CAT    <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100")
DIV_LO <- "#2a78d6"; DIV_HI <- "#d03b3b"
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
grid_v <- function(at) abline(v = at, col = GRID, lwd = 0.6)

# ------------------------------------------------------------------ the model ----
## Assumed order, upstream to downstream: composition, then rate and the spread of
## rates over sites, then the signal those produce (informative sites, resolution
## probability), then gene-tree discordance with the species tree. Every variable
## of the hub selection is listed so that a different hub choice still gets an order.
ORDER <- c("gc", "titv", "rate", "rate_cv", "het_gain", "fast_rate", "k", "comphet", "entropy",
           "unclnl", "sn", "scf_sig", "stem", "rttcov", "qirp", "qihp", "scf", "delta", "cid")
EXO_LAB <- c(log_chicken_chr_len = "chicken chr. length", log_contig_span = "contig span",
             log_dist_end = "distance to contig end", softmask = "soft-masked fraction")
hs <- read.delim(file.path(OUT, "hub_selection.tsv"), stringsAsFactors = FALSE)
ENDO <- intersect(ORDER, hs$variable[hs$hub])
stopifnot(setequal(ENDO, hs$variable[hs$hub]))
LAB <- c(EXO_LAB, VLAB[ENDO])

message("reading ...")
d0 <- load_link_data()
if (!any(is.finite(d0$softmask))) stop("results/stage0/locus_softmask.tsv missing: run scripts/slurm/stage0_softmask.sh")
d0 <- d0[d0$clade %in% c("BirdsAnc0", "SquamataAnc0", "CartilaginousFishesAnc0", "VertebratesAnc5") &
         d0$n_loci_contig >= 10, ]

FITS <- list(
  pooled  = list(lab = "Pooled (4 clades, within clade)", clades = c("BirdsAnc0", "SquamataAnc0",
                 "CartilaginousFishesAnc0", "VertebratesAnc5"), exo = c("log_contig_span", "log_dist_end", "softmask")),
  birds   = list(lab = "Birds", clades = "BirdsAnc0",
                 exo = c("log_chicken_chr_len", "log_contig_span", "log_dist_end", "softmask")),
  squam   = list(lab = "Squamates", clades = "SquamataAnc0", exo = c("log_contig_span", "log_dist_end", "softmask")),
  cfish   = list(lab = "Cartilaginous fishes", clades = "CartilaginousFishesAnc0",
                 exo = c("log_contig_span", "log_dist_end", "softmask")))

## Direct paths B[to, from] by OLS on standardised variables; total effects
## (I - B)^-1 - I; indirect = total - direct.
fit_paths <- function(Z, exo) {
  v <- c(exo, ENDO); B <- matrix(0, length(v), length(v), dimnames = list(v, v))
  r2 <- setNames(numeric(length(ENDO)), ENDO)
  for (j in seq_along(ENDO)) {
    pr <- c(exo, ENDO[seq_len(j - 1)])
    f <- lm.fit(cbind(1, Z[, pr, drop = FALSE]), Z[, ENDO[j]])
    B[ENDO[j], pr] <- f$coefficients[-1]
    r2[ENDO[j]] <- 1 - sum(f$residuals^2) / sum((Z[, ENDO[j]] - mean(Z[, ENDO[j]]))^2)
  }
  Tt <- solve(diag(length(v)) - B) - diag(length(v)); dimnames(Tt) <- dimnames(B)
  list(B = B, T = Tt, r2 = r2)
}
flat <- function(f) {
  g <- expand.grid(to = rownames(f$B), from = colnames(f$B), stringsAsFactors = FALSE)
  g$direct <- f$B[cbind(g$to, g$from)]; g$total <- f$T[cbind(g$to, g$from)]
  g$indirect <- g$total - g$direct
  g
}

res <- list(); keep <- list()
for (nm in names(FITS)) {
  F <- FITS[[nm]]
  dd <- d0[d0$clade %in% F$clades, ]
  dd <- dd[complete.cases(dd[, c(F$exo, ENDO)]), ]
  Z <- sapply(c(F$exo, ENDO), function(v) resid_clade_size(dd[[v]], dd))
  est <- flat(fit_paths(Z, F$exo))
  ## contig bootstrap within clade
  grp <- split(seq_len(nrow(dd)), paste(dd$clade, dd$contig))
  gcl <- sub(" .*", "", names(grp))
  bs <- replicate(NBOOT, {
    pick <- unlist(lapply(split(seq_along(grp), gcl), function(ix) sample(ix, length(ix), replace = TRUE)))
    ix <- unlist(grp[pick], use.names = FALSE)
    b <- flat(fit_paths(scale(Z[ix, ]), F$exo)); cbind(b$direct, b$total, b$indirect)
  })
  q <- function(k, p) apply(bs[, k, ], 1, quantile, p)
  est$direct_lo <- q(1, 0.025); est$direct_hi <- q(1, 0.975)
  est$total_lo <- q(2, 0.025); est$total_hi <- q(2, 0.975)
  est$indirect_lo <- q(3, 0.025); est$indirect_hi <- q(3, 0.975)
  ft <- fit_paths(Z, F$exo)
  est <- est[est$total != 0 | est$direct_hi != 0, ]
  est <- est[est$to %in% ENDO & est$from != est$to, ]
  est$r2_to <- ft$r2[est$to]
  est$fit <- nm; est$n <- nrow(dd); est$n_contigs <- length(grp)
  res[[nm]] <- est; keep[[nm]] <- list(fit = ft, Z = Z, exo = F$exo, n = nrow(dd))
  message(sprintf("  %-7s %6d loci, %5d contigs; R2: %s", nm, nrow(dd), length(grp),
                  paste(sprintf("%s %.2f", ENDO, ft$r2), collapse = ", ")))
}
tab <- do.call(rbind, res)
num <- vapply(tab, is.numeric, TRUE); tab[num] <- lapply(tab[num], function(x) signif(x, 4))
write.table(tab[, c("fit", "from", "to", "direct", "direct_lo", "direct_hi", "indirect", "indirect_lo",
                    "indirect_hi", "total", "total_lo", "total_hi", "r2_to", "n", "n_contigs")],
            file.path(OUT, "path_coefficients.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ------------------------------------------------------------ P1: diagrams ----
## Hubs in a row in the assumed causal order (left to right); genomic predictors in a
## row above, with straight paths down to the hubs; hub-to-hub paths as arcs below
## the row, deeper for longer jumps. No path can pass through a node. Edges: direct
## paths with |beta| >= 0.1 whose bootstrap interval excludes 0; width ~ |beta|; red
## positive, blue negative. Under each hub: its R2 from everything upstream.
THR <- 0.1
diagram <- function(nm) {
  k <- keep[[nm]]; e <- res[[nm]]
  ne <- length(ENDO); nx <- length(k$exo)
  pos <- rbind(cbind(x = seq(0, 1, length.out = ne), y = 0),
               cbind(x = seq(0.5 - 0.17 * (nx - 1), 0.5 + 0.17 * (nx - 1), length.out = nx), y = 0.62))
  rownames(pos) <- c(ENDO, k$exo)
  base_par(mar = c(0.3, 0.3, 2, 0.3))
  plot(NA, xlim = c(-0.1, 1.1), ylim = c(-0.62, 0.72), axes = FALSE, xlab = "", ylab = "")
  sig <- e[abs(e$direct) >= THR & sign(e$direct_lo) == sign(e$direct_hi), ]
  sig <- sig[order(abs(sig$direct)), ]
  rx <- 0.055; ry <- 0.075                     # node half-width and half-height
  for (i in seq_len(nrow(sig))) {
    a <- pos[sig$from[i], ]; b <- pos[sig$to[i], ]
    col <- adjustcolor(if (sig$direct[i] > 0) DIV_HI else DIV_LO, 0.8); lwd <- 0.6 + 7 * abs(sig$direct[i])
    if (sig$from[i] %in% ENDO) {               # arc below the row
      depth <- 0.09 * (match(sig$to[i], ENDO) - match(sig$from[i], ENDO))
      t <- seq(pi, 0, length.out = 60)
      xs <- (a[1] + b[1]) / 2 + (b[1] - a[1]) / 2 * cos(t); ys <- -ry - depth * sin(t)
      lines(xs[1:57], ys[1:57], col = col, lwd = lwd)
      arrows(xs[57], ys[57], xs[60], ys[60], length = 0.07, col = col, lwd = lwd)
      text(mean(xs[29:32]), min(ys) - 0.025, sprintf("%.2f", sig$direct[i]), cex = 0.5, col = INK, font = 2)
    } else {                                   # straight down from a predictor
      dv <- b - a; edge <- function(px, py) min(px / abs(dv[1]), py / abs(dv[2]))   # where the line leaves a box
      f0 <- edge(rx, ry); f1 <- 1 - edge(rx + 0.01, ry + 0.01)
      arrows(a[1] + dv[1] * f0, a[2] + dv[2] * f0, a[1] + dv[1] * f1, a[2] + dv[2] * f1, length = 0.07,
             col = col, lwd = lwd)
      mp <- a + dv * 0.4
      text(mp[1], mp[2], sprintf("%.2f", sig$direct[i]), cex = 0.5, col = INK, font = 2)
    }
  }
  isE <- rownames(pos) %in% ENDO
  rect(pos[, 1] - rx, pos[, 2] - ry, pos[, 1] + rx, pos[, 2] + ry, col = ifelse(isE, "#f0efec", "#e1e0d9"),
       border = ifelse(isE, INK2, MUTED), lwd = 1.1)
  text(pos[, 1], pos[, 2] + ifelse(isE, 0.012, 0), sub("(.{8,}?) ", "\\1\n", LAB[rownames(pos)]),
       cex = 0.52, col = INK, xpd = NA)
  text(pos[isE, 1], pos[isE, 2] - 0.05, sprintf("R2 %.2f", k$fit$r2[ENDO]), cex = 0.45, col = MUTED)
  title(main = sprintf("%s  (%s loci)", FITS[[nm]]$lab, format(k$n, big.mark = ",")), adj = 0, line = 0.6)
}
render("P1_path_diagrams", 13, 10, function() {
  par(mfrow = c(2, 2), oma = c(1.6, 0, 1.6, 0))
  for (nm in names(FITS)) diagram(nm)
  mtext(sprintf(paste("Direct paths |beta| >= %.1f with a contig-bootstrap 95%% interval excluding 0; red +, blue -; width ~ |beta|.",
                      "Genomic predictors (top row, free to correlate) and hubs in the assumed causal order (left to right)."), THR),
        side = 3, outer = TRUE, line = 0.2, adj = 0.01, cex = 0.62, col = MUTED)
  mtext("All variables within clade and size-free (SD units). Saturated recursive model: a decomposition under the assumed order, not a test of it.",
        side = 1, outer = TRUE, line = 0.4, adj = 0.01, cex = 0.62, col = MUTED)
})

# ------------------------------------------------- P2: effects compared ----
## a: direct paths between hubs, per fit; b: total effect of each genomic predictor
## on each hub, split into its direct and indirect part (pooled), with the lineage
## fits' totals as points.
FCOL <- c(pooled = INK, birds = CAT[1], squam = CAT[2], cfish = CAT[3])
render("P2_path_effects", 13, 8.8, function() {
  par(oma = c(0, 0, 1, 0))
  layout(matrix(1:2, 1), widths = c(1, 1.15))
  ## a
  pe <- res$pooled; pe <- pe[pe$from %in% ENDO, ]; pe <- pe[order(abs(pe$direct)), ]
  pe <- tail(pe, 15)
  rows <- paste(pe$from, pe$to)
  base_par(mar = c(3.4, 14, 3, 1))
  xl <- range(unlist(lapply(res, function(e) { e <- e[paste(e$from, e$to) %in% rows, ]; c(e$direct_lo, e$direct_hi) })), 0)
  plot(NA, xlim = xl, ylim = c(0.5, length(rows) + 0.5), axes = FALSE, xlab = "", ylab = "", yaxs = "i")
  grid_v(pretty(xl)); abline(v = 0, col = AXIS)
  off <- setNames(seq(0.24, -0.24, length.out = length(FITS)), names(FITS))
  for (nm in names(FITS)) {
    e <- res[[nm]]; i <- match(rows, paste(e$from, e$to)); y <- seq_along(rows) + off[[nm]]
    segments(e$direct_lo[i], y, e$direct_hi[i], y, col = FCOL[[nm]], lwd = 1.2)
    points(e$direct[i], y, pch = 21, bg = FCOL[[nm]], col = SURF, cex = 0.9)
  }
  axis(1, at = pretty(xl), col = AXIS)
  axis(2, at = seq_along(rows), labels = paste(VLAB[pe$from], "\u2192", VLAB[pe$to]), col = NA,
       col.axis = INK2, cex.axis = 0.66)
  title(xlab = "direct path coefficient (SD/SD)", line = 2.1)
  title(main = "a  Direct paths between hubs (15 largest)", adj = 0, line = 1.9)
  legend("top", inset = c(0, -0.045), xpd = NA, horiz = TRUE, bty = "n", legend = sub(" \\(.*", "", sapply(FITS, `[[`, "lab")),
         pt.bg = FCOL, pch = 21, col = SURF, cex = 0.62, text.col = INK2)
  ## b
  pb <- res$pooled; pb <- pb[pb$from %in% names(EXO_LAB), ]
  pb <- pb[order(match(pb$from, names(EXO_LAB)), -match(pb$to, ENDO)), ]
  rows <- paste(pb$from, pb$to)
  bb <- res$birds; bb <- bb[bb$from == "log_chicken_chr_len", ]; bb <- bb[order(-match(bb$to, ENDO)), ]
  rows <- c(paste(bb$from, bb$to), rows)
  base_par(mar = c(3.4, 15, 3, 1))
  all <- do.call(rbind, res)
  xl <- range(all$total[paste(all$from, all$to) %in% rows], all$direct[paste(all$from, all$to) %in% rows], 0)
  xl <- xl + c(-1, 1) * 0.05 * diff(xl)
  plot(NA, xlim = xl, ylim = c(0.5, length(rows) + 0.5), axes = FALSE, xlab = "", ylab = "", yaxs = "i")
  grid_v(pretty(xl)); abline(v = 0, col = AXIS)
  for (r in seq_along(rows)) {
    src <- if (startsWith(rows[r], "log_chicken")) res$birds else res$pooled
    e <- src[paste(src$from, src$to) == rows[r], ]
    rect(0, r - 0.32, e$direct, r - 0.02, col = "#9ec5f4", border = NA)
    rect(0, r + 0.02, e$indirect, r + 0.32, col = "#c3c2b7", border = NA)
    for (nm in names(FITS)) {
      e2 <- res[[nm]][paste(res[[nm]]$from, res[[nm]]$to) == rows[r], ]
      if (nrow(e2)) {
        segments(e2$total_lo, r, e2$total_hi, r, col = FCOL[[nm]], lwd = 1)
        points(e2$total, r, pch = 21, bg = FCOL[[nm]], col = SURF, cex = 0.85)
      }
    }
  }
  fr <- sub(" .*", "", rows); to <- sub(".* ", "", rows)
  axis(2, at = seq_along(rows), labels = paste(EXO_LAB[fr], "\u2192", VLAB[to]), col = NA,
       col.axis = INK2, cex.axis = 0.62)
  br <- which(diff(match(fr, names(EXO_LAB))) != 0) + 0.5; abline(h = br, col = AXIS, lwd = 0.8)
  axis(1, at = pretty(xl), col = AXIS)
  title(xlab = "effect of the genomic predictor (SD/SD)", line = 2.1)
  title(main = "b  Genomic predictors: total (points) = direct + indirect", adj = 0, line = 1.9)
  legend("top", inset = c(0, -0.045), xpd = NA, horiz = TRUE, bty = "n",
         legend = c("direct (pooled; Birds for chr. length)", "indirect via hubs"), fill = c("#9ec5f4", "#c3c2b7"),
         border = SURF, cex = 0.62, text.col = INK2)
})
message("done")
