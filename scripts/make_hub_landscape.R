#!/usr/bin/env Rscript
## Figure set E: the hub variables (results/links/hub_selection.tsv, from
## scripts/make_link_figures.R) along the ancestral genome, as set B does for the
## best-model parameters, plus their relation to reference-species chromosomes.
##
##   E1  correlogram along ancestral contigs vs a within-contig permutation null
##   E2  tracks along the most locus-dense contigs (rug: chicken chromosome class)
##   E3  against distance to the outermost locus of the contig
##   E4  by chromosome: Birds on chicken chromosomes, deep vertebrate clades on human
##       chromosomes, each against chromosome length
##
## Hub values are within clade and size-free (rank residual on clade + log tips +
## log sites, SD units; scripts/lib/link_data.R): rate per branch, for one, tracks
## taxon sampling (docs/PLAN, Step 1), and sampling itself varies along contigs.
## Soft-masked fraction is treated the same way; log tips is shown raw, as the
## check on that confound.
##
## Writes figures/model/E1..E4 (PDF + PNG). Usage: Rscript scripts/make_hub_landscape.R

root <- Sys.getenv("VGP_PHYLO_ROOT", unset = normalizePath(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "..")))
source(file.path(root, "config", "project.R"))
source(file.path(PATHS$root, "scripts", "lib", "link_data.R"))
FIG <- file.path(PATHS$figures, "model")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
set.seed(1)

# ------------------------------------------------------------------ palette ----
SURF   <- "#fcfcfb"; INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#898781"
GRID   <- "#e1e0d9"; AXIS <- "#c3c2b7"
CAT    <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100")   # adjacent use only, except slots 1-3
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
grid_v <- function(at) abline(v = at, col = GRID, lwd = 0.6)
source(file.path(PATHS$root, "scripts", "lib", "landscape.R"))

# --------------------------------------------------------------------- data ----
message("reading ...")
hs <- read.delim(file.path(PATHS$root, "results", "links", "hub_selection.tsv"), stringsAsFactors = FALSE)
HUBS <- hs$variable[hs$hub]
d <- load_link_data()
d <- d[complete.cases(d[, HUBS]), ]
for (v in HUBS) d[[v]] <- resid_clade_size(d[[v]], d)
has_sm <- any(is.finite(d$softmask))
if (has_sm) d$softmask <- resid_clade_size(d$softmask, d) else message("  no soft-mask table yet; skipped")
d$log_ntips <- log10(d$ntips)

PARAM <- c(setNames(paste(VLAB[HUBS], "(z)"), HUBS),
           if (has_sm) c(softmask = "soft-masked fraction (z)"), log_ntips = "log10 tips (raw)")

cg <- landscape_correlograms(d, PARAM)
draw_correlograms("E1_hubs_correlogram_along_contigs", d, PARAM, cg)
CLASS <- c("macro", "intermediate", "micro")
draw_tracks("E2_hubs_tracks_along_contigs", d, PARAM[names(PARAM) != "log_ntips"], rug = "chicken_class",
            rug_levels = CLASS, rug_label = "chicken chromosome",
            note = "Points: loci (within clade, size-free z); line: running median of 7. Rug: class of the chicken chromosome the locus lies on.")
draw_dist_ends("E3_hubs_distance_to_contig_ends", d, PARAM)

# ----------------------------------------------------- E4: by chromosome ----
## One point per reference chromosome: median z of its loci against log chromosome
## length; point area ~ loci. rho: Spearman across chromosomes, weighting ignored
## (chromosomes with < 30 loci dropped).
by_chr <- function(dd, chr, len, v) {
  ok <- !is.na(dd[[chr]]) & is.finite(dd[[v]])
  g <- split(dd[[v]][ok], dd[[chr]][ok]); first <- function(col) tapply(dd[[col]][ok], dd[[chr]][ok], `[`, 1)[names(g)]
  data.frame(chr = names(g), med = sapply(g, median), n = lengths(g), len = first(len),
             cls = if ("chicken_class" %in% names(dd)) first("chicken_class") else NA, stringsAsFactors = FALSE)
}
E4_V <- names(PARAM)
ROWS <- list(
  list(lab = "Birds\non chicken\nchromosomes", dd = d[d$clade == "BirdsAnc0", ], chr = "chicken_chr",
       len = "log_chicken_chr_len", by_class = TRUE),
  list(lab = "Deep vertebrate\nclades on human\nchromosomes", dd = d[d$lineage == "Deep vertebrates", ],
       chr = "human_chr", len = "log_human_chr_len", by_class = FALSE))
render("E4_hubs_by_chromosome", 13 * length(E4_V) / 6, 5.6, function() {
  par(mfrow = c(2, length(E4_V)), oma = c(1, 8.5, 2.4, 0.5))
  for (ri in seq_along(ROWS)) for (v in E4_V) {
    rw <- ROWS[[ri]]; t <- by_chr(rw$dd, rw$chr, rw$len, v); t <- t[t$n >= 30, ]
    base_par(mar = c(3.2, 2.6, if (ri == 1) 2 else 1, 0.6))
    yl <- range(c(t$med, 0)); yl <- yl + c(-1, 1) * 0.08 * diff(yl)
    plot(NA, xlim = range(t$len) + c(-0.1, 0.1), ylim = yl, axes = FALSE, xlab = "", ylab = "")
    grid_h(pretty(yl)); abline(h = 0, col = AXIS, lwd = 0.8)
    col <- if (rw$by_class) CAT[match(t$cls, CLASS)] else SEQ[8]
    symbols(t$len, t$med, circles = sqrt(t$n), inches = 0.09, add = TRUE, bg = adjustcolor(col, 0.75),
            fg = SURF)
    at <- c(0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200)
    at <- at[log10(at * 1e6) >= min(t$len) - 0.1 & log10(at * 1e6) <= max(t$len) + 0.1]
    axis(1, at = log10(at * 1e6), labels = at, col = AXIS); axis(2, at = pretty(yl), col = AXIS)
    title(xlab = "chromosome length (Mb, log scale)", line = 1.9, cex.lab = 0.65)
    if (ri == 1) title(main = PARAM[[v]], adj = 0, line = 0.6, cex.main = 0.72)
    mtext(sprintf("rho %.2f", cor(t$len, t$med, method = "spearman")), side = 3, line = -1, adj = 0.04,
          cex = 0.5, col = INK2)
    if (v == E4_V[1]) mtext(sprintf("%s\n%s loci", rw$lab, format(sum(!is.na(rw$dd[[rw$chr]])), big.mark = ",")),
                            side = 2, line = 2.4, las = 1, adj = 1, cex = 0.55, col = INK2)
  }
  legend(grconvertX(0.62, "ndc", "user"), grconvertY(0.995, "ndc", "user"), xpd = NA, horiz = TRUE, bty = "n",
         legend = paste("chicken", CLASS), pt.bg = CAT[1:3], pch = 21, col = SURF, pt.cex = 1.4, cex = 0.66,
         text.col = INK2, yjust = 1)
  mtext("One point per chromosome (>= 30 loci): median within-clade, size-free z of its loci; area ~ loci.",
        side = 3, outer = TRUE, line = 1, adj = 0, cex = 0.6, col = MUTED)
})
message("done")
