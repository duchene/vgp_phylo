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
## Squamates, Cartilaginous fishes and VertebratesAnc5 (see docs/PLAN). The B1-B3
## drawing code lives in scripts/lib/landscape.R.
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
SUBST <- c("GTR", "TN", "HKY", "K2P")

## B1-B3 are drawn by scripts/lib/landscape.R, shared with the hub-variable set E.
source(file.path(PATHS$root, "scripts", "lib", "landscape.R"))
cg <- landscape_correlograms(d, PARAM)
draw_correlograms("B1_correlogram_along_contigs", d, PARAM, cg)
draw_tracks("B2_tracks_along_contigs", d, PARAM, rug = "subst", rug_levels = SUBST,
            rug_label = "best subst. model",
            note = "Points: loci; line: running median of 7. Each locus shows its own best model's values.")
draw_dist_ends("B3_distance_to_contig_ends", d, PARAM)

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
