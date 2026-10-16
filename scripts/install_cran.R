#!/usr/bin/env Rscript
## Install the two packages conda-forge does not provide into the pixi
## environment's library. Run via `pixi run install-cran` (scripts/setup_env.sh).
##
## - PBSmodelling: current CRAN release (pure R + Tcl/Tk; needs BWidget, which
##   pixi.toml provides).
## - PhyInformR 1.0: archived from CRAN; pinned by URL and md5 so the build is
##   exactly the version the published results used.

lib <- .libPaths()[1]
repo <- "https://cloud.r-project.org"

if (!requireNamespace("PBSmodelling", quietly = TRUE))
  install.packages("PBSmodelling", lib = lib, repos = repo)

if (!requireNamespace("PhyInformR", quietly = TRUE)) {
  url <- paste0(repo, "/src/contrib/Archive/PhyInformR/PhyInformR_1.0.tar.gz")
  tgz <- file.path(tempdir(), basename(url))
  download.file(url, tgz, mode = "wb")
  md5 <- unname(tools::md5sum(tgz))
  if (md5 != "ed5150aaf92fb0c5f0b491332300b196") stop("PhyInformR tarball md5 mismatch: ", md5)
  install.packages(tgz, lib = lib, repos = NULL, type = "source")
}

for (p in c("ape", "phangorn", "TreeDist", "PBSmodelling", "PhyInformR"))
  cat(sprintf("%-14s %s\n", p, as.character(packageVersion(p))))
stopifnot(capabilities("cairo"))
