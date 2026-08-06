## Shared helpers: locus IO, label handling, and the metric registry.

suppressPackageStartupMessages(library(ape))

# ---------------------------------------------------------------- locus IO ----

#' Accession from a gene-tree / alignment label.
#' Labels are "<accession>.<chr>|copy<N>"; the accession carries its own version
#' dot (GCA_000001405.15) and <chr> may contain dots too, so split on the
#' GCA_/GCF_ pattern rather than on "." positions.
accession_from_label <- function(x) {
  m <- regexpr("^GC[AF]_[0-9]+\\.[0-9]+", x)
  out <- sub("\\|.*$", "", x)          # fallback for anything not accession-shaped
  out[m > 0] <- regmatches(x, m)      # regexpr+regmatches returns matches only
  out
}

#' Read a gene tree. IQ-TREE writes node labels as "/<support>" (e.g. "/0.734"),
#' which read.tree keeps as a character label -- strip the slash so support is
#' numeric, and stash it on the tree as $support.
read_gene_tree <- function(path) {
  tr <- ape::read.tree(path)
  if (!is.null(tr$node.label)) tr$support <- suppressWarnings(as.numeric(sub("^/", "", tr$node.label)))
  tr
}

#' Read a locus alignment as a DNAbin matrix.
read_locus_alignment <- function(path) {
  ape::read.dna(path, format = "fasta", as.matrix = TRUE)
}

#' Relabel a gene tree's tips to bare accessions, dropping any species that
#' appears more than once (defensive -- these loci are single-copy, verified in
#' the manifest, but a metric should not silently average over duplicates).
to_species_labels <- function(tr) {
  acc <- accession_from_label(tr$tip.label)
  dup <- acc %in% acc[duplicated(acc)]
  if (any(dup)) tr <- ape::drop.tip(tr, tr$tip.label[dup])
  tr$tip.label <- accession_from_label(tr$tip.label)
  tr
}

#' Load a reference species tree by key (see SPECIES_TREES in config/project.R).
#' Cached per process so repeated calls inside an mclapply worker are free.
.species_tree_cache <- new.env(parent = emptyenv())
load_species_tree <- function(key = "astralpro") {
  if (!exists(key, envir = .species_tree_cache)) {
    if (!key %in% names(SPECIES_TREES)) stop("unknown species tree: ", key)
    tr <- ape::read.tree(SPECIES_TREES[[key]])
    hit <- match(tr$tip.label, names(TIP_ALIASES))
    tr$tip.label[!is.na(hit)] <- TIP_ALIASES[hit[!is.na(hit)]]
    if (anyDuplicated(tr$tip.label)) stop("TIP_ALIASES created duplicate tips in ", key)
    assign(key, tr, envir = .species_tree_cache)
  }
  get(key, envir = .species_tree_cache)
}

# ------------------------------------------------------------ metric registry ----

.metric_registry <- new.env(parent = emptyenv())

#' Declare a metric. Called at the top level of a file in one of the METRIC_PATH
#' directories. See metrics/README.md for the full contract.
#'
#' @param name        identifier; also the output filename (results/metrics/<name>.tsv)
#' @param needs       which inputs to load: any of "tree", "alignment", "species_tree"
#' @param columns     the exact scalar columns `compute` returns, in order.
#'                    Required -- it fixes the output schema so a failed locus can
#'                    still emit a correctly-shaped all-NA row.
#' @param compute     function(loc) returning a named list/vector of scalars
#' @param description one line, shown by `run_metrics.sh --list`
#' @param species_tree which reference tree to attach when needs includes "species_tree"
register_metric <- function(name, needs, columns, compute, description = "",
                            species_tree = "astralpro") {
  stopifnot(is.character(name), length(name) == 1L, nzchar(name))
  stopifnot(is.character(columns), length(columns) > 0L, !anyDuplicated(columns))
  stopifnot(is.function(compute))
  bad <- setdiff(needs, c("tree", "alignment", "species_tree"))
  if (length(bad)) stop("metric '", name, "' declares unknown needs: ", paste(bad, collapse = ", "))
  if (any(columns %in% c("locus_id", "clade", "error")))
    stop("metric '", name, "': columns must not shadow locus_id/clade/error")
  assign(name, list(name = name, needs = needs, columns = columns, compute = compute,
                    description = description, species_tree = species_tree),
         envir = .metric_registry)
  invisible(NULL)
}

#' Locate a metric definition file on METRIC_PATH.
find_metric_file <- function(name) {
  for (d in METRIC_PATH) {
    p <- file.path(d, paste0(name, ".R"))
    if (file.exists(p)) return(p)
  }
  stop("no definition for metric '", name, "' on METRIC_PATH:\n  ",
       paste(METRIC_PATH, collapse = "\n  "))
}

#' Source a metric file and return the registered metric object.
load_metric <- function(name) {
  f <- find_metric_file(name)
  sys.source(f, envir = globalenv())
  if (!exists(name, envir = .metric_registry))
    stop("file ", f, " did not register_metric(name = '", name, "')")
  m <- get(name, envir = .metric_registry)
  m$file <- f
  m
}

#' All metric names discoverable on METRIC_PATH (files starting with "_" are
#' templates/helpers and are skipped).
available_metrics <- function() {
  fs <- unlist(lapply(METRIC_PATH, function(d)
    if (dir.exists(d)) list.files(d, pattern = "\\.R$", full.names = TRUE) else character()))
  nm <- sub("\\.R$", "", basename(fs))
  unique(nm[!startsWith(nm, "_")])
}

# ---------------------------------------------------------------- manifest ----

read_manifest <- function(path = PATHS$manifest) {
  if (!file.exists(path)) stop("manifest missing: ", path, " -- run scripts/build_manifest.sh")
  m <- utils::read.delim(path, stringsAsFactors = FALSE)
  m$fasta    <- file.path(PATHS$loci, m$clade, paste0(m$locus_id, ".fa"))
  m$treefile <- file.path(PATHS$gene_trees, m$clade, paste0(m$locus_id, ".fa.treefile"))
  m
}
