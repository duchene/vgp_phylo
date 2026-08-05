## Source every library file in a fixed order. Both run_metric.R and
## merge_metrics.R go through here so metrics can rely on all helpers existing.
local({
  d <- file.path(VGP_ROOT, "scripts", "lib")
  for (f in c("utils.R", "quartets.R", "iqtree.R", "treegeom.R")) sys.source(file.path(d, f), envir = globalenv())
})
