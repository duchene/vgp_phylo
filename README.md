# vgp_phylo

Per-locus phylogenetic metrics across 112,199 single-copy vertebrate loci.

**Location.** The repository lives at `/faststorage/project/DeiC-KU-L04/vgp_phylo` on
GenomeDK (moved from a 32-core workstation in October 2026). **All of its files —
input data, results, logs, job scripts and the R environment — stay inside this
directory**, not under the project's `<username>/` folders (an exception recorded in
the project manual, `../README.md`). Every path below is relative to the repo root.

## Setting up a fresh checkout

```bash
sbatch scripts/restore_data.sh   # data/ + working results/ from ../vgp_phylo_untracked.tar.gz (~8 min)
sbatch scripts/setup_env.sh      # R environment in .pixi/ from pixi.toml + pixi.lock (~10 min)
```

Both run on a compute node (no heavy I/O or installs on the login node). The restore
uses `--skip-old-files`, so it never overwrites anything already in the checkout —
in particular the tracked, tip-label-fixed `data/species_trees/astralpro.l6.tre`.
Afterwards, run any R script through the environment:

```bash
export PATH=/faststorage/project/DeiC-KU-L04/shared_software/pixi/current/bin:$PATH
pixi run Rscript scripts/compare_panel.R
```

## Layout

```
data/loci/<clade>/*.fa                 112,199 alignments
data/gene_trees/<clade>/*.treefile     112,199 IQ-TREE gene trees (1:1 with loci, verified)
data/species_trees/                    astralpro.l6.tre, wastral.merged.nwk (564 tips each)
data/metadata/                         astralinput.l6.mc.trees.namemap
data/manifest.tsv                      locus_id, clade, ntips, nspecies, ncopies_max, nsites
data/panel_187.tsv                     fixed 187-locus pilot panel (scripts/make_panel.R)
pixi.toml, pixi.lock                   R environment (installed into .pixi/, untracked)
config/project.R                       all paths, defaults, IQ-TREE settings
metrics/                               one file per metric  -- see metrics/README.md
scripts/                               drivers; scripts/lib/ holds shared helpers
scripts/submit_metrics.sh              Slurm submission: one job array per metric + gather
scripts/slurm/                         array task and gather bodies
scripts/submit_modelfinder.sh          ModelFinder screen (Slurm array); dedup_loci.R, parse_modelfinder.R
results/metrics/<metric>.tsv           one output file per metric
results/scf_tables/<tag>_m<model>_q<n>/  cached IQ-TREE runs (the dominant cost)
results/modelfinder/                   per-locus ModelFinder outputs + best_model.tsv, all_models.tsv.gz
results/locus_metrics.tsv              merged wide table (112,199 x 64), full precision, LOCAL
results/locus_quality.tsv              per-locus composite quality + ranks, LOCAL
results/published/                     the committed result set -- see below
scripts/make_figures.R                 summary figures for locus selection
figures/fig1..fig5                     PDF (publication) + PNG (viewing)
logs/
```

Loci are grouped by the **ancestral reference genome** they were called against
(22 clades, from `PetromyzontidaeAnc0` with 2 loci to `BirdsAnc0` with 47,126).
Alignment names encode provenance: `<Anc>_<Anc>refChr<N>_<start>_<end>.fa`.

Rebuild the manifest after any change to `data/loci`:

```bash
scripts/build_manifest.sh 32
```

## Dataset facts worth knowing

- **Strictly single-copy.** `ncopies_max == 1` and `ntips == nspecies` for all
  112,199 loci. The `|copyN` suffix in sequence labels is an artifact of the
  ASTRAL-Pro pipeline, not within-locus duplication.
- Sequence labels are `<accession>.<chr>|copy<N>`. The accession carries its own
  version dot (`GCA_000001405.15`) and `<chr>` may contain dots too, so parse with
  `accession_from_label()` rather than splitting on `.`.
- **8–557 tips, 1200–2500 sites** (239.6 Mb total). This range is why every metric
  needs an explicit comparability argument — see `metrics/README.md`.
- Alignments are soft-masked (mixed case) and contain gaps and `N`.
- Gene-tree node labels are IQ-TREE's `/<support>` form; `read_gene_tree()` strips
  the slash into a numeric `$support`.
- `data/metadata/singlecopy_gene.emptyfile` is a zero-byte file that came with the
  archive, preserved rather than deleted. It carries no data.

## Species tree tip label fix

`astralpro.l6.tre` originally labelled the human assembly **`GCA_000001405.5`**,
whereas the loci, the gene trees, the namemap and `wastral.merged.nwk` all use
**`GCA_000001405.15`** — the same genome under a different assembly-version suffix.

Left alone, human silently dropped out of every gene-vs-species comparison against
`astralpro`. The label was therefore changed **in the file** to `GCA_000001405.15`
(single occurrence; original preserved as `astralpro.l6.tre.orig`).

This was the only tip differing between the two reference trees. After the fix both
have identical 564-tip sets, and every namemap accession is present in both.
`TIP_ALIASES` in `config/project.R` is now empty but kept as a hook.

## Running metrics

All compute goes through Slurm on the project account (`DeiC-KU-L04`). Each metric
is one **job array**, one task per shard of loci, followed by a dependent **gather**
job that checks every shard is present, headers agree and no locus is missing or
duplicated, then writes `<outdir>/<metric>.tsv` atomically.

```bash
scripts/submit_metrics.sh --all --merge            # production: all loci, 2000 per task
scripts/submit_metrics.sh scf_own titv             # some metrics only
scripts/submit_metrics.sh --all --loci data/panel_187.tsv   # -> results/panel_187/
scripts/submit_metrics.sh --all --dry-run          # print the sbatch commands only
squeue --me                                        # progress; logs/m_<metric>_<array>_<task>.out
```

Each task writes its own shard file (`<outdir>/shards/<metric>/shard_NNN.tsv`), so
tasks never share an output file. Shards are round-robin over the manifest, and
`run_metric.R` resumes within a shard at block granularity, so resubmitting a failed
task (or the whole array) only redoes the missing loci. `scf_own` and `subst_model`
share one cached IQ-TREE run per locus; when both are submitted, `subst_model` waits
for `scf_own` so that work is not duplicated. Cached runs replay at ~440 loci/s.

Defaults are 4 cores, 8 GB and 4 h per task (`--cores`, `--mem`, `--time`).
IQ-TREE likelihood here is memory-bandwidth-bound: on the old workstation 12 -> 32
workers raised throughput only 1.16x. Prefer more array tasks over more cores per
task; going past ~8 cores per task mostly wastes the shared billing hours.

Inside an interactive allocation (`srun --account=DeiC-KU-L04 ...`) the original
drivers still work directly:

```bash
pixi run scripts/run_metrics.sh --list                     # what's registered
pixi run Rscript scripts/run_metric.R --metric scf_own --cores 4 --clade AnuraAnc1
pixi run Rscript scripts/merge_metrics.R
```

`VGP_PHYLO_ROOT` overrides the repo root (default: the checkout containing
`config/project.R`), `IQTREE` the binary, and `VGP_SCF_CACHE` the IQ-TREE cache
directory — used for the reproducibility check below, so the panel gets fresh
IQ-TREE runs instead of replaying the production cache.

`run_full.sh`, `run_detached.sh`, `run_status.sh` and `finish_run.sh` are the
workstation-era drivers (`setsid`/`nohup`, a 32-core budget). They are superseded
by `submit_metrics.sh` and must not be used on the GenomeDK login nodes.

## ModelFinder screen (model selection only)

Per-locus BIC model selection on the fixed gene tree: JC/F81/K2P/HKY/TN/GTR ×
{E, +R2..15, +H2..15}, after collapsing near-identical sequences. Design, pilot and
cost are in `docs/PLAN_genome_landscapes.md` (Step 0).

```bash
scripts/submit_modelfinder.sh --loci data/panel_187.tsv --shards 12   # -> results/modelfinder_panel_187/
scripts/submit_modelfinder.sh --shards 200 --cores 8 --mem 8G --time 06:00:00   # all loci -> results/modelfinder/
pixi run Rscript scripts/parse_modelfinder.R results/modelfinder 8    # best_model.tsv, all_models.tsv.gz
scripts/modelfinder_status.sh [--full] [results/modelfinder]          # progress, failures, CPU-hours, ETA
sbatch --account=DeiC-KU-L04 --dependency=afterany:<array> --cpus-per-task=8 --mem=16G --time=06:00:00 \
  --output=logs/mf_finalize_%j.out scripts/slurm/modelfinder_finalize.sh results/modelfinder   # parse + pack
```

`modelfinder_finalize.sh` parses the run, then packs each clade's per-locus files into
`results/modelfinder/archives/<clade>.tar` (~450k files -> 22), deleting a clade directory
only after its tar holds all 4 files per locus. It packs nothing if any locus is missing.
Read one locus without unpacking:

```bash
tar -xOf results/modelfinder/archives/BirdsAnc0.tar BirdsAnc0/<locus>.treefile.gz | zcat
```

## Best-model landscape figures

```bash
pixi run Rscript scripts/model_summary.R results/modelfinder/best_model.tsv results/modelfinder/model_params.tsv
pixi run Rscript scripts/stage0_positions.R          # results/stage0/locus_position.tsv (ancestral coordinates)
pixi run Rscript scripts/make_model_figures.R        # figures/model/A1-A3
pixi run Rscript scripts/make_landscape_figures.R    # figures/model/B1-B3, D1 (run on a compute node, ~16 GB)
```

What each figure shows, and the first observations, are in `docs/PLAN_genome_landscapes.md` (Step 1).

## Metric <-> model links, hub variables and path analysis

```bash
sbatch --account=DeiC-KU-L04 --cpus-per-task=8 --mem=4G --time=02:00:00 \
  --output=logs/stage0_softmask_%j.out scripts/slurm/stage0_softmask.sh   # results/stage0/locus_softmask.tsv (~25 min)
pixi run Rscript scripts/stage0_chromosomes.R                             # results/stage0/locus_chromosome.tsv
sbatch --account=DeiC-KU-L04 --cpus-per-task=2 --mem=24G --time=01:00:00 \
  --output=logs/link_figures_%j.out scripts/slurm/link_figures.sh         # C1-C4, E1-E4, P1-P2 (~20 min)
```

| Figure (`figures/model/`) | Shows |
|---|---|
| `C1_variance_partitioning` | per model output (and per metric): variance explained by clade, data size, the other block |
| `C2_rda_correlation_circle` | redundancy analysis of model outputs on metrics, within clade; repeated per lineage |
| `C3_gc_links_by_lineage` | GC against sCF, ti/tv, rate-heterogeneity gain, compositional heterogeneity |
| `C4_hub_variables` | the 6 hub variables chosen from 19 (rotated factors), and how much of each variable they cover |
| `E1`-`E3` | hubs along ancestral contigs, as B1-B3 |
| `E4_hubs_by_chromosome` | hubs per chicken chromosome (Birds) and per human chromosome (deep clades) against chromosome length |
| `P1_path_diagrams`, `P2_path_effects` | recursive path model: genomic predictors -> hubs, pooled and per lineage, contig-bootstrap CIs |

Shared code: `scripts/lib/link_data.R` (variables, transforms, within-clade residuals) and
`scripts/lib/landscape.R` (the along-contig figures, also used by `make_landscape_figures.R`).
Tables go to `results/links/` (local). Design and first observations: `docs/PLAN_genome_landscapes.md`,
Step 1b.

## Assessing which loci are best

```bash
Rscript scripts/make_figures.R      # ~1.5 min; writes figures/ and results/locus_quality.tsv
```

| Figure | Shows |
|---|---|
| `fig1_clade_confound` | variance in each metric explained by clade — **read this first** |
| `fig2_metric_redundancy` | clustered Spearman correlations; which metrics are near-duplicates |
| `fig3_pca_axes` | scree + loadings on PC1–PC4 |
| `fig4_tradeoff` | the divergence trade-off within clade, per lineage: branch length and informative sites vs sCF, QIRP vs branch length (other pairs logged) |
| `fig5_locus_quality` | composite score, by lineage, standardisation check, selection |

**The governing result: clade explains 21–69% of each metric, and R² = 0.71 for PC1.**
Loci therefore cannot be ranked globally — a global ranking mostly recovers which
ancestral reference the locus came from. Everything in `locus_quality.tsv` is
standardised **within clade** for this reason.

**The trade-off is real, not noise, and it survives within clade.** Long branches
buy informative sites (ρ = +0.62) and predicted power (QIRP, ρ = +0.70) but cost
per-site concordance (sCF, ρ = −0.72), pooled. Figure 4 shows them
within clade, per lineage (robust z per clade; CID as a within-clade percentile, since it is
discrete on small trees). There, branch length vs sCF is ρ = −0.55 (−0.46 after
also removing log tips and log sites), between −0.45 and −0.66 in every lineage.
Informative sites vs sCF is −0.45. Short branches give clean sites but unresolved
topologies. No single metric identifies "the best loci"; the composite is one
defensible compromise, not the answer.

**QIRP is mostly divergence, within clade.** Pooled, QIRP vs CID is ρ = −0.60;
within clade it is only −0.27 (−0.17 in Birds, sharks and the deep clades, about
−0.43 in Squamates and Amphibians). The rest was clade differences. Within clade
QIRP tracks mean branch length (+0.39), so it carries the divergence trade-off:
slightly lower sCF (−0.26) and slightly higher compositional heterogeneity
(+0.16). It is nearly independent of Holland δ (−0.10). Read it as expected
resolving power from the gene tree's own branch lengths, not as accuracy.
`make_figures.R` logs the within-clade ρ of every panel per lineage.

`locus_quality.tsv` columns: `quality` (mean of 8 signed within-clade robust
z-scores), `rank_within_clade`, `good_on_all_axes` (beats its clade median on all 8
— 3,111 loci), plus each component z-score. **Weights are equal by default; change
them in `scripts/make_figures.R`.**

⚠ One axis is circular for some purposes: `CID vs species tree` rewards agreement
with the ASTRAL tree. If the goal is to *infer* that tree, drop it from the
composite — the code makes this a one-line change.

## Published results (in this repository)

`results/published/` is committed, so the results are usable without re-running
anything:

| File | Size | Contents |
|---|---|---|
| `locus_metrics.tsv.gz` | 21 MB | all 112,199 loci x 64 metric columns |
| `locus_quality.tsv.gz` | 5.5 MB | composite quality, within-clade ranks, component z-scores |
| `data_dictionary.tsv` | 6 KB | every column: meaning, which direction is better, caveats |
| `summary_stats.tsv` | 6 KB | n / NAs / mean / sd / quantiles per column |

Read them straight from the gzip — no need to decompress:

```r
d <- read.delim(gzfile("results/published/locus_metrics.tsv.gz"))
```
```python
d = pd.read_csv("results/published/locus_metrics.tsv.gz", sep="\t")
```

**Start with `data_dictionary.tsv`.** Several columns are deliberately *not*
comparable across loci (`total_length`, `n_patterns`) or have no defensible "better"
direction (`scf_mean`, `entropy_mean`), and the dictionary says which.

Regenerate with `Rscript scripts/export_results.R` (~30 s). Two things it does:

- **Rounds to 6 significant digits.** R writes 15, which for concordance factors and
  probabilities is noise past the seventh. This halves the gzipped size (39 -> 21 MB)
  for a maximum relative change of 5x10^-6 (verified).
- **Separates published from working.** The full-precision `results/*.tsv` stay
  local, so regenerating them never produces a noisy diff.

⚠ These files sit in git history permanently. Re-committing a regenerated copy adds
another ~27 MB each time. If the metrics get revised repeatedly, attach the tables to
a GitHub Release instead and drop them from the tree.

## What is versioned, and what is not

Tracked: all code, config, both READMEs, and `data/species_trees/` (76 KB — it
records the tip-label fix above).

Not tracked, and why:

| Path | Size | How to restore |
|---|---|---|
| `data/loci/`, `data/gene_trees/` | 31 GB | `sbatch scripts/restore_data.sh` |
| `data/metadata/` | 3.4 MB | same archive |
| `data/manifest.tsv` | 8.4 MB | same archive, or `scripts/build_manifest.sh 32` (~16s) |
| `results/metrics/`, `results/*.tsv` | ~290 MB | same archive, or `scripts/submit_metrics.sh --all --merge` then `scripts/make_figures.R` |
| `results/scf_tables/` | 1.2 GB | same archive (cached IQ-TREE 3.1.3 runs, 112k files) |
| `results/panel_187/` | small | `scripts/submit_metrics.sh --all --loci data/panel_187.tsv` |
| `results/modelfinder/` (tables + `archives/`) | ~7 GB | `scripts/submit_modelfinder.sh --shards 200 --cores 8 --mem 8G --time 06:00:00` (~1,800 core-hours) |
| `results/modelfinder_panel_187{,_cmax5}/` | 12 MB | pilots: `scripts/submit_modelfinder.sh --loci data/panel_187.tsv --shards 12` |
| `.pixi/` | ~1 GB | `sbatch scripts/setup_env.sh` |
| `logs/` (except `run_full.log`) | — | regenerated per run |

Tracked results: `results/published/` (27 MB), `figures/` as PNG + PDF (1.2 MB), and
`logs/run_full.log` as the provenance record of the production run.

The archive is `../vgp_phylo_untracked.tar.gz` (7 GB; repo-relative paths). It is
the only copy of the input data on this cluster and nothing here is backed up —
keep the off-cluster original.

**The repository is self-contained.** Nothing is read from outside it at run time
except the IQ-TREE binary. Three helper functions that used to be sourced live from
the unversioned `~/scripts` (`pathnode`, `rtt.cov`, `stemmy`) are now vendored
verbatim in `scripts/lib/vendored.R`, with their upstream md5s recorded in the header.
`~/scripts` is untouched and remains the upstream — re-copy and update the md5 if you
improve one there.

Verified by making `~/scripts` unreadable: all 11 metrics still load and run, and the
three affected metrics reproduce their pre-vendoring values exactly.

## Reproducibility check after the move

Run on 2026-10-16 after moving to GenomeDK (R 4.4.3 via pixi, IQ-TREE 3.1.3 -> 3.1.4),
all 11 metrics recomputed on `data/panel_187.tsv` with a fresh IQ-TREE cache, and
compared with the original full-precision `results/metrics/`
(`scripts/compare_panel.R`; output in `results/panel_187/comparison.tsv`):

- **10 of 11 metrics are bit-for-bit identical**, including `subst_model` (κ and α
  re-estimated by IQ-TREE 3.1.4).
- **`scf_own` is not, and cannot be.** IQ-TREE is called without `-seed`, so
  `--scfl 100` draws different random quartets every run, and the old cache stores
  no seed. Two independent 3.1.4 runs differ from *each other* as much as from the
  3.1.3 original (scf_mean ρ = 0.992 between 3.1.4 replicates vs 0.997 against
  3.1.3; max |Δ| 0.040 vs 0.015). The difference is sampling noise, not the version.
  `n_branches` varies too because it counts branches with `sN > 0`.

The published results therefore stand.

**Seeding, from 2026-10-16 on.** Every IQ-TREE run now gets a deterministic per-locus
`-seed` (`locus_seed()` in `scripts/lib/iqtree.R`, from the locus_id and the base
`IQTREE_SEED`, default 1), and the base is part of the cache key
(`results/scf_tables/own_mHKYG4_q100_s1/`). Two independent seeded panel runs are
bit-for-bit identical. Set `IQTREE_SEED=""` to reproduce the legacy unseeded
behaviour and reuse the restored cache (`own_mHKYG4_q100/`), which is what the
published results came from.

## Dependencies

- R 4.4 with ape, phangorn, TreeDist and **PhyInformR 1.0**, all pinned in
  `pixi.toml` / `pixi.lock` (conda-forge). PhyInformR is archived from CRAN and its
  dependency PBSmodelling is not on conda-forge; `scripts/install_cran.R` builds
  both from source into the environment (PhyInformR checked by md5).
- IQ-TREE 3.1.4, the project's shared install
  (`shared_software/iqtree3/`) — the only external runtime dependency. Override with
  `IQTREE=/path/to/iqtree`. Requires `--seqtype DNA`; auto-detection fails on these
  alignments. The published results were computed with **3.1.3** on the old
  workstation; see "Reproducibility check after the move".

## Status

**Complete.** All 11 metrics computed for all 112,199 loci with **zero errors**;
merged table and summary figures written. Wall clock was ~9h on 12 cores, dominated
by `scf_own`. See `metrics/README.md` for per-metric detail and the remaining review
flags.

Dropped after review, each for a measured reason: `scf_species` and `qirp_species`
(needed branch lengths fitted to the ASTRAL topology, which would be a poor tree);
traditional RF, Nye similarity and MSID; and `phylomad_adequacy` — of which only
`unc_lnl_per_site` was retained, folded into `site_entropy`.
