# vgp_phylo

Per-locus phylogenetic metrics across 112,199 single-copy vertebrate loci.

## Layout

```
data/loci/<clade>/*.fa                 112,199 alignments
data/gene_trees/<clade>/*.treefile     112,199 IQ-TREE gene trees (1:1 with loci, verified)
data/species_trees/                    astralpro.l6.tre, wastral.merged.nwk (564 tips each)
data/metadata/                         astralinput.l6.mc.trees.namemap
data/manifest.tsv                      locus_id, clade, ntips, nspecies, ncopies_max, nsites
config/project.R                       all paths, defaults, IQ-TREE settings
metrics/                               one file per metric  -- see metrics/README.md
scripts/                               drivers; scripts/lib/ holds shared helpers
scripts/run_full.sh                    production run: everything, two phases
results/metrics/<metric>.tsv           one output file per metric
results/scf_tables/<tag>_m<model>_q<n>/  cached IQ-TREE runs (the dominant cost)
results/locus_metrics.tsv              merged wide table (112,199 x 64)
results/locus_quality.tsv              per-locus composite quality + within-clade ranks
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

**Production run — use the detached launcher.** It survives SSH disconnects,
terminal closes, and any editor/agent session ending, because `setsid` puts it in a
new session and process group and `nohup` ignores SIGHUP:

```bash
scripts/run_detached.sh 12        # cores; writes logs/run_full.pid
scripts/run_status.sh             # progress, any time, from any shell
kill $(cat logs/run_full.pid)     # stop it
```

`run_full.sh` holds an `flock` single-instance lock, so a second launch refuses
rather than corrupting output — two runs appending to the same `results/metrics/*.tsv`
would interleave mid-line. Everything is resumable at block granularity, so if a run
does die, relaunching picks up where it stopped and cached IQ-TREE results replay at
~440 loci/s.

Run it in the foreground only if you will stay connected:

```bash
scripts/run_full.sh 12
```

Or piecemeal:

```bash
scripts/run_metrics.sh --list                      # what's registered
scripts/run_metrics.sh --all --jobs 8 --cores 4    # metrics in parallel, N cores each
Rscript scripts/run_metric.R --metric scf_own --cores 32 --clade AnuraAnc1
Rscript scripts/merge_metrics.R
```

Parallelism has two axes: `run_metric.R` forks over loci within one metric, and
`run_metrics.sh` runs several metrics as independent processes. Metrics share no
state and no output file, so adding one is free. Total load is roughly
`jobs x cores`.

**This machine is shared** and `nproc` is misleading. It reports 64, but the CPU is a
32-**physical**-core Threadripper PRO 5975WX with 2-way SMT; compute-bound work gains
little from an SMT sibling. With ~25 logged-in users and a measured ~6.5 cores of
other people's load (java/iqtree2, codeml), the usable budget is well under 32.

Note that a large share of the competing load is often **your own** other analyses
(`codeml`, `power_bootstrap`), not collaborators — check with `run_status.sh` and
`top` before assuming otherwise.

**More cores buys much less than you would expect.** IQ-TREE likelihood is
memory-bandwidth-bound: a 156-tip x 2500-site HKY+G4 partial-likelihood array is
~48 MB, already larger than one 32 MiB L3 instance (128 MiB total across 4
instances, 8 physical cores sharing each). Concurrent jobs thrash L3 and fall back to
DRAM. Measured: scaling 12 -> 32 workers raised throughput only **1.16x**
(size-normalised) — per-core efficiency fell 57% — and the *same* locus went from
1.0s to 5.3s. **16 or fewer cores is the sweet spot**; beyond that you mostly slow
down your colleagues. RAM is never the constraint (226 GB free; each worker holds one
alignment and tree).

Because `scf_own` is ~60% of the total work, give it everything first rather than
splitting evenly; see the recommended schedule in `metrics/README.md`.

## Assessing which loci are best

```bash
Rscript scripts/make_figures.R      # ~1.5 min; writes figures/ and results/locus_quality.tsv
```

| Figure | Shows |
|---|---|
| `fig1_clade_confound` | variance in each metric explained by clade — **read this first** |
| `fig2_metric_redundancy` | clustered Spearman correlations; which metrics are near-duplicates |
| `fig3_pca_axes` | scree + loadings on PC1–PC4 |
| `fig4_tradeoff` | six density panels showing metrics disagree about what "good" means |
| `fig5_locus_quality` | composite score, by lineage, standardisation check, selection |

**The governing result: clade explains 21–69% of each metric, and R² = 0.71 for PC1.**
Loci therefore cannot be ranked globally — a global ranking mostly recovers which
ancestral reference the locus came from. Everything in `locus_quality.tsv` is
standardised **within clade** for this reason.

**The trade-off is real, not noise.** Long branches buy informative sites
(ρ = +0.62) and predicted power (QIRP, ρ = +0.70) but cost per-site concordance
(sCF, ρ = −0.72). Short branches give clean sites but unresolved topologies
(higher CID). No single metric identifies "the best loci"; the composite is one
defensible compromise, not the answer.

`locus_quality.tsv` columns: `quality` (mean of 8 signed within-clade robust
z-scores), `rank_within_clade`, `good_on_all_axes` (beats its clade median on all 8
— 3,111 loci), plus each component z-score. **Weights are equal by default; change
them in `scripts/make_figures.R`.**

⚠ One axis is circular for some purposes: `CID vs species tree` rewards agreement
with the ASTRAL tree. If the goal is to *infer* that tree, drop it from the
composite — the code makes this a one-line change.

## What is versioned, and what is not

Tracked: all code, config, both READMEs, and `data/species_trees/` (76 KB — it
records the tip-label fix above).

Not tracked, and why:

| Path | Size | How to restore |
|---|---|---|
| `data/loci/`, `data/gene_trees/` | 31 GB | extract `singlecopy.tar.gz` |
| `data/metadata/` | 3.4 MB | ships with the archive |
| `data/manifest.tsv` | 8.4 MB | `scripts/build_manifest.sh 32` (~16s) |
| `results/` | ~210 MB | `scripts/run_detached.sh 12` then `scripts/make_figures.R` |
| `logs/` | — | regenerated per run |

**The repository is self-contained.** Nothing is read from outside it at run time
except the IQ-TREE binary. Three helper functions that used to be sourced live from
the unversioned `~/scripts` (`pathnode`, `rtt.cov`, `stemmy`) are now vendored
verbatim in `scripts/lib/vendored.R`, with their upstream md5s recorded in the header.
`~/scripts` is untouched and remains the upstream — re-copy and update the md5 if you
improve one there.

Verified by making `~/scripts` unreadable: all 11 metrics still load and run, and the
three affected metrics reproduce their pre-vendoring values exactly.

## Dependencies

- R packages: ape, phangorn, TreeDist, phytools, seqinr, **PhyInformR 1.0**.
  PhyInformR is archived from CRAN; it was installed from the CRAN archive plus
  gplots, geiger, gridExtra, hexbin, PBSmodelling.
- `~/iqtree3` (IQ-TREE 3.1.3) — the only external runtime dependency. Override with
  `IQTREE=/path/to/iqtree`. Requires `--seqtype DNA`; auto-detection fails on these
  alignments.

## Status

**Complete.** All 11 metrics computed for all 112,199 loci with **zero errors**;
merged table and summary figures written. Wall clock was ~9h on 12 cores, dominated
by `scf_own`. See `metrics/README.md` for per-metric detail and the remaining review
flags.

Dropped after review, each for a measured reason: `scf_species` and `qirp_species`
(needed branch lengths fitted to the ASTRAL topology, which would be a poor tree);
traditional RF, Nye similarity and MSID; and `phylomad_adequacy` — of which only
`unc_lnl_per_site` was retained, folded into `site_entropy`.
