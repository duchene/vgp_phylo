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
results/locus_metrics.tsv              merged wide table
logs/  figures/
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

## What is versioned, and what is not

Tracked: all code, config, both READMEs, and `data/species_trees/` (76 KB — it
records the tip-label fix above).

Not tracked, and why:

| Path | Size | How to restore |
|---|---|---|
| `data/loci/`, `data/gene_trees/` | 31 GB | extract `singlecopy.tar.gz` |
| `data/metadata/` | 3.4 MB | ships with the archive |
| `data/manifest.tsv` | 8.4 MB | `scripts/build_manifest.sh 32` (~16s) |
| `results/` | ~86 MB and growing | `scripts/run_detached.sh 12` |
| `logs/` | — | regenerated per run |

⚠ **`~/scripts` is not under version control.** Four metrics source it via
`source_shared()` (`stemmy`, `rtt.cov`, `pathnode`), so an unversioned edit there
silently changes results here. Worth making it a repo of its own, given it is shared
across projects.

## Dependencies

- R packages: ape, phangorn, TreeDist, phytools, seqinr, **PhyInformR 1.0**.
  PhyInformR is archived from CRAN; it was installed from the CRAN archive plus
  gplots, geiger, gridExtra, hexbin, PBSmodelling.
- `~/iqtree3` (IQ-TREE 3.1.3). Requires `--seqtype DNA`; auto-detection fails on
  these alignments.
- `~/scripts` — reusable single-function scripts shared with other projects,
  loaded by name via `source_shared()`.

## Status

11 metrics, all running clean on a 187-locus sample spanning every clade (0 errors),
including trees up to 557 tips. Nothing has been run at full scale yet; the measured
budget is ~7h at 12 cores / ~3h at 32 for everything. See `metrics/README.md` for
per-metric detail, the budget table, and open review flags.

Dropped after review, each for a measured reason: `scf_species` and `qirp_species`
(needed branch lengths fitted to the ASTRAL topology, which would be a poor tree);
traditional RF, Nye similarity and MSID; and `phylomad_adequacy` — of which only
`unc_lnl_per_site` was retained, folded into `site_entropy`.
