# Per-locus metrics

One file per metric in this directory. Each calls `register_metric()` and is picked
up automatically — `scripts/run_metrics.sh --list` shows what's registered, and
metrics run as independent processes so adding one costs nothing.

```r
register_metric(
  name        = "my_metric",          # also the output file: results/metrics/my_metric.tsv
  needs       = c("tree","alignment"),# also "species_tree"; only these are loaded
  columns     = c("a","b"),           # fixes the schema, so a failed locus emits a shaped NA row
  compute     = function(loc) list(a = ..., b = ...),
  description = "one line"
)
```

`loc` carries `locus_id`, `clade`, `fasta`, `treefile`, `ntips`, `nsites`, plus
`tree` / `alignment` / `species_tree` if declared in `needs`.

---

## The comparability rule

Every metric must be comparable across **8–557 tips** and **1200–2500 sites**. That
constraint is not automatic, and it is the single most common way a standard
phylogenetic statistic goes wrong here. Four failure modes recur:

| Failure mode | Example | Treatment used |
|---|---|---|
| Scales with taxon count | RF max is 2(n−3); total tree length sums 2n−3 branches | normalise by the maximum, or report per-branch |
| Scales with site count | unconstrained lnL; any *t*-statistic grows as √n_sites | per-site, or a standardised effect size |
| Estimator bias depends on n | plug-in entropy is biased by ≈(k−1)/2n | Miller–Madow correction |
| Sampling effort depends on n | quartets number C(n,4): 70 at n=8, 4×10⁹ at n=557 | fixed sample of 2000 quartets |

Where a raw value is kept alongside a corrected one, it is **only** so the
correction can be audited. Columns named `*_raw`, `total_length`, and `chisq_raw`
must not be compared across loci.

Two standing conventions:

- **No re-inference.** Gene-tree topology *and* branch lengths are taken as given
  from the existing IQ-TREE trees, held via `-blfix`. Only substitution-model
  parameters are estimated. Every metric requiring branch lengths fitted to a
  *different* topology has been dropped — see the coverage table.
- **`--scfl` requires a model**, since it is likelihood-based. Fixed at `HKY+G4`
  for all loci (never selected per locus, which would make sCF partly a function of
  which model won). Measured on 187 loci: HKY+G4 vs GTR+G4 gives ρ=0.9973 and a
  mean difference of 0.0003 at half the cost; JC gives ρ=0.9896 at a quarter.
- **One IQ-TREE run per locus**, cached in `results/scf_tables/`. The model and
  quartet count are part of the cache key, so changing either invalidates cleanly.

---

## Coverage of the requested list

| Requested | Implemented as | Note |
|---|---|---|
| sCF mean & variance, own tree | `scf_own` | |
| sCF mean & variance, species tree | — | **dropped**: needs branch lengths fitted to the ASTRAL topology, which would be a poor tree |
| QIRP mean & variance, own tree | `qirp_own` | per-internode, as you expected |
| QIRP mean & variance, species tree | — | **dropped**: same reason — QIRP is entirely a function of branch lengths |
| Holland's δ, reticulation | `delta_plot` | |
| Probabilistic RF to species tree, weighted & unweighted | `treedist_species` | now CID only; see below |
| CoV root-to-tip, own tree | `rtt_cov` | your `rtt.cov.R` |
| Stemminess, own tree | `stemminess` | your `stemmy.R` |
| Mean TiTv | `subst_model.kappa` (preferred) + `titv` | ML kappa is free from the sCF run; titv is the model-free version — see below |
| Aitchison compositional distances | `aitchison_basefreq` | base frequencies, per your correction |
| Site entropy mean & variance | `site_entropy` | |
| Other phylomad adequacy metrics, size-comparable | `site_entropy.unc_lnl_per_site` | only this one retained; the rest were measurably redundant — see below |
| — | `tree_length` | added: overall rate proxy from existing trees |
| — | `subst_model` | added: HKY kappa and gamma alpha, free from the sCF run |

---

# Alignment-only metrics

## `site_entropy` — site-pattern information

Shannon entropy of the ACGT distribution at each site, summarised over sites, plus
the unconstrained multinomial log-likelihood. Low mean entropy means conserved; high
mean entropy approaches the saturation expectation. For the entropy columns, gaps, N
and IUPAC codes are missing data and sites with <2 unambiguous bases are dropped.

| Column | Meaning |
|---|---|
| `entropy_mean`, `entropy_var` | plug-in entropy, natural log, over sites |
| `entropy_mean_mm`, `entropy_var_mm` | **use these** — Miller–Madow bias-corrected |
| `entropy_mean_norm` | mean ÷ log 4, so [0,1] = fraction of maximal disorder |
| `prop_invariant` | share of sites with one observed state |
| `unc_lnl_per_site` | unconstrained (multinomial) lnL ÷ n_sites — see below |
| `n_patterns` | distinct site patterns. ⚠ raw diagnostic, **not** comparable: rises with both taxon count and length |
| `nsites_used`, `ntips_used` | denominators |

**Comparability.** Entropy is bounded by log 4 regardless of n, but the *plug-in
estimator* is biased downward by ≈(k−1)/2n_obs. Since n runs 8→557, that bias alone
would manufacture a correlation between entropy and locus size — hence the
Miller–Madow columns.

### `unc_lnl_per_site` — the one retained PhyloMAd statistic

Σ n_i·log(n_i) over distinct site patterns, minus n_sites·log(n_sites), divided by
n_sites. This is the likelihood ceiling no substitution model can exceed, and
equivalently a measure of how repetitive the site patterns are.

Algebraically equivalent to `~/phylomad/codeFolder/testStatistics/get.unconstrained.R`
— **verified identical** (`all.equal` TRUE) on a 509×2500 locus — but reimplemented
inline so this metric carries no dependency on the phylomad checkout, and running
~1.6× faster (0.199s vs 0.321s). Computed on the full alignment, gaps and N included
as states, matching phylomad rather than the ACGT-filtered subset used for entropy.

Retained because it was the only PhyloMAd statistic carrying independent
information: ρ = +0.22 against `entropy_mean_mm` (here +0.34 on a smaller sample),
where the others hit −1.00 against `stemminess` and +0.99 against
`aitchison_basefreq`. Per-site because the raw sum scales with alignment length.

## `titv` — transition/transversion ratio

`mean(TS distance) / mean(TV distance)` over all sequence pairs, following your
`meanTiTv.R`, using `ape::dist.dna(pairwise.deletion = TRUE)`.

| Column | Meaning |
|---|---|
| `titv` | ratio of means |
| `ts_mean`, `tv_mean` | the two components |
| `pdist_mean` | mean raw p-distance (divergence context) |

**Why ratio-of-means, not mean-of-ratios.** The per-pair ratio is undefined whenever
a pair has zero transversions, which is common for close tips in short alignments.

**This is not κ**, and is not scaled for the 2:1 excess of possible transversion
types, so values below 1 are expected and do not mean transversions outnumber
transitions per opportunity. Observed range 0.78–1.94, median 1.25.

**Prefer `subst_model.kappa`.** It is the ML transition/transversion rate ratio,
obtained free from the run `scf_own` already pays for. See that section for the
measured comparison — the two agree at ρ = 0.834, but `titv` is materially more
saturation-biased. `titv` is retained as a model-free cross-check.

## `aitchison_basefreq` — compositional heterogeneity among taxa

Each sequence gives a 4-part composition (A,C,G,T proportions). These sum to 1, so
Euclidean distance between them is not meaningful — a gain in one base is
necessarily offset by others. Aitchison distance (Euclidean between centred
log-ratio transforms) depends only on the *ratios* between bases and is the correct
geometry. Measures the same signal as the χ² composition test, in proper
compositional space.

| Column | Meaning |
|---|---|
| `ait_mean`, `ait_var`, `ait_max` | over all sequence pairs |
| `ait_totvar` | compositional total variance = mean squared distance to the compositional centre; the standard one-number dispersion summary |
| `ait_to_centroid_mean` | mean distance to the centre |
| `gc_mean`, `gc_sd` | GC content, for interpretation |
| `nsites_used`, `ntips_used` | see below |

Zeros are smoothed with a **Jeffreys pseudocount on the raw base counts**, so a
2500-base sequence is barely perturbed while a short one is smoothed more. The usual
"fraction of the smallest positive part" replacement is wrong here — it blows up on
compositions like (1,0,0). Sequences with <20 unambiguous bases are dropped.

**Comparability.** Scale-invariant and dimensionless. Two second-order effects,
both reportable: composition noise falls as 1/√n_sites (small over 1200–2500), and
the number of pairs is C(n,2), so the *precision* of the mean varies with n even
though its expectation does not.

## `delta_plot` — Holland's δ (reticulation signal)

For a quartet {i,j,k,l} the three pairings give distance sums; sorted descending as
m₁≥m₂≥m₃, δ = (m₁−m₂)/(m₁−m₃). δ=0 is perfectly treelike, δ=1 fully "boxy". Averaged
over a fixed sample of quartets.

| Column | Meaning |
|---|---|
| `delta_mean`, `delta_var` | over sampled quartets |
| `delta_q90` | 90th percentile — catches a reticulate minority a mean would hide |
| `q_used`, `ntips_used` | quartets actually resolvable |

Distances are **TN93, not raw p-distance**: saturation compresses raw distances
non-linearly, which by itself manufactures four-point violations and would confound
reticulation with divergence. TN93 falls back to raw for the few pairs where it is
undefined.

**Comparability.** δ is a ratio of distance differences — invariant to overall rate
and alignment length. Fixed 2000-quartet sample keeps sampling effort (and so the
sampling variance of the mean) constant instead of varying with C(n,4).

Holland, Huber, Dress & Moulton (2002) *MBE* 19:2051–2059.

---

# Tree-only metrics

## `tree_length` — overall rate proxy

Read straight off the existing IQ-TREE gene tree; nothing re-estimated.

| Column | Meaning |
|---|---|
| `mean_branch_length` | **the rate proxy to use** |
| `median_branch_length` | robust alternative |
| `mean_root_to_tip` | mean tip depth, midpoint-rooted — per-lineage accumulated substitutions |
| `total_length` | ⚠ **not comparable across loci** — see below |
| `n_branches` | 2n−3 |
| `prop_zero_length` | share of ~zero branches = effectively unresolved topology |

**⚠ Comparability.** Total tree length grows roughly linearly with taxon count
(2n−3 branches), so a 557-tip tree out-totals an 8-tip tree at any rate. It is
reported only to audit the normalisation.

## `stemminess` — internal vs total tree length

Your `stemmy.R`: share of total tree length on internal branches.

| Column | Meaning |
|---|---|
| `stemminess` | raw proportion |
| `stemminess_per_branch` | **size-corrected** — divided by the fraction of branches that are internal |
| `n_internal`, `n_branches` | counts |

**Comparability.** The raw proportion is free of rate and length effects but *not*
of taxon count: the share of edges that are internal rises from 5/13 at n=8 toward
1/2 as n grows, so larger trees drift upward for purely combinatorial reasons. The
corrected column is mean internal branch length relative to mean branch length.

## `rtt_cov` — root-to-tip rate variation

Your `rtt.cov.R` + `pathnode.R`: CoV of root-to-tip distances, midpoint-rooted
(phangorn) when the tree is unrooted, which these are.

| Column | Meaning |
|---|---|
| `rtt_cov` | sd/mean of root-to-tip distances |
| `rtt_mean` | mean root-to-tip distance |
| `ntips_used` | |

**Comparability.** CoV is scale-free, so unaffected by rate or length. Midpoint
rooting depends only on the tree, so no outgroup is assumed and loci with different
taxon sets stay comparable. Note the *spread* of any CoV estimate depends on n.

---

# Metrics using IQ-TREE

`scf_own` uses `--scfl` (likelihood-based, IQ-TREE's recommendation) with a fixed
100 quartets per branch. Branches with `sN == 0` are **dropped, not counted as
zero** — no informative sites is missing data, and averaging them in would bias sCF
downward precisely for the shortest branches.

**Why not `--scf`.** The parsimony version needs no branch lengths or model at all
(verified: a topology-only tree gives 52.26 vs 52.21 with lengths) and costs half.
But it reads ~16 points lower than `--scfl` (52.2 vs 68.1) — its known downward
bias, since it only counts sites decisive under parsimony. Cost is not the binding
constraint here, so the less biased estimator wins.

## `scf_own` — site concordance with the locus's own tree

Internal consistency: how much of the site-level signal supports each branch the
gene tree itself asserts. Uses `-blfix`, so IQ-TREE estimates nothing at all.

| Column | Meaning |
|---|---|
| `scf_mean`, `scf_var` | over branches, as proportions in [0,1] |
| `sdf1_mean`, `sdf2_mean` | the two discordant resolutions |
| `sn_mean` | mean informative sites per branch |
| `n_branches` | branches with informative sites |
| `scf_prop_sig` | **proportion** of branches beating p=1/3 by binomial test at 0.001, per your `binom.scf.test.R` |

**Comparability.** sCF is a proportion, and quartet sampling per branch is fixed.
`scf_prop_sig` is a proportion rather than a count, so it does not scale with size.

## `qirp_own` — quartet internode resolution probability

Townsend, Su & Tekle (2012) via `PhyInformR::Approximator(t, t0, rates, 4)`, which
returns (QIRP, QIPP, QIHP) = probability the quartet spanning an internode is
resolved **correctly**, left **polytomous**, or resolved **incorrectly**
(homoplasy). Computed per internode, so summarised as mean and variance.

`t0` is the internode's own length; `t` is the mean tip depth below it (mean, not
max, because these trees are not ultrametric and QIRP concerns the typical quartet).
Trees are midpoint-rooted to get a consistent time axis.

| Column | Meaning |
|---|---|
| `qirp_mean` | mean predicted probability of **correct** resolution — the locus's predicted power |
| `qihp_mean` | mean predicted probability of **incorrect** resolution — a locus can be informative *and* actively misleading. No other metric in the set captures this |
| `qirp_var_rel` | dispersion of QIRP across internodes, ÷ its mechanical maximum m(1−m) — see below |
| `n_internodes` | internodes actually evaluated (capped, see below) |
| `n_internodes_total` | internodes available before capping |
| `rate_var` | rates are scaled to mean 1, so this **is** the locus's rate heterogeneity |
| `prop_invariant_rate` | share of zero-rate (invariant) sites |

Geometry is read from the existing gene tree — nothing is inferred, and no IQ-TREE
call is made.

**Capping (`QIRP_MAX_INTERNODES = 100`, `QIRP_MAX_RATES = 500`).** `Approximator`
dominates the runtime: profiled on a 509-tip / 2500-site locus it took 21.8s of
22.5s (507 internodes × 2500 rates = 1.27M evaluations), against 0.06s for the site
rates and 0.18s for node depths. Both dimensions are subsampled at random, which
cut the metric from 49.3s to 9.7s per 187 loci — a 5.1× speedup. This is also a
comparability fix on the same logic as the fixed quartet sample: uncapped, the mean
would be estimated from 5 internodes at n=8 and 554 at n=557, so its sampling
precision would scale with tree size. Set either constant to `Inf` to use everything.

**Site rates.** Per-site parsimony steps on the fixed tree, scaled to mean 1
(0.01s for 557 tips × 2190 sites). Scaling to mean 1 isolates rate *heterogeneity*
— what QIRP responds to — from overall rate, which enters separately through `t`
and `t0`. `as.phyDat` compresses to site patterns, so scores are expanded back
through the pattern index; otherwise the rate distribution would be reweighted by
pattern uniqueness.

**Comparability.** All three outputs are probabilities in [0,1].

**⚠ For your review — two things.** (1) Parsimony undercounts substitutions at fast
sites, so the rate distribution is mildly compressed relative to ML; a deliberate
cost trade, since ML site rates would mean fitting a rate model per locus.
(2) `qirp_mean ~ mean_branch_length` = +0.62 while `qirp_mean ~ rate_var` = only
+0.15. QIRP responds far more to branch lengths than to the rate heterogeneity it is
fed. That is intrinsic rather than a bug — with rates normalised to mean 1, the
scale of rate variation is absorbed into the branch lengths and only the *shape* of
the rate distribution remains — but it means QIRP is substantially a nonlinear
re-expression of `tree_length.mean_branch_length`. It is retained because +0.62 is
moderate rather than redundant, because it encodes internode length *relative to its
depth* (which a mean branch length cannot), and because `qihp_mean` has no substitute
in the set.

### Why only two probabilities, and a normalised variance

`Approximator` returns (QIRP, QIPP, QIHP), which **sum to exactly 1** — verified, max
deviation 1.1×10⁻¹⁵. It is a 3-part composition with only **two free dimensions**, so
reporting all three is redundant by construction; that constraint, not a coincidence,
is why `qirp_mean` correlated −0.86 with `qihp_mean`. `qipp_mean` is therefore dropped
as the remainder.

The raw `qirp_var` was also dropped. Variance of a [0,1] quantity is mechanically
bounded by m(1−m), so it tracks the mean for reasons carrying no information:
ρ(mean, var) = −0.83 raw, falling to **−0.54** once divided by that bound. The
normalised `qirp_var_rel` is reported instead.

## `subst_model` — fitted substitution-model parameters

Read off the IQ-TREE run `scf_own` already pays for, so this metric costs **no extra
IQ-TREE invocations** (verified: after `subst_model` populates the cache, `scf_own`
runs at 440 loci/s instead of 7).

| Column | Meaning |
|---|---|
| `kappa` | ML transition/transversion rate ratio under HKY. Observed 2.16–4.89 |
| `gamma_alpha` | +G4 gamma shape. Low α = strong across-site rate heterogeneity. Observed 0.85–18.4 |

Under HKY, IQ-TREE's "Rate parameter R" block fixes all transversions at 1.0 and
reports the shared transition rate on the `A-G` (== `C-T`) line — that value is κ.

**Comparability.** Both are dimensionless — a rate ratio and a shape parameter — so
neither depends on alignment length or taxon count by construction.

### κ vs the model-free `titv`, measured on 187 loci

They agree at **ρ = 0.834**: clearly the same signal, but not interchangeable.

| | κ | `titv` |
|---|---|---|
| range (median) | 2.16–4.89 (3.76) | 0.78–1.94 (1.25) |
| ρ with mean p-distance | **−0.55** | **−0.76** |

- **Scale.** κ/`titv` averages 3.00 (sd 0.30, range 2.46–4.20), so `titv` understates
  the rate ratio ~threefold and by a *varying* amount — `lm(κ ~ titv)` gives only
  R² = 0.70, so it is not a clean rescaling.
- **Saturation.** The κ/`titv` gap widens monotonically with divergence — 2.81, 2.94,
  3.04, 3.21 across divergence quartiles (ρ = +0.52 with mean p-distance). That is
  multiple hits erasing observed transitions: exactly the bias a likelihood model
  corrects for and a raw distance ratio cannot. It is also why `titv` is the more
  divergence-confounded of the two (−0.76 vs −0.55).

**Recommendation: use `kappa`.** `titv` is retained as a cheap model-free
cross-check, useful where the model fit is doubtful.

### `gamma_alpha` vs `qirp_own.rate_var`

Both measure across-site rate heterogeneity, from different estimators. They agree in
the expected direction — ρ(`gamma_alpha`, `rate_var`) = **−0.441**, negative because
low α means high heterogeneity — but only moderately. That independently corroborates
the caveat under `qirp_own`: parsimony site rates compress the rate distribution
relative to ML. If a single rate-heterogeneity covariate is wanted, `gamma_alpha` is
the better-founded one, and it is free.

## `treedist_species` — gene tree vs species tree

**Traditional RF is deliberately excluded**, on two grounds. It gives no partial
credit — one rogue taxon destroys every split it participates in, so RF reports
near-maximal distance for otherwise-identical trees, which over up to 557 taxa is
the normal case. And it takes only even integers up to 2(n−3), so both its
resolution and its range depend on taxon count.

| Column | Meaning |
|---|---|
| `cid` | **Clustering Information Distance** — treats each split as a *clustering* of tips and measures the mutual clustering information between the two trees. Splits need not correspond exactly, so near-matches get partial credit. A distance: **0 = identical** |
| `discord_support_weighted` | share of the gene tree's total branch support sitting on splits **absent** from the species tree |
| `ntips_shared` | |

`normalize = TRUE` scales by the information content of the trees compared, so 0 =
identical and 1 = as different as trees of this size can be — which is what makes
the value valid across the 8–557 tip range. Smith (2020) *Bioinformatics*
36:5007–5013.

**Why CID and not MSID or Nye.** Nye similarity is dropped because 1, not 0, meant
identical trees. Between the two distances the empirical choice barely matters — on
187 loci they correlate at Spearman **+0.987** and both hit exactly 0 for the same 5
loci — but CID uses more of the available range here (0–0.305 vs MSID's 0–0.188),
giving ~1.6× the dynamic range for downstream modelling. Both correlate with taxon
count almost identically (+0.57 vs +0.58).

**No branch-length weighted distance is reported.** The gene trees are in
substitutions/site while the ASTRAL species tree is in coalescent units, so any
length-weighted comparison between them mixes incommensurable quantities. The
earlier `wrf_norm` column was removed for this reason.

`discord_support_weighted` weights each gene-tree split by its own aLRT support
rather than counting it 1/0, so a weakly supported conflict is penalised less than a
well-supported one — the dominant noise source in single-locus gene trees. Returns
NA rather than guessing whenever tip-dropping would break the support-to-node
mapping.

## PhyloMAd statistics — excluded, except one

`phylomad_adequacy` has been **removed**. Only `unc_lnl_per_site` survives, folded
into `site_entropy` (documented there).

**The conceptual problem.** PhyloMAd's statistics are *test statistics for a
parametric bootstrap*. Their value comes from comparison against a null distribution
simulated for that one alignment, which yields an adequacy verdict. Producing that
null costs `Nsims` simulations plus an IQ-TREE run per locus — prohibitive at
112,199 loci. Stripped of the null they are bare descriptive statistics.

**The empirical problem.** Correlations against what we already compute, on 187 loci
spanning all clades (Spearman):

| phylomad statistic | duplicates | ρ | kept? |
|---|---|---|---|
| `df` / tree length | `stemminess.stemminess_per_branch` | **−1.00** — the same quantity, sign-flipped | no |
| χ² → Cramér's V | `aitchison_basefreq.ait_totvar` | **+0.99** | no |
| entropy test → Cohen's d | `site_entropy.entropy_mean_mm` | −0.60 | no |
| biodiversity → rarefied | `site_entropy.entropy_mean_mm` | +0.58 | no |
| unconstrained lnL → per-site | `site_entropy.entropy_mean_mm` | **+0.22** | **yes** |

**Also not included:** `get.phylosignal.metrics.R` (quartet entropy, internode
certainty, d- and kc-statistics). It operates on quartet CF vectors, which `scf_own`
now produces per branch, so it could be added — but it would inherit the same
redundancy risk against `scf_own`'s own columns.

The rescaling arguments for the dropped statistics are preserved in git history if
any need revisiting.

---

## Running

```bash
scripts/run_metrics.sh --list                      # what's registered
scripts/run_metrics.sh --all --jobs 5 --cores 6    # 5 metrics at once, 6 cores each
scripts/run_metrics.sh scf_own -- --clade AnuraAnc1 --limit 50   # scope down
Rscript scripts/merge_metrics.R                    # join into results/locus_metrics.tsv
```

Runs are resumable (per-locus rows, block-flushed) and per-locus seeds are
deterministic, so quartet sampling does not depend on sharding. A locus that fails
gets an all-NA row plus its error message rather than killing the run.

### Checking output before you trust it

```bash
# per-metric row and error counts
for f in results/metrics/*.tsv; do
  printf '%-24s rows=%-7s errors=%s\n' "$(basename $f .tsv)" \
    "$(tail -n +2 $f | wc -l)" "$(awk -F'\t' 'NR>1 && $NF!="NA"' $f | wc -l)"
done

# what went wrong, grouped
awk -F'\t' 'FNR>1 && $NF!="NA"{print FILENAME": "$NF}' results/metrics/*.tsv |
  sort | uniq -c | sort -rn | head
```

Every metric should report 112,199 rows. Fewer means the run was interrupted — just
re-run, resume will fill the gap. A nonzero error count is not necessarily a bug:
small loci legitimately fail the minimum-taxa guards (`delta_plot` and `qirp_own`
need ≥4 tips, `treedist_species` needs ≥4 taxa shared with the species tree). Read
the messages before assuming otherwise.

**Status.** All 11 metrics run clean on a 187-locus sample spanning every clade
(0 errors), including trees up to 557 tips. A full run over all 112,199 loci is
**in progress** — see `logs/run_full.log`. Update this line with the outcome.

**Measured budget** for all 112,199 loci, extrapolated from that sample:

| metric | @12 cores | @32 cores |
|---|---|---|
| `scf_own` + `subst_model` (shared IQ-TREE run) | ~4.3h | ~1.6h |
| `qirp_own` (capped) | ~1.2h | ~0.5h |
| `treedist_species` | ~0.7h | ~0.3h |
| all others combined | ~1.2h | ~0.5h |

**Recommended schedule**, wrapped up in `scripts/run_full.sh`. `scf_own` is ~60% of
the total, so an even `--all` split starves it; run it alone with everything first,
then the rest:

```bash
scripts/run_detached.sh 12    # both phases, merge, summary; survives disconnects
scripts/run_status.sh         # check progress at any time
```

Phase 2 re-runs `scf_own` as a no-op (resume skips completed loci) and `subst_model`
reads phase 1's cache, so neither repeats work. The script deliberately does **not**
use `set -e`: one failing metric must not abort the other ten.

`scf_own` dominates; switching `IQTREE_MODEL` to `JC` would cut it to ~2.0h at 12
cores at ρ=0.99 fidelity, if that trade becomes worthwhile.
