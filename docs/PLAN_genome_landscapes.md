# Plan: model, recombination and reticulation landscapes across 112k VGP loci

Four stages. Each describes a family of locus-level signals *along the genome*, and the
last asks how those signals evolve across lineages. Stage 0 is a prerequisite that
defines what "along the genome" means for this dataset.

```
Stage 0  genomic frame ──┬──> Stage 1  substitution-model landscape (+R, +H, +T over k)
                         │        └──> Stage 2  recombination from +T
                         ├──> Stage 3  introgression across depths (Holland's delta)
dated species tree ──────┴──────────> Stage 4  evolution of all parameters across lineages
```

Compute runs on GenomeDK (Slurm arrays), on **all 112,199 loci**. It reuses the existing
conventions: one cached IQ-TREE run per locus, blocks of 2000 loci, and the
`register_metric()` framework.

**Topology is always fixed.** No fit in any stage runs a tree search. Every model is
fitted on the locus's existing IQ-TREE gene tree. Where that tree is unavailable or too
slow to use, a BIONJ tree from the same alignment is the fallback. Only model
parameters, and branch lengths where a model requires them, are estimated.

---

## Stage 0 — Genomic frame (prerequisite)

Shuffled synteny in living species rules out "position in species X". There are two
frames that rearrangements don't affect:

1. **Ancestral-reference coordinates.** These are already in every locus name:
   `<Anc>_<Anc>refChr<N>_<start>_<end>`. `refChr<N>` is a contig of the reconstructed
   ancestral genome. `BirdsAnc0` has 1,939 such contigs, up to 5.75 Mb long, with up to
   272 loci each. Within a clade, distance along a contig is a real physical distance
   in the ancestor.
2. **Conserved linkage groups.** Each sequence label carries that species'
   chromosome/scaffold accession. Contigs that land on the same chromosome in most
   species belong to the same macrosynteny unit. Clustering contigs on that
   co-occurrence matrix rebuilds the linkage groups. Labels come from human chromosome
   (10,109 loci contain GRCh38) and from a micro/macrochromosome proxy (loci per
   species-chromosome, as a stand-in for chromosome size).

**Output:** `results/locus_position.tsv` with columns locus, clade, contig, start, end,
rank along contig, linkage group, human chr, micro/macro call, plus covariates:
GC (already computed), soft-masked fraction, and a coding score (third-position
periodicity of variability, gap lengths that are multiples of 3, open reading frames).

**Spatial toolkit**, reused by Stages 1–3:
- autocorrelation of any locus metric against distance along a contig (variogram, or
  Moran's I by distance bin), compared with a within-contig permutation null;
- segmentation into domains (changepoint or HMM along contigs);
- summaries by linkage group and by micro/macro class.

Caveat: the frame is clade-specific. The 22 references are different ancestors, so
along-genome comparisons are made within a clade and only compared across clades
through linkage groups.

---

## Stage 1 — Substitution-model landscape along the genome

**Model grid.** The exchangeability and frequency part is held fixed at GTR+FO, and
only the heterogeneity family and its size *k* vary:

| Family | IQ-TREE | k | Captures |
|---|---|---|---|
| baseline | `+G4`, `+I+G4` | — | reference |
| FreeRate | `+R k` | 2–8 | across-site rate heterogeneity, free shape |
| GHOST heterotachy | `+H k` / `*H k` (linked/unlinked) | 2–6 | across-site × across-lineage rate variation |
| MAST tree mixture | `+T`, k candidate trees via `-te` | 1–4 | topological heterogeneity within the locus |

**Per locus:** lnL, df, AIC/BIC, best *k* within each family, best family overall,
and parameters (FreeRate rates and weights, GHOST class weights and per-class tree
lengths, MAST tree weights).

**Comparability.** This is the main methodological risk, given the project's
comparability rule:
- Information criteria grow with site count, so ΔBIC is reported **per site**.
- The df of `+H` and `+T` grow with tip count, since each class adds 2n−3 branch
  lengths. The BIC penalty therefore differs between an 8-tip and a 557-tip locus.
- Calibration uses a **parametric bootstrap** on a stratified subset (≈500 loci):
  simulate under the best `+R` model and record how often `+H` or `+T` is wrongly
  preferred, as a function of ntips, nsites and tree depth. The decision thresholds
  come from that error surface rather than raw BIC.

**Topology fixed, branch lengths where needed.** `+G`, `+I+G` and `+R` keep the current
`-blfix` convention on the gene tree. `+H` needs per-class branch lengths and `+T` needs
per-tree branch lengths. Both are estimated on **fixed topologies** (`-te`) with no
tree search. This is the only departure from `-blfix`.

**Cost.** About 20 fits per locus across 112k loci, roughly 2M IQ-TREE runs, all of
them on fixed topologies. The plan is to run every locus, spreading the work across
GenomeDK Slurm arrays. The 187-locus panel is used first to time each family and to
check the MAST and GHOST syntax in IQ-TREE 3.x.

**Deliverables:** `results/metrics/model_landscape.tsv`. Figures: best-family share
and best *k* along contigs, autocorrelation range of each parameter, and relationships
with GC, coding score and micro/macro class.

---

## Stage 2 — Recombination signal from +T

A locus that is really a mosaic of genealogies is fitted better by a tree mixture. On
one locus, more than one genealogy requires intralocus recombination; the other
explanations are alignment or paralogy errors.

**Per locus:**
- best k under MAST, and the tree weights;
- whether `+T` beats **both** `+R` and `+H` from Stage 1, so heterotachy doesn't
  masquerade as topology mixing;
- topological distance between the mixed trees, and *which* splits differ, so that
  the signal maps onto species-tree branches. Stage 4 uses this.

**Telling recombination apart from misspecification.** Recombination predicts that
sites belonging to each tree come in **contiguous blocks**. Misspecification predicts
they are scattered. The test uses per-site posterior tree membership (site-wise
mixture probabilities; the exact output option for MAST to be verified) and checks
for runs along the alignment against a permutation null. A subset is cross-checked
with independent tests: Phi (PhiPack), and 3SEQ or GARD.

**Power and false positives.** Two simulation sets matched to real loci (ntips,
nsites, depth):
- without recombination (Stage 1 best models), for the false-positive rate;
- with recombination (coalescent simulations with known breakpoints), for power as a
  function of tip count, breakpoint position and divergence.

**Along the genome:** rate of recombination signal along contigs and by linkage
group. Predictions to test:
- higher on bird microchromosomes;
- positively related to GC, via GC-biased gene conversion;
- position relative to contig ends, as a weak stand-in for distance to telomeres
  (contigs are not chromosomes).

---

## Stage 3 — Introgression/hybridisation across evolutionary depths

The existing `delta_plot` metric averages Holland's δ over 2000 random quartets per
locus, which hides depth. The revision:

- **δ per quartet, tagged with depth:** the age of the quartet's MRCA on the dated
  tree. An interim version can use substitution-scale height on the species tree.
- **Model-corrected distances:** ML distances under the Stage 1 best model, not
  p-distances. Saturation inflates δ at depth, which would create a false trend with
  age.
- **Null per depth bin:** δ is non-zero under ILS alone, and ILS is concentrated at
  short internodes. Expected δ comes from multispecies-coalescent simulations on the
  species tree without gene flow. The signal is the **excess over this null**.
- **Recommended cheap complement**, using the existing gene trees: asymmetry between
  the two minor quartet topologies around each species-tree branch. ILS predicts
  symmetry, introgression predicts asymmetry. It is directional and branch-specific,
  which δ is not.

**Along the genome:** δ excess and asymmetry along contigs and by linkage group,
cross-referenced with Stage 2. Introgressed ancestry is often retained preferentially
in high-recombination regions.

---

## Stage 4 — Evolution of the parameters across lineages

**Prerequisite: a dated species tree.** The current ASTRAL-Pro tree (564 tips) is not
ultrametric. Collaborators are building the dated tree (see Decisions).

**From locus-level to lineage-level traits.** Most Stage 1–3 outputs are per locus.
Each gets mapped onto species-tree branches or tips:

| Trait | Per branch / tip |
|---|---|
| Biochemical composition | per-species GC (and GC3 where coding), ancestral reconstruction of composition shifts |
| Substitution rate | gene-tree branch length ÷ branch duration, summarised across loci |
| Rate heterogeneity / heterotachy | per-branch spread across GHOST classes; per-lineage share of loci favouring `+H` |
| Recombination emergence | per branch, share of loci whose `+T` trees disagree at that branch (from Stage 2) |
| Reticulation | per branch, δ excess and quartet asymmetry (from Stage 3) |

**Comparative models:**
- continuous traits: BM vs OU vs variable-rates models (e.g. BayesTraits variable
  rates, or RRphylo) to locate rate shifts;
- discrete states, e.g. "recombination signal present": Mk with hidden rates (corHMM),
  which gives the **rates of emergence and loss**;
- correlated evolution: phylogenetic GLS. Life-history covariates (generation time,
  body size) are an optional external addition.

**Coverage caveat.** Loci are tied to clade references. Deep branches are informed
mainly by the ~6,000 `VertebratesAnc1–5` loci (≈350–400 tips each); `VertebratesAnc0`
has only 2 loci. Within-clade branches are informed by that clade's set. Lineage
comparisons therefore mix locus sets, and per-branch traits need the locus count kept
alongside them as a weight.

---

## Decisions

Settled (2026-10-16):

1. **Topology always fixed:** gene tree, with BIONJ as the fast fallback. No tree search
   in any stage.
2. **`+T` candidate trees:** (a) the gene tree, (b) the species tree pruned to the locus's
   tips, and (c) NNI alternatives of the gene tree at its least-supported branches. Support
   values are read from the existing gene-tree node labels, so no bootstrap search is
   needed. k = 1–4 draws from this set in that order.
3. **Scope:** all 112,199 loci for every family, spread across GenomeDK Slurm arrays.

Open:

4. **Dated tree source.** Collaborators are building a dated tree. Stage 4, and the
   depth-binning in Stage 3, wait for it. Until then, Stage 3 can use substitution-scale
   height on the species tree as an interim depth.
