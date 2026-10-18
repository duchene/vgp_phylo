#!/bin/bash
## One Slurm array task of the ModelFinder run: the loci of shard
## $SLURM_ARRAY_TASK_ID (round-robin over the job table), run CORES at a time,
## one single-threaded IQ-TREE per locus. Submitted by scripts/submit_modelfinder.sh.
##
## Each locus first goes through scripts/dedup_loci.R, which collapses identical and
## near-identical sequences (p <= MF_NEARIDENT, default 0.001) and prunes the gene tree
## to match; it streams the reduced loci straight into the IQ-TREE workers.
##
## Model selection only, on the locus's FIXED gene tree (-te, no tree search):
##   base models  JC F81 K2P HKY TN GTR   (-mfreq F: F81/HKY/TN/GTR get +F; JC/K2P
##                                         stay equal-frequency -> no SYM/TNe/... twins)
##   rates        E (= 1 class), +R2..+R15, +H2..+H15 (GHOST, linked parameters)
##   -> up to 174 candidate models; BIC picks the best. Branch lengths are re-optimised per
##      model on the fixed topology (+H needs one set per class). ModelFinder stops
##      raising k within a family once BIC stops improving, so high k is only
##      fitted where it pays (on the panel +H never got past H3 and +R past R5 at k <= 5).
## -keep-ident is required: without it ModelFinder drops identical sequences, which
## then no longer match the fixed gene tree ("Tree taxon ... does not appear").
##
## Per locus, kept in <outdir>/<clade>/ (gzipped): .iqtree (best-model parameters),
## .model.gz (lnL/df/tree length of every model tested), .treefile (best-model
## branch lengths; per-class lengths in [..] for +H), .ckp.gz (full checkpoint).
## A locus counts as done once its .iqtree.gz exists, so resubmitting resumes.
## One timing row per locus (wall/CPU seconds, max RSS, whether the -safe retry was
## needed) is appended to <outdir>/timing/shard_NNN.tsv.
##
## Usage: modelfinder_shard.sh <jobs.tsv> <nshards> <cores> <outdir>
set -euo pipefail
JOBS="$1"; NSHARDS="$2"; CORES="$3"; OUTDIR="$4"
I="${SLURM_ARRAY_TASK_ID:?must run as an array task}"
export IQTREE="${IQTREE:-/faststorage/project/DeiC-KU-L04/shared_software/iqtree3/iqtree-3.1.4-Linux/bin/iqtree3}"
export OUTDIR
export SCRATCH="${TMPDIR:-/tmp}/mf_${SLURM_JOB_ID:-$$}_$I"
export TIMING="$OUTDIR/timing/shard_$(printf '%03d' "$I").tsv"
mkdir -p "$SCRATCH" "$OUTDIR/timing"
trap 'rm -rf "$SCRATCH"' EXIT

run_one() {   # locus_id clade ntips_kept nsites fasta treefile seed
  local id="$1" clade="$2" ntips="$3" nsites="$4" fa="$5" tre="$6" seed="$7"
  local dest="$OUTDIR/$clade" pre="$SCRATCH/$id"
  [ -s "$dest/$id.iqtree.gz" ] && return 0
  mkdir -p "$dest"
  local t0 t1 st=0 safe=0
  mf() {   # extra IQ-TREE args, e.g. -safe
    /usr/bin/time -f "%U %S %M" -o "$pre.rus" \
      "$IQTREE" -s "$fa" -te "$tre" --seqtype DNA \
        -m MF -mset JC,F81,K2P,HKY,TN,GTR -mfreq F -mrate E,R,H -cmin 2 -cmax 15 \
        -keep-ident -T 1 -seed "$seed" --prefix "$pre" -quiet -redo "$@" > "$pre.out" 2>&1
  }
  t0=$(date +%s.%N)
  mf || st=$?
  # Many-class +H fits on the largest trees can underflow in the fast kernel; IQ-TREE's
  # own advice is to rerun with the safe (slower, scaled) likelihood kernel.
  if [ "$st" -ne 0 ] && grep -q "Numerical underflow" "$pre.log" 2>/dev/null; then
    safe=1; st=0; mf -safe || st=$?
  fi
  t1=$(date +%s.%N)
  local best=NA
  [ -s "$pre.iqtree" ] && best=$(grep -m1 "^Best-fit model according to BIC:" "$pre.iqtree" | awk '{print $NF}')
  read -r usr sys rss < <(tail -n 1 "$pre.rus") || true
  if [ "$st" -eq 0 ] && [ -s "$pre.iqtree" ]; then
    gzip -f "$pre.iqtree" "$pre.treefile"
    mv "$pre.model.gz" "$pre.ckp.gz" "$pre.treefile.gz" "$dest/"
    mv "$pre.iqtree.gz" "$dest/"   # last: its presence marks the locus done
    # Explicit names, never a glob: a glob makes bash read the whole clade directory
    # (BirdsAnc0 holds ~190k files; one listing takes ~5 s and loads the metadata server).
    rm -f "$dest/$id.FAILED.log.gz" "$dest/$id.FAILED.out"
  else
    gzip -c "$pre.log" > "$dest/$id.FAILED.log.gz" 2>/dev/null || cp "$pre.out" "$dest/$id.FAILED.out"
  fi
  # One short line, appended in a single write -> atomic across the parallel workers.
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%.2f\t%s\t%s\t%s\t%s\n' "$id" "$clade" "$ntips" "$nsites" \
    "$st" "$safe" "$(awk -v a="$t0" -v b="$t1" 'BEGIN{print b-a}')" "${usr:-NA}" "${sys:-NA}" "${rss:-NA}" "$best" >> "$TIMING"
  rm -f "$pre".*
}
export -f run_one

echo "[$(date -Is)] modelfinder shard $I/$NSHARDS on $(hostname), $CORES workers -> $OUTDIR"
[ -s "$TIMING" ] || printf 'locus_id\tclade\tntips\tnsites\texit\tsafe\twall_s\tuser_s\tsys_s\tmaxrss_kb\tbest_bic\n' > "$TIMING"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd -P)"
[ -f "$ROOT/pixi.toml" ] || ROOT="$SLURM_SUBMIT_DIR"   # sbatch runs a spooled copy
export PATH=/faststorage/project/DeiC-KU-L04/shared_software/pixi/current/bin:$PATH
eval "$(pixi shell-hook --frozen --manifest-path "$ROOT/pixi.toml")"
Rscript "$ROOT/scripts/dedup_loci.R" "$JOBS" "$I" "$NSHARDS" "$SCRATCH" "$OUTDIR" |
  xargs -P "$CORES" -L 1 bash -c 'run_one "$@"' _
echo "[$(date -Is)] done"
