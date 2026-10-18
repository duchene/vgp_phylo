#!/bin/bash
## Submit the per-locus ModelFinder run (model selection only, fixed gene tree) as
## one Slurm job array. What is fitted is documented in scripts/slurm/modelfinder_shard.sh.
##
## Usage (from anywhere):
##   scripts/submit_modelfinder.sh --loci data/panel_187.tsv --shards 12   # timing pilot
##   scripts/submit_modelfinder.sh --shards 500 --cores 8 --time 12:00:00 # all loci
##
## Options:
##   --loci FILE   restrict to these locus_ids (default: all loci in the manifest)
##   --outdir DIR  default results/modelfinder (all loci) or results/modelfinder_<loci basename>
##   --shards N    array tasks (default 1)
##   --cores N     concurrent single-threaded IQ-TREE runs per task (default 4)
##   --mem M       memory per task (default 8G)
##   --time T      wall time per task (default 04:00:00)
##   --dry-run     write the job table, print the sbatch command, submit nothing
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
export PATH=/faststorage/project/DeiC-KU-L04/shared_software/pixi/current/bin:$PATH

LOCI=""; OUTDIR=""; SHARDS=1; CORES=4; MEM=8G; TIME=04:00:00; DRY=0
while [ $# -gt 0 ]; do
  case "$1" in
    --loci)    LOCI="$(realpath "$2")"; shift ;;
    --outdir)  OUTDIR="$2"; shift ;;
    --shards)  SHARDS="$2"; shift ;;
    --cores)   CORES="$2"; shift ;;
    --mem)     MEM="$2"; shift ;;
    --time)    TIME="$2"; shift ;;
    --dry-run) DRY=1 ;;
    *)         echo "unknown option $1" >&2; exit 2 ;;
  esac
  shift
done

cd "$ROOT"
if [ -n "$LOCI" ]; then OUTDIR="${OUTDIR:-results/modelfinder_$(basename "$LOCI" .tsv)}"
else OUTDIR="${OUTDIR:-results/modelfinder}"; fi
OUTDIR="$(realpath -m "$OUTDIR")"
mkdir -p logs "$OUTDIR"

JOBS="$OUTDIR/jobs.tsv"
pixi run --frozen Rscript scripts/modelfinder_jobs.R "$JOBS" ${LOCI:+"$LOCI"}
NJOBS=$(( $(wc -l < "$JOBS") - 1 ))
[ "$SHARDS" -le "$NJOBS" ] || SHARDS="$NJOBS"

CMD=(sbatch --parsable --account=DeiC-KU-L04 --job-name=modelfinder
     --array="1-$SHARDS" --cpus-per-task="$CORES" --mem="$MEM" --time="$TIME"
     --output="logs/mf_%A_%a.out"
     scripts/slurm/modelfinder_shard.sh "$JOBS" "$SHARDS" "$CORES" "$OUTDIR")
echo "$NJOBS loci, $SHARDS task(s) x $CORES cores, -> $OUTDIR"
if [ "$DRY" -eq 1 ]; then echo "${CMD[*]}"; exit 0; fi
JOBID=$("${CMD[@]}")
echo "$JOBID" >> "$OUTDIR/slurm_jobs"   # read by scripts/modelfinder_status.sh
echo "array job $JOBID"
