#!/bin/bash
## Submit metrics to Slurm on GenomeDK: one job array per metric (one task per
## shard), then a dependent gather job that checks and concatenates the shards.
## Replaces the workstation runners (run_full.sh / run_detached.sh).
##
## scf_own populates the IQ-TREE cache that subst_model reads, so when both are
## submitted together, subst_model's array waits for scf_own's gather. Every other
## metric is independent and starts at once.
##
## Usage (from anywhere):
##   scripts/submit_metrics.sh --all                         # full production run
##   scripts/submit_metrics.sh --all --loci data/panel_187.tsv --outdir results/panel_187
##   scripts/submit_metrics.sh scf_own titv --shards 10 --cores 8
##
## Options:
##   --all             every registered metric
##   --loci FILE       restrict to these locus_ids (one shard by default)
##   --outdir DIR      default results/metrics (full run) or results/<loci basename>
##   --shards N        default: ceil(nloci / 2000), i.e. one block per task
##   --cores N         cores per task (default 4; IQ-TREE gains little beyond ~8)
##   --mem M           memory per task (default 8G)
##   --time T          wall time per task (default 04:00:00)
##   --merge           after all gathers, run merge_metrics.R (full run only)
##   --dry-run         print the sbatch commands, submit nothing
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
export PATH=/faststorage/project/DeiC-KU-L04/shared_software/pixi/current/bin:$PATH
ACCOUNT=DeiC-KU-L04

METRICS=(); LOCI=""; OUTDIR=""; SHARDS=""; CORES=4; MEM=8G; TIME=04:00:00
MERGE=0; DRY=0; ALL=0
while [ $# -gt 0 ]; do
  case "$1" in
    --all)     ALL=1 ;;
    --loci)    LOCI="$(realpath "$2")"; shift ;;
    --outdir)  OUTDIR="$2"; shift ;;
    --shards)  SHARDS="$2"; shift ;;
    --cores)   CORES="$2"; shift ;;
    --mem)     MEM="$2"; shift ;;
    --time)    TIME="$2"; shift ;;
    --merge)   MERGE=1 ;;
    --dry-run) DRY=1 ;;
    -*)        echo "unknown option $1" >&2; exit 2 ;;
    *)         METRICS+=("$1") ;;
  esac
  shift
done

cd "$ROOT"
[ -f data/manifest.tsv ] || { echo "data/manifest.tsv missing: run scripts/build_manifest.sh 32" >&2; exit 1; }
if [ "$ALL" -eq 1 ]; then
  mapfile -t METRICS < <(pixi run --frozen Rscript -e \
    "source('config/project.R'); source('scripts/lib/load.R'); cat(sort(available_metrics()), sep='\n')")
fi
[ ${#METRICS[@]} -gt 0 ] || { echo "no metrics given (name them, or --all)" >&2; exit 2; }

if [ -n "$LOCI" ]; then
  NLOCI=$(grep -cv '^locus_id' "$LOCI")
  OUTDIR="${OUTDIR:-results/$(basename "$LOCI" .tsv)}"
  SHARDS="${SHARDS:-1}"
else
  NLOCI=$(( $(wc -l < data/manifest.tsv) - 1 ))
  OUTDIR="${OUTDIR:-results/metrics}"
  SHARDS="${SHARDS:-$(( (NLOCI + 1999) / 2000 ))}"
fi
OUTDIR="$(realpath -m "$OUTDIR")"
[ "$SHARDS" -le "$NLOCI" ] || SHARDS="$NLOCI"
mkdir -p logs "$OUTDIR"

sub() {   # sbatch wrapper; prints the job id
  if [ "$DRY" -eq 1 ]; then echo "sbatch $*" >&2; echo "DRY$RANDOM"; else sbatch --parsable "$@"; fi
}

echo "$NLOCI loci, $SHARDS shard(s) x $CORES cores, -> $OUTDIR"
declare -A GATHER
order=()
for m in "${METRICS[@]}"; do [ "$m" == scf_own ] && order=(scf_own "${order[@]}") || order+=("$m"); done

for m in "${order[@]}"; do
  dep=()
  if [ "$m" == subst_model ] && [ -n "${GATHER[scf_own]:-}" ]; then
    dep=(--dependency="afterok:${GATHER[scf_own]}")
  fi
  arr=$(sub --account="$ACCOUNT" --job-name="m_$m" --array="1-$SHARDS" \
          --cpus-per-task="$CORES" --mem="$MEM" --time="$TIME" \
          --output="logs/%x_%A_%a.out" "${dep[@]}" \
          scripts/slurm/metric_shard.sh "$m" "$SHARDS" "$CORES" "$OUTDIR" "$LOCI")
  gat=$(sub --account="$ACCOUNT" --job-name="g_$m" --dependency="afterok:$arr" \
          --cpus-per-task=1 --mem=2G --time=00:20:00 --output="logs/%x_%j.out" \
          scripts/slurm/gather_metric.sh "$m" "$SHARDS" "$OUTDIR" "$NLOCI")
  GATHER[$m]=$gat
  printf '  %-22s array %-10s gather %s\n' "$m" "$arr" "$gat"
done

if [ "$MERGE" -eq 1 ]; then
  [ -z "$LOCI" ] || { echo "--merge only applies to the full run" >&2; exit 2; }
  deps=$(IFS=:; echo "${GATHER[*]}")
  mrg=$(sub --account="$ACCOUNT" --job-name=merge_metrics --dependency="afterok:$deps" \
          --cpus-per-task=1 --mem=16G --time=01:00:00 --output="logs/%x_%j.out" \
          --wrap="cd '$ROOT' && pixi run --frozen Rscript scripts/merge_metrics.R")
  echo "  merge                  $mrg"
fi
echo "monitor: squeue --me   |  logs/m_<metric>_<array>_<task>.out"
