#!/bin/bash
## One Slurm array task: shard $SLURM_ARRAY_TASK_ID of NSHARDS for one metric.
## Submitted by scripts/submit_metrics.sh, which sets the #SBATCH resources.
##
## Each task writes its OWN file, <outdir>/shards/<metric>/shard_<i>.tsv, so no two
## tasks ever append to the same TSV (the reason the old workstation runner needed
## an flock). Shards are round-robin over the manifest order (run_metric.R --shard),
## so a resubmitted task resumes its own shard file at block granularity.
##
## Usage: metric_shard.sh <metric> <nshards> <cores> <outdir> [loci_file]
set -euo pipefail
METRIC="$1"; NSHARDS="$2"; CORES="$3"; OUTDIR="$4"; LOCI="${5:-}"
I="${SLURM_ARRAY_TASK_ID:?must run as an array task}"

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd -P)"
# sbatch copies the script to a spool dir; fall back to the submit dir.
[ -f "$ROOT/pixi.toml" ] || ROOT="$SLURM_SUBMIT_DIR"
export VGP_PHYLO_ROOT="$ROOT"
cd "$ROOT"

export PATH=/faststorage/project/DeiC-KU-L04/shared_software/pixi/current/bin:$PATH
eval "$(pixi shell-hook --frozen --manifest-path "$ROOT/pixi.toml")"

OUT="$OUTDIR/shards/$METRIC/shard_$(printf '%03d' "$I").tsv"
mkdir -p "$(dirname "$OUT")"

ARGS=(--metric "$METRIC" --cores "$CORES" --shard "$I/$NSHARDS" --out "$OUT")
[ -n "$LOCI" ] && ARGS+=(--loci "$LOCI")

echo "[$(date -Is)] $METRIC shard $I/$NSHARDS on $(hostname), $CORES cores -> $OUT"
Rscript "$ROOT/scripts/run_metric.R" "${ARGS[@]}"
echo "[$(date -Is)] done"
