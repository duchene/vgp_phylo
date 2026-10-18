#!/bin/bash
## Progress of a ModelFinder run (scripts/submit_modelfinder.sh). Read-only; fine to
## run on the login node.
##
## Usage: scripts/modelfinder_status.sh [--full] [outdir]   (default results/modelfinder)
##
## Reports Slurm array state, loci finished / failed, -safe retries, CPU-hours spent,
## and a rough ETA. Counts come from the per-shard timing tables (last row per locus
## wins), never from listing the ~450k output files.
##
## Speed: the shared filesystem slows down when the cluster is busy (even a plain ls
## can take ~15 s), and a file being appended to from a compute node costs ~0.5 s to
## open however many are read in parallel. So an ACTIVE run is sampled by default:
## every 10th shard's timing file is read and counts are scaled up. Shards are
## round-robin over the manifest, so each is a representative sample. --full reads
## all of them (exact; slow while running, fast once the run is over).
set -uo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
FULL=0; OUT=""
for a in "$@"; do case "$a" in --full) FULL=1 ;; *) OUT="$a" ;; esac; done
OUT="$(realpath -m "${OUT:-$ROOT/results/modelfinder}")"
[ -s "$OUT/jobs.tsv" ] || { echo "no run in $OUT (jobs.tsv missing)"; exit 1; }
TOTAL=$(( $(wc -l < "$OUT/jobs.tsv") - 1 ))
echo "ModelFinder run: $OUT ($TOTAL loci)"

# ---- Slurm: submit_modelfinder.sh appends each array job id to <outdir>/slurm_jobs
ACTIVE=0
JOB=$(tail -n 1 "$OUT/slurm_jobs" 2>/dev/null)
if [ -n "$JOB" ] && command -v squeue >/dev/null; then
  R=$(squeue -h -j "$JOB" -t R 2>/dev/null | wc -l)
  P=$(squeue -h -r -j "$JOB" -t PD 2>/dev/null | wc -l)
  C=$(sacct -n -X -j "$JOB" -s CD --format=JobID 2>/dev/null | wc -l)
  F=$(sacct -n -X -j "$JOB" -s F,TO,OOM,CA,NF --format=JobID 2>/dev/null | wc -l)
  [ $((R + P)) -gt 0 ] && ACTIVE=1
  echo "  Slurm array $JOB: $R running, $P pending, $C completed, $F failed/timeout/cancelled"
  START=$(sacct -n -X -j "$JOB" --format=Start 2>/dev/null | grep -v Unknown | sort | head -1)
fi

# ---- loci, from the timing tables
mapfile -t FILES < <(ls "$OUT"/timing/shard_*.tsv 2>/dev/null)
NF_ALL=${#FILES[@]}
[ "$NF_ALL" -gt 0 ] || { echo "  no timing files yet"; exit 0; }
SCALE=1; NOTE="exact"
if [ "$ACTIVE" -eq 1 ] && [ "$FULL" -eq 0 ] && [ "$NF_ALL" -ge 20 ]; then
  mapfile -t FILES < <(printf '%s\n' "${FILES[@]}" | awk 'NR % 10 == 1')
  SCALE=$(awk -v a="$NF_ALL" -v b="${#FILES[@]}" 'BEGIN { print a / b }')
  NOTE="estimated from ${#FILES[@]} of $NF_ALL shards; --full for exact"
fi
read -r DONE FAILED SAFE CPU_H <<<"$(cat "${FILES[@]}" |
  awk -F'\t' -v k="$SCALE" '$1 != "locus_id" { row[$1] = $0 }
    END { for (id in row) { split(row[id], f, "\t")
            if (f[5] == 0) { d++; c += f[8] + f[9]; s += f[6] } else x++ }
          printf "%d %d %d %.1f\n", d * k, x * k, s * k, c * k / 3600 }')"

pct() { awk -v a="$1" -v b="$TOTAL" 'BEGIN { printf "%.1f", 100 * a / b }'; }
echo "  loci done     $DONE ($(pct "$DONE")%)   [$NOTE]"
echo "  loci failed   $FAILED   (logs: <clade>/<locus>.FAILED.*; resubmitting reruns them)"
echo "  -safe retries $SAFE"
echo "  CPU-hours     $CPU_H   (IQ-TREE user+sys of finished loci)"

# Rough ETA from the completion rate since the first task started. Ignores queue
# changes and the slower large loci, so treat it as a guide only.
if [ "$ACTIVE" -eq 1 ] && [ -n "${START:-}" ] && [ "$DONE" -gt 0 ]; then
  EL=$(( $(date +%s) - $(date -d "$START" +%s) ))
  awk -v d="$DONE" -v t="$TOTAL" -v e="$EL" 'BEGIN {
    printf "  elapsed %.1f h; ~%.1f h to go at the current rate\n", e / 3600, (t - d) / (d / e) / 3600 }'
fi
