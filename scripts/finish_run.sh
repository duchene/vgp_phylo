#!/usr/bin/env bash
# Finish an incomplete run: bring every metric to full coverage, then merge.
#
# Metrics are run SEQUENTIALLY, each with the full core budget -- deliberately, and
# unlike run_full.sh's phase 2, which runs them concurrently with a small share each.
#
# Why: metric costs differ by more than an order of magnitude, so a concurrent split
# ends in a long tail where one metric plods along on 2 cores while the rest of the
# budget sits idle. That is exactly what happened on the first full run -- ten metrics
# finished, then qirp_own needed 6 more hours on 2 cores when 12 were free.
# Sequential-with-all-cores has no tail. It also suits the two cost profiles here:
# IQ-TREE metrics are memory-bandwidth-bound and gain little from extra cores, while
# qirp_own is CPU-bound and scales linearly.
#
# Already-complete metrics cost seconds (resume skips every finished locus), so this
# is safe and cheap to re-run at any time.
#
# Usage: scripts/finish_run.sh [cores]
set -uo pipefail            # NOT -e: one failing metric must not abort the rest

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export VGP_PHYLO_ROOT="$ROOT"
CORES="${1:-12}"
LOG="$ROOT/logs/run_full.log"
mkdir -p "$ROOT/logs"

exec 9>"$ROOT/logs/.run_full.lock"
if ! flock -n 9; then
  echo "ERROR: a run is already in progress (PID $(cat "$ROOT/logs/run_full.pid" 2>/dev/null || echo unknown))." >&2
  exit 1
fi
echo $$ > "$ROOT/logs/run_full.pid"
trap 'rm -f "$ROOT/logs/run_full.pid"' EXIT

say() { echo "[$(date -Is)] $*" | tee -a "$LOG"; }
total=$(( $(wc -l < "$ROOT/data/manifest.tsv") - 1 ))

say "--- FINISH RUN: sequential, ${CORES} cores each, target ${total} loci ---"

mapfile -t METRICS < <(Rscript -e "source('$ROOT/config/project.R'); source('$ROOT/scripts/lib/load.R'); cat(sort(available_metrics()), sep='\n')")

for m in "${METRICS[@]}"; do
  f="$ROOT/results/metrics/$m.tsv"
  have=$([ -f "$f" ] && tail -n +2 "$f" | wc -l || echo 0)
  if [ "$have" -ge "$total" ]; then
    say "  $m: already complete ($have rows), skipping"
    continue
  fi
  say "  $m: $have/$total -- running on ${CORES} cores"
  Rscript "$ROOT/scripts/run_metric.R" --metric "$m" --cores "$CORES" >> "$ROOT/logs/$m.log" 2>&1
  say "  $m: done, now $(tail -n +2 "$f" 2>/dev/null | wc -l) rows"
done

say "--- MERGE ---"
Rscript "$ROOT/scripts/merge_metrics.R" >> "$LOG" 2>&1
say "MERGE done"

say "=== PER-METRIC SUMMARY ==="
for f in "$ROOT"/results/metrics/*.tsv; do
  n=$(tail -n +2 "$f" | wc -l)
  e=$(awk -F'\t' 'NR>1 && $NF!="NA"' "$f" | wc -l)
  printf '%-24s rows=%-7s errors=%s\n' "$(basename "$f" .tsv)" "$n" "$e" | tee -a "$LOG"
done
say "=== FULL RUN COMPLETE ==="
