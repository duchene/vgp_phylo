#!/usr/bin/env bash
# Full production run: every metric over all 112,199 loci, on a 32-core budget.
#
# Two phases, because scf_own is ~60% of the total work and an even --all split
# would starve it:
#   Phase 1  scf_own alone with all 32 cores. This also populates the IQ-TREE cache.
#   Phase 2  everything else, 8 metrics x 4 cores. scf_own re-runs as a no-op
#            (resume skips completed loci) and subst_model reads phase 1's cache.
#
# Usage: scripts/run_full.sh [total_cores]
set -uo pipefail            # NOT -e: a single failing metric must not abort the run

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export VGP_PHYLO_ROOT="$ROOT"
CORES="${1:-32}"
LOG="$ROOT/logs/run_full.log"
mkdir -p "$ROOT/logs"

say() { echo "[$(date -Is)] $*" | tee -a "$LOG"; }

say "=== FULL RUN START: $(wc -l < "$ROOT/data/manifest.tsv" | awk '{print $1-1}') loci, ${CORES} cores ==="
say "model: $(Rscript -e "source('$ROOT/config/project.R'); cat(IQTREE_MODEL)")"

say "--- PHASE 1: scf_own on ${CORES} cores ---"
Rscript "$ROOT/scripts/run_metric.R" --metric scf_own --cores "$CORES" >> "$ROOT/logs/scf_own.log" 2>&1
rc=$?
say "PHASE 1 done (exit $rc); rows=$(tail -n +2 "$ROOT/results/metrics/scf_own.tsv" 2>/dev/null | wc -l)"

say "--- PHASE 2: all remaining metrics, 8 jobs x $((CORES/8)) cores ---"
"$ROOT/scripts/run_metrics.sh" --all --jobs 8 --cores "$((CORES/8))" >> "$LOG" 2>&1
say "PHASE 2 done"

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
