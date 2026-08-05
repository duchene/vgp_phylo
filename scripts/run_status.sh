#!/usr/bin/env bash
# Progress report for a detached run. Safe to call any time, by anyone.
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
total=$(( $(wc -l < "$ROOT/data/manifest.tsv") - 1 ))

pid=$(cat "$ROOT/logs/run_full.pid" 2>/dev/null || true)
if [ -n "${pid:-}" ] && kill -0 "$pid" 2>/dev/null; then
  echo "RUNNING (PID $pid, started $(ps -o lstart= -p "$pid" 2>/dev/null | xargs))"
else
  echo "NOT RUNNING$([ -n "${pid:-}" ] && echo " (stale PID file)")"
fi
echo
printf '%-24s %8s %8s %8s\n' metric rows pct errors
for f in "$ROOT"/results/metrics/*.tsv; do
  [ -e "$f" ] || continue
  n=$(tail -n +2 "$f" | wc -l)
  e=$(awk -F'\t' 'NR>1 && $NF!="NA"' "$f" | wc -l)
  printf '%-24s %8s %7.1f%% %8s\n' "$(basename "$f" .tsv)" "$n" "$(echo "scale=1;100*$n/$total"|bc)" "$e"
done
echo
echo "phase log:"; grep -E 'PHASE|MERGE|COMPLETE|START' "$ROOT/logs/run_full.log" 2>/dev/null | tail -4
echo "latest block:"; tail -1 "$ROOT/logs/scf_own.log" 2>/dev/null
echo; echo "load: $(cut -d' ' -f1-3 /proc/loadavg)  (32 physical cores)"
