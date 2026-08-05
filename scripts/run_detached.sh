#!/usr/bin/env bash
# Launch the full run fully detached, so it survives SSH disconnects, terminal
# closes, and this session ending.
#
#   setsid    puts it in a NEW session and process group, so anything that kills
#             the launching shell's process group cannot reach it
#   nohup     ignores SIGHUP
#   </dev/null it never blocks waiting on a terminal that is gone
#
# Usage: scripts/run_detached.sh [cores]
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CORES="${1:-12}"
mkdir -p "$ROOT/logs"
setsid nohup "$ROOT/scripts/run_full.sh" "$CORES" \
  > "$ROOT/logs/run_full.stdout" 2>&1 < /dev/null &
disown || true
sleep 2
if [ -f "$ROOT/logs/run_full.pid" ]; then
  echo "launched detached, PID $(cat "$ROOT/logs/run_full.pid"), ${CORES} cores"
  echo "  progress: scripts/run_status.sh"
  echo "  stop:     kill \$(cat logs/run_full.pid)"
else
  echo "WARNING: no PID file after 2s -- check logs/run_full.stdout" >&2
  tail -5 "$ROOT/logs/run_full.stdout" >&2 || true
fi
