#!/usr/bin/env bash
# Run several metrics concurrently, one process per metric.
#
#   scripts/run_metrics.sh --list                       # what's available
#   scripts/run_metrics.sh scf_own tree_shape           # named metrics
#   scripts/run_metrics.sh --all                        # everything on METRIC_PATH
#   scripts/run_metrics.sh --all -- --clade AnuraAnc1   # args after -- go to run_metric.R
#
# Total load is roughly (--jobs x --cores). This box is shared, so both default low.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export VGP_PHYLO_ROOT="$ROOT"

JOBS=""
CORES=""
METRICS=()
PASSTHRU=()

while [ $# -gt 0 ]; do
  case "$1" in
    --list)  Rscript -e "source('$ROOT/config/project.R'); source('$ROOT/scripts/lib/load.R');
                for (n in sort(available_metrics())) {
                  m <- tryCatch(load_metric(n), error = function(e) NULL)
                  cat(sprintf('%-26s %s\n', n, if (is.null(m)) '<<failed to load>>' else m\$description))
                }"; exit 0 ;;
    --all)   mapfile -t METRICS < <(Rscript -e "source('$ROOT/config/project.R'); source('$ROOT/scripts/lib/load.R'); cat(sort(available_metrics()), sep='\n')"); shift ;;
    --jobs)  JOBS="$2"; shift 2 ;;
    --cores) CORES="$2"; shift 2 ;;
    --)      shift; PASSTHRU=("$@"); break ;;
    -*)      echo "unknown option: $1" >&2; exit 2 ;;
    *)       METRICS+=("$1"); shift ;;
  esac
done

[ ${#METRICS[@]} -gt 0 ] || { echo "no metrics given (try --list or --all)" >&2; exit 2; }

# Defaults come from config/project.R so there is one place to tune them.
# The `|| true` matters: read returns non-zero at EOF-without-newline, which under
# `set -e` would abort the script even though both variables were assigned.
read -r DEF_JOBS DEF_CORES < <(Rscript -e "source('$ROOT/config/project.R'); cat(DEFAULTS\$metric_jobs, DEFAULTS\$cores, '\n')") || true
JOBS="${JOBS:-$DEF_JOBS}"
CORES="${CORES:-$DEF_CORES}"

mkdir -p "$ROOT/logs"
echo "running ${#METRICS[@]} metric(s), $JOBS at a time, $CORES cores each" >&2

# One log per metric. xargs -P gives the across-metric parallelism; each child
# gets its own log so concurrent progress lines never interleave.
printf '%s\n' "${METRICS[@]}" | xargs -I{} -P "$JOBS" bash -c '
  set -o pipefail
  log="'"$ROOT"'/logs/{}.log"
  echo "=== $(date -Is) start {} ===" >> "$log"
  if Rscript "'"$ROOT"'/scripts/run_metric.R" --metric {} --cores '"$CORES"' '"${PASSTHRU[*]:-}"' >> "$log" 2>&1; then
    echo "OK   {}" >&2
  else
    echo "FAIL {}  (see logs/{}.log)" >&2
  fi
  echo "=== $(date -Is) end {} ===" >> "$log"
'

echo "logs in $ROOT/logs/" >&2
