#!/usr/bin/env bash
# Build data/manifest.tsv -- the single source of truth for which loci exist and
# how big they are. Every metric driver reads this instead of re-scanning 112k files.
#
# Usage: scripts/build_manifest.sh [n_parallel]
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
JOBS="${1:-16}"
OUT="$ROOT/data/manifest.tsv"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/manifest.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

echo "scanning $ROOT/data/loci with $JOBS parallel workers ..." >&2
find "$ROOT/data/loci" -name '*.fa' | sort > "$WORK/filelist"
split -n "l/$JOBS" "$WORK/filelist" "$WORK/chunk."

# One serial xargs stream per chunk, each appending to its OWN file. Never let
# concurrent awk processes share a stdout -- their buffer flushes interleave
# mid-line and silently corrupt the output.
for chunk in "$WORK"/chunk.*; do
  ( xargs -a "$chunk" -n 400 awk -f "$ROOT/scripts/manifest.awk" > "$chunk.out" ) &
done
wait

{
  printf 'locus_id\tclade\tntips\tnspecies\tncopies_max\tnsites\n'
  cat "$WORK"/chunk.*.out | sort -k2,2 -k1,1
} > "$OUT"

n=$(( $(wc -l < "$OUT") - 1 ))
echo "wrote $OUT ($n loci)" >&2
[ "$n" -eq "$(wc -l < "$WORK/filelist")" ] || { echo "ERROR: row count != fasta count" >&2; exit 1; }

# Cross-check: every locus must have a matching gene tree.
missing=0
while IFS=$'\t' read -r locus clade _; do
  [ -f "$ROOT/data/gene_trees/$clade/$locus.fa.treefile" ] || { echo "MISSING TREE: $clade/$locus" >&2; missing=$((missing+1)); }
done < <(tail -n +2 "$OUT")
echo "loci without a gene tree: $missing" >&2
[ "$missing" -eq 0 ] || exit 1
