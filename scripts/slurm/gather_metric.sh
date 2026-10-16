#!/bin/bash
## Concatenate one metric's shard files into <outdir>/<metric>.tsv, after checking
## that every shard is present, all headers agree, and no locus is missing or
## duplicated. Runs as a dependent job (afterok on the array), or by hand.
##
## Usage: gather_metric.sh <metric> <nshards> <outdir> <expected_loci>
set -euo pipefail
METRIC="$1"; NSHARDS="$2"; OUTDIR="$3"; EXPECTED="$4"
SDIR="$OUTDIR/shards/$METRIC"
OUT="$OUTDIR/$METRIC.tsv"

shards=()
for i in $(seq 1 "$NSHARDS"); do
  f="$SDIR/shard_$(printf '%03d' "$i").tsv"
  [ -s "$f" ] || { echo "ERROR: missing or empty $f" >&2; exit 1; }
  shards+=("$f")
done

hdr=$(head -1 "${shards[0]}")
for f in "${shards[@]}"; do
  [ "$(head -1 "$f")" == "$hdr" ] || { echo "ERROR: header mismatch in $f" >&2; exit 1; }
done

tmp="$OUT.tmp.$$"
awk 'FNR == 1 && NR != 1 { next } { print }' "${shards[@]}" > "$tmp"

rows=$(( $(wc -l < "$tmp") - 1 ))
uniq=$(tail -n +2 "$tmp" | cut -f1 | sort -u | wc -l)
errs=$(awk -F'\t' 'NR > 1 && $NF != "NA"' "$tmp" | wc -l)
if [ "$rows" -ne "$EXPECTED" ] || [ "$uniq" -ne "$rows" ]; then
  echo "ERROR: $METRIC has $rows rows / $uniq unique loci, expected $EXPECTED" >&2
  rm -f "$tmp"; exit 1
fi
mv "$tmp" "$OUT"   # atomic: readers never see a partial table
printf '%-24s rows=%-7s errors=%s -> %s\n' "$METRIC" "$rows" "$errs" "$OUT"
