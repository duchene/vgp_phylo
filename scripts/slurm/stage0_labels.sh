#!/bin/bash
## Stage 0, step 1: one row per (locus, sequence) with the species assembly and the
## chromosome/scaffold accession that the sequence label carries:
##   >GCA_008658365.1.CM018088.1|copy1  ->  GCA_008658365.1  CM018088.1
## The assembly accession is GC[AF]_<digits>.<version>; the rest up to "|" is the
## sequence accession (it may itself contain dots). ~17M rows, gzipped.
##
## Opens all 112k alignments, so it runs on a compute node:
##   sbatch --account=DeiC-KU-L04 --cpus-per-task=4 --mem=4G --time=02:00:00 \
##     --output=logs/stage0_labels_%j.out scripts/slurm/stage0_labels.sh
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd -P)"
[ -f "$ROOT/pixi.toml" ] || ROOT="$SLURM_SUBMIT_DIR"
cd "$ROOT"
OUT=results/stage0/locus_labels.tsv.gz
mkdir -p results/stage0

# One awk pass per clade, CPUS clades at a time; each writes its own part file.
TMP=$(mktemp -d results/stage0/labels_parts.XXXX)
trap 'rm -rf "$TMP"' EXIT
cat > "$TMP/labels.awk" <<'AWK'
FNR == 1 { id = FILENAME; sub(/.*\//, "", id); sub(/\.fa$/, "", id) }
/^>/ {
  l = substr($0, 2); sub(/[ \t].*/, "", l); sub(/\|.*/, "", l)
  if (match(l, /^GC[AF]_[0-9]+\.[0-9]+\./))
    print id "\t" clade "\t" substr(l, 1, RLENGTH - 1) "\t" substr(l, RLENGTH + 1)
  else
    print id "\t" clade "\t" l "\tNA"
}
AWK
N=0
for clade in $(awk -F'\t' 'NR > 1 { print $2 }' data/manifest.tsv | sort -u); do
  find "data/loci/$clade" -name '*.fa' -print0 |
    xargs -0 awk -v clade="$clade" -f "$TMP/labels.awk" > "$TMP/$clade.tsv" &
  N=$((N + 1)); [ $((N % ${SLURM_CPUS_PER_TASK:-1})) -eq 0 ] && wait
done
wait

{ printf 'locus_id\tclade\tassembly\tseq_accession\n'; cat "$TMP"/*.tsv; } | gzip > "$OUT.part"
mv "$OUT.part" "$OUT"
echo "rows: $(zcat "$OUT" | tail -n +2 | wc -l); unparsed labels: $(zcat "$OUT" | awk -F'\t' '$4 == "NA"' | wc -l)"
