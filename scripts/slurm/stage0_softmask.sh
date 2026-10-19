#!/bin/bash
## Stage 0 covariate: soft-masked (repeat) fraction per locus. The alignments keep
## the assemblies' soft-masking, so lower-case a/c/g/t are bases that the
## species' repeat annotation masked. Per locus, pooled over all sequences:
##   softmask = lower-case acgt / all acgt (gaps and N excluded)
## Also writes the mean of the per-sequence fractions, so one heavily masked
## species cannot dominate unnoticed.
##
## Opens all 112k alignments, so it runs on a compute node:
##   sbatch --account=DeiC-KU-L04 --cpus-per-task=8 --mem=4G --time=02:00:00 \
##     --output=logs/stage0_softmask_%j.out scripts/slurm/stage0_softmask.sh
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd -P)"
[ -f "$ROOT/pixi.toml" ] || ROOT="$SLURM_SUBMIT_DIR"
cd "$ROOT"
OUT=results/stage0/locus_softmask.tsv
mkdir -p results/stage0

# One awk pass per clade, CPUS clades at a time; each writes its own part file.
TMP=$(mktemp -d results/stage0/softmask_parts.XXXX)
trap 'rm -rf "$TMP"' EXIT
cat > "$TMP/softmask.awk" <<'AWK'
function endseq() { if (nb > 0) { sf += nl / nb; ns++ } nl = 0; nb = 0 }
function endlocus() {
  endseq()
  if (id != "") printf "%s\t%s\t%d\t%d\t%.6g\t%.6g\n", id, clade, L, B, (B ? L / B : -1), (ns ? sf / ns : -1)
}
FNR == 1 { endlocus(); id = FILENAME; sub(/.*\//, "", id); sub(/\.fa$/, "", id); L = 0; B = 0; sf = 0; ns = 0 }
/^>/ { endseq(); next }
{ l = gsub(/[acgt]/, "&"); u = gsub(/[ACGT]/, "&"); nl += l; nb += l + u; L += l; B += l + u }
END { endlocus() }
AWK
N=0
for clade in $(awk -F'\t' 'NR > 1 { print $2 }' data/manifest.tsv | sort -u); do
  find "data/loci/$clade" -name '*.fa' -print0 |
    xargs -0 -n 5000 awk -v clade="$clade" -f "$TMP/softmask.awk" > "$TMP/$clade.tsv" &
  N=$((N + 1)); [ $((N % ${SLURM_CPUS_PER_TASK:-1})) -eq 0 ] && wait
done
wait

{ printf 'locus_id\tclade\tsoftmasked_bases\tacgt_bases\tsoftmask\tsoftmask_seqmean\n'; cat "$TMP"/*.tsv; } > "$OUT.part"
n=$(($(wc -l < "$OUT.part") - 1)); want=$(($(wc -l < data/manifest.tsv) - 1))
[ "$n" -eq "$want" ] || { echo "softmask: $n loci, manifest has $want" >&2; exit 1; }
mv "$OUT.part" "$OUT"
echo "wrote $OUT ($n loci)"
