#!/bin/bash
## Finish a ModelFinder run: parse it, then pack the per-locus files into one
## uncompressed tar per clade (the members are already gzipped), so ~450k small
## files become ~22. Run as a Slurm job after the array, e.g.
##
##   sbatch --account=DeiC-KU-L04 --dependency=afterany:<array> --cpus-per-task=8 \
##     --mem=16G --time=06:00:00 --output=logs/mf_finalize_%j.out \
##     scripts/slurm/modelfinder_finalize.sh results/modelfinder
##
## Stops before packing anything if any locus is missing or unparseable
## (parse_modelfinder.R exits non-zero): resubmit the array first.
## A clade's directory is deleted only after its tar lists exactly the expected
## 4 files per locus (.iqtree.gz .model.gz .treefile.gz .ckp.gz).
##
## Reading one locus afterwards, without unpacking:
##   tar -xOf <outdir>/archives/<clade>.tar <clade>/<locus>.iqtree.gz | zcat
## Re-parsing needs the directories back: tar -xf <outdir>/archives/<clade>.tar -C <outdir>
set -euo pipefail
OUT="$(realpath "${1:?usage: modelfinder_finalize.sh <outdir>}")"
CORES="${SLURM_CPUS_PER_TASK:-1}"

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd -P)"
[ -f "$ROOT/pixi.toml" ] || ROOT="$SLURM_SUBMIT_DIR"   # sbatch runs a spooled copy
export PATH=/faststorage/project/DeiC-KU-L04/shared_software/pixi/current/bin:$PATH
eval "$(pixi shell-hook --frozen --manifest-path "$ROOT/pixi.toml")"

echo "[$(date -Is)] parsing $OUT with $CORES cores"
Rscript "$ROOT/scripts/parse_modelfinder.R" "$OUT" "$CORES"   # exits 1 if incomplete

mkdir -p "$OUT/archives"
for clade in $(awk -F'\t' 'NR > 1 { print $2 }' "$OUT/jobs.tsv" | sort -u); do
  [ -d "$OUT/$clade" ] || { echo "$clade: already packed"; continue; }
  want=$(( 4 * $(awk -F'\t' -v c="$clade" 'NR > 1 && $2 == c' "$OUT/jobs.tsv" | wc -l) ))
  tar -cf "$OUT/archives/$clade.tar.part" -C "$OUT" "$clade"
  got=$(tar -tf "$OUT/archives/$clade.tar.part" | grep -cv '/$')
  if [ "$got" -ne "$want" ]; then
    echo "$clade: tar has $got files, expected $want; leaving the directory in place" >&2
    exit 1
  fi
  mv "$OUT/archives/$clade.tar.part" "$OUT/archives/$clade.tar"
  rm -rf "${OUT:?}/$clade"
  echo "[$(date -Is)] $clade: $got files -> archives/$clade.tar ($(du -h --apparent-size "$OUT/archives/$clade.tar" | cut -f1))"
done
echo "[$(date -Is)] done"
