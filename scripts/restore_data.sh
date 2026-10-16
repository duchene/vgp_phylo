#!/bin/bash
#SBATCH --account=DeiC-KU-L04
#SBATCH --job-name=vgp_restore_data
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=03:00:00
#SBATCH --output=logs/%x_%j.out
## Restore the untracked input data (data/loci, data/gene_trees, data/metadata)
## from the archive into the repository root. Archive paths are repo-relative.
##
## --skip-old-files: never overwrite anything already present, so tracked files
## (e.g. the tip-label-fixed data/species_trees/astralpro.l6.tre) always win.
##
## Usage (from the repo root):  sbatch scripts/restore_data.sh [archive]
set -euo pipefail

ARCHIVE="${1:-../vgp_phylo_untracked.tar.gz}"
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/..}"

echo "[$(date)] extracting $ARCHIVE into $PWD"
tar --skip-old-files --no-same-owner -xzf "$ARCHIVE"
echo "[$(date)] done"

for d in data/loci data/gene_trees; do
  printf '%-16s %7d clades %8d files\n' "$d" \
    "$(find "$d" -mindepth 1 -maxdepth 1 -type d | wc -l)" \
    "$(find "$d" -type f | wc -l)"
done
ls -la data/metadata
