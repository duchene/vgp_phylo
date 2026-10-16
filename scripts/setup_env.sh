#!/bin/bash
#SBATCH --account=DeiC-KU-L04
#SBATCH --job-name=vgp_setup_env
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x_%j.out
## Build the pixi R environment (.pixi/) on a compute node, as the shared pixi
## README asks, then install the CRAN-only packages. Idempotent.
##
## Usage (from the repo root):  sbatch scripts/setup_env.sh
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/..}"

export PATH=/faststorage/project/DeiC-KU-L04/shared_software/pixi/current/bin:$PATH
## Prefer the shared package cache, but pixi writes its cache files mode 0600, so
## another member's entries can be unreadable to you ("Permission denied (os error
## 13)" while solving). Fall back to pixi's per-user default cache in that case.
SHARED_CACHE=/faststorage/project/DeiC-KU-L04/shared_software/pixi/cache
if [ -z "$(find "$SHARED_CACHE" ! -readable -print -quit 2>/dev/null)" ]; then
  export PIXI_CACHE_DIR="$SHARED_CACHE"
else
  echo "shared pixi cache has files you cannot read; using the per-user cache" >&2
  unset PIXI_CACHE_DIR
fi

pixi install
pixi run install-cran
pixi run Rscript -e 'sessionInfo()'
