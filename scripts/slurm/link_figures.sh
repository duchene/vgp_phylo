#!/bin/bash
## Figure sets C, E and P (metric <-> model links, hub variables along the genome,
## path analysis), in dependency order. Needs results/modelfinder/, results/stage0/
## (locus_position, locus_chromosome, locus_softmask) and results/published/.
##   sbatch --account=DeiC-KU-L04 --cpus-per-task=2 --mem=24G --time=01:00:00 \
##     --output=logs/link_figures_%j.out scripts/slurm/link_figures.sh
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd -P)"
[ -f "$ROOT/pixi.toml" ] || ROOT="$SLURM_SUBMIT_DIR"
cd "$ROOT"
export PATH=/faststorage/project/DeiC-KU-L04/shared_software/pixi/current/bin:$PATH
pixi run Rscript scripts/make_link_figures.R      # C1-C4, results/links/hub_selection.tsv
pixi run Rscript scripts/make_hub_landscape.R     # E1-E4
pixi run Rscript scripts/make_path_analysis.R     # P1-P2, results/links/path_coefficients.tsv
