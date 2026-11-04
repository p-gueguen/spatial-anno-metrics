#!/bin/bash
#SBATCH --job-name=sam_scdesign3_fit
#SBATCH --partition=employee
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G
#SBATCH --time=06:00:00
#SBATCH --output=/srv/GT/analysis/pgueguen/sam_dev/runs/sim-scdesign3/fit_%j.log
# The one expensive step: the scDesign3 fit (generate.R). Run `python build.py prepare` first and
# `python build.py pack` after; submit from this folder. CACHE, K_GP and HERE can be overridden.
set -euo pipefail
CACHE=${CACHE:-/srv/GT/analysis/pgueguen/sam_dev/data/sim_worlds/scdesign3}
K_GP=${K_GP:-100}
HERE=${HERE:-$SLURM_SUBMIT_DIR}   # submit from benchmarks/worlds/scdesign3
source /etc/profile.d/lmod.sh 2>/dev/null || true
module load Dev/R/4.6.0
export R_LIBS_USER=/srv/GT/analysis/pgueguen/sam_dev/Rlib
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
/usr/bin/time -v Rscript "$HERE/generate.R" "$CACHE" "$SLURM_CPUS_PER_TASK" "$K_GP"
