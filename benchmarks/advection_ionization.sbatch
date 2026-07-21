#!/bin/bash
#SBATCH --job-name=cmac-advection
#SBATCH --time=00:10:00
#SBATCH --mem=2G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --output=advection-ionization-%j.out

set -euo pipefail

threads="${SLURM_CPUS_PER_TASK:-1}"
srun --ntasks=1 ../../CMacIonize --task-based-rhd \
  --params advection_ionization.param --threads "${threads}"
python3 advection_ionization.py
