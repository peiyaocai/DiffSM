#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 10:00:00
#SBATCH --account=change_to_your_account
#SBATCH --mem-per-cpu=8000m
#SBATCH --output=gaunull_%j.out
#SBATCH --error=gaunull_%j.err

set -euo pipefail

module purge
module load R
module load julia

echo "JobID: $SLURM_JOB_ID on host $(hostname)"
echo "Args: $arg1 $arg2 $arg3 $arg4 $arg5"

echo "Which julia: $(which julia)"
julia --version

echo "Running R..."
Rscript /change_to_your_directory/Gaussian_Compare/Gau_alt.R "$arg1" "$arg2" "$arg3" "$arg4" "$arg5"

echo "Running Julia..."
julia /change_to_your_directory/Gaussian_Compare/GauAlt.jl "$arg1" "$arg2" "$arg3" "$arg4" "$arg5"

echo "Done."
