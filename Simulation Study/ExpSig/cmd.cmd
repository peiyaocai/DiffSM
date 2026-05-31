#!/bin/bash
# parallel job using 1 processor and runs for 8 hours (max)
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 2:00:00
#SBATCH --account=change_to_your_account
#SBATCH --mem-per-cpu=8000m

# Execute commands
module load R
R CMD BATCH "--args $arg1 $arg2 $arg3 $arg4 $arg5" '/change_to_your_directory/ExpSig/simu.R'

