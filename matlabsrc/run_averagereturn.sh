#!/bin/bash
#SBATCH --nodes=1
#SBATCH --partition=gpu
#SBATCH --time=30-00:00:00
#SBATCH --exclusive
#SBATCH --job-name=AVRET
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=256
#SBATCH --error=batch-average-error.log
#SBATCH --output=batch-average-output.log


module purge
module load matlab/R2021a
module list

matlab -batch averagereturnprobabilities
