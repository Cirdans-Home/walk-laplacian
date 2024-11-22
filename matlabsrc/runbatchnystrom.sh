#!/bin/bash
#SBATCH --job-name="NYSTAVG"
#SBATCH --partition=gpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=256
#SBATCH --time=2-00:00:00
#SBATCH -o "nystrom.log"
#SBATCH -e "nystrom.err"

module purge
module load matlab
module list

matlab -batch "batchnystrom"
