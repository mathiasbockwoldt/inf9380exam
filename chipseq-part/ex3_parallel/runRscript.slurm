#!/bin/bash

#SBATCH --job-name=r
#SBATCH --account=ln0002k
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=00:10:00
#SBATCH --ntasks=1

# get environment and helpers
source /cluster/bin/jobsetup

# Load modules
module load R

# Normalize and plot in R
Rscript normalizeAndPlot.R
