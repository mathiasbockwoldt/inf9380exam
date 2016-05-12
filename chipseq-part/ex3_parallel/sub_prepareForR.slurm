#!/bin/bash

#SBATCH --job-name=prepr
#SBATCH --account=ln0002k
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=00:10:00
#SBATCH --ntasks=1

# get environment and helpers
source /cluster/bin/jobsetup

# Get the columns for chr, start, and end for each peak. Sort and count them to find the number of reads that overlap them
# Format to please R (no spaces)
arrayrun 0-3 work_intersectBed.slurm
