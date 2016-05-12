#!/bin/bash

#SBATCH --job-name=excs2_10
#SBATCH --account=ln0002k
#SBATCH --mem-per-cpu=3500M
#SBATCH --time=01:00:00			########################################
#SBATCH --ntasks=1

# # # # # # # # # # # # # # # # # # # #
# This SLURM script starts the preparation step for R (step 10) of the ChIP-seq
# variant calling workflow. This is the SUBMIT script.
# # # # # # # # # # # # # # # # # # # #

# get environment and helpers
source /cluster/bin/jobsetup

# Get the columns for chr, start, and end for each peak. Sort and count them to find the number of reads that overlap them
# Format to please R (no spaces)
arrayrun 0-3 r09w_intersectBed.slurm
