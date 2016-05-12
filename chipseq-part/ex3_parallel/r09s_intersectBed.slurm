#!/bin/bash

#SBATCH --job-name=excs2_09
#SBATCH --account=ln0002k
#SBATCH --mem-per-cpu=3500M
#SBATCH --time=01:00:00			########################################
#SBATCH --ntasks=1

# # # # # # # # # # # # # # # # # # # #
# This SLURM script starts the intersection step (step 9) of the ChIP-seq
# variant calling workflow. This is the SUBMIT script.
# # # # # # # # # # # # # # # # # # # #

# get environment and helpers
source /cluster/bin/jobsetup

# Run intersectBed to find all reads that intersect with the peaks
arrayrun 0-3 r09w_intersectBed.slurm
