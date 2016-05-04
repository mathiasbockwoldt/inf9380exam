#!/bin/bash

#SBATCH --job-name=excs2_04
#SBATCH --account=ln0002k
#SBATCH --mem-per-cpu=3500M
#SBATCH --time=01:00:00			########################################
#SBATCH --ntasks=1

# # # # # # # # # # # # # # # # # # # #
# This SLURM script starts the conversion step (step 4) of the ChIP-seq variant
# calling workflow. This is the SUBMIT script.
# # # # # # # # # # # # # # # # # # # #

# get environment and helpers
source /cluster/bin/jobsetup

# Convert sorted bam to bedGraph files
arrayrun 0-3 r04w_makeBedGraphs.slurm
