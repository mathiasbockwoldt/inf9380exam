#!/bin/bash

#SBATCH --job-name=excs2_06
#SBATCH --account=ln0002k
#SBATCH --mem-per-cpu=3500M
#SBATCH --time=01:00:00			########################################
#SBATCH --ntasks=1

# # # # # # # # # # # # # # # # # # # #
# This SLURM script starts the conversion step (step 6) of the ChIP-seq variant
# calling workflow. This is the SUBMIT script.
# # # # # # # # # # # # # # # # # # # #

# get environment and helpers
source /cluster/bin/jobsetup

# Convert bedGraph to bwig files
arrayrun 0-3 r06w_makeBigWigs.slurm
