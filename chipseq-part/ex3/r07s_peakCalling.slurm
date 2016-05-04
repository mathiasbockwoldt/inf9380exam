#!/bin/bash

#SBATCH --job-name=excs2_07
#SBATCH --account=ln0002k
#SBATCH --mem-per-cpu=3500M
#SBATCH --time=01:00:00			########################################
#SBATCH --ntasks=1

# # # # # # # # # # # # # # # # # # # #
# This SLURM script starts the peak calling step (step 7) of the ChIP-seq
# variant calling workflow. This is the SUBMIT script.
# # # # # # # # # # # # # # # # # # # #

# get environment and helpers
source /cluster/bin/jobsetup

# Call peaks using macs
arrayrun 0-1 r07w_pealCalling.slurm
