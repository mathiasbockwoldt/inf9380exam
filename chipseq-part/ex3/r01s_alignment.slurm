#!/bin/bash

#SBATCH --job-name=excs2_01
#SBATCH --account=ln0002k
#SBATCH --mem-per-cpu=3500M
#SBATCH --time=01:00:00			########################################
#SBATCH --ntasks=1

# # # # # # # # # # # # # # # # # # # #
# This SLURM script starts the alignment step (step 1) of the ChIP-seq variant
# calling workflow. This is the SUBMIT script.
# # # # # # # # # # # # # # # # # # # #

# get environment and helpers
source /cluster/bin/jobsetup

# Aligning reads with Bowtie
arrayrun 0-3 r01e_alignment.slurm
