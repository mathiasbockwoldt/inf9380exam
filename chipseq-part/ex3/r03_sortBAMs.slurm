#!/bin/bash

#SBATCH --job-name=excs2_03
#SBATCH --account=ln0002k
#SBATCH --mem-per-cpu=3500M
#SBATCH --time=01:00:00			########################################
#SBATCH --ntasks=1

# # # # # # # # # # # # # # # # # # # #
# This SLURM script starts the sorting step (step 3) of the ChIP-seq variant
# calling workflow.
# # # # # # # # # # # # # # # # # # # #

# get environment and helpers
source /cluster/bin/jobsetup

# Load modules
module load samtools

# Sort bam files
samtools sort N1_INPUT.bam N1_INPUT_sorted
samtools sort T1_INPUT.bam T1_INPUT_sorted
samtools sort T1_CHIP.bam T1_CHIP_sorted
samtools sort N1_CHIP.bam N1_CHIP_sorted
