#!/bin/bash

#SBATCH --job-name=excs2_04
#SBATCH --account=ln0002k
#SBATCH --mem-per-cpu=3500M
#SBATCH --time=01:00:00			########################################
#SBATCH --ntasks=1

# # # # # # # # # # # # # # # # # # # #
# This SLURM script starts the conversion step (step 4) of the ChIP-seq variant
# calling workflow.
# # # # # # # # # # # # # # # # # # # #

# get environment and helpers
source /cluster/bin/jobsetup

# Load modules
module load bedtools

# Convert sorted bam to bedGraph files
genomeCoverageBed -ibam N1_CHIP_sorted.bam -bg -g data/sacisizes.txt > N1_CHIP.bedGraph
genomeCoverageBed -ibam T1_CHIP_sorted.bam -bg -g data/sacisizes.txt > T1_CHIP.bedGraph
genomeCoverageBed -ibam T1_INPUT_sorted.bam -bg -g data/sacisizes.txt > T1_INPUT.bedGraph
genomeCoverageBed -ibam N1_INPUT_sorted.bam -bg -g data/sacisizes.txt > N1_INPUT.bedGraph
