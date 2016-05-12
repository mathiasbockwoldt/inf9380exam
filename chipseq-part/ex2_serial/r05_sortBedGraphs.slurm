#!/bin/bash

#SBATCH --job-name=excs2_05
#SBATCH --account=ln0002k
#SBATCH --mem-per-cpu=3500M
#SBATCH --time=01:00:00			########################################
#SBATCH --ntasks=1

# # # # # # # # # # # # # # # # # # # #
# This SLURM script starts the bedgraph sorting step (step 5) of the ChIP-seq
# variant calling workflow.
# # # # # # # # # # # # # # # # # # # #

# get environment and helpers
source /cluster/bin/jobsetup

# Sort bedgraph files
sort -k1,1 -k2,2n N1_INPUT.bedGraph > N1_INPUT_sorted.bedGraph
sort -k1,1 -k2,2n T1_INPUT.bedGraph > T1_INPUT_sorted.bedGraph
sort -k1,1 -k2,2n T1_CHIP.bedGraph > T1_CHIP_sorted.bedGraph
sort -k1,1 -k2,2n N1_CHIP.bedGraph > N1_CHIP_sorted.bedGraph
