#!/bin/bash

#SBATCH --job-name=excs2_04_w
#SBATCH --account=ln0002k
#SBATCH --mem-per-cpu=3500M
#SBATCH --time=01:00:00			########################################
#SBATCH --ntasks=1

# # # # # # # # # # # # # # # # # # # #
# This SLURM script starts the conversion step (step 4) of the ChIP-seq variant
# calling workflow. This is the WORKER script.
# # # # # # # # # # # # # # # # # # # #

# get environment and helpers
source /cluster/bin/jobsetup

# Load modules
module load bedtools

declare -a fnames=(N1_CHIP N1_INPUT T1_CHIP T1_INPUT)

# Convert sorted bam to bedGraph files
genomeCoverageBed -ibam ${fnames[$TASK_ID]}_sorted.bam -bg -g data/sacisizes.txt > ${fnames[$TASK_ID]}.bedGraph
