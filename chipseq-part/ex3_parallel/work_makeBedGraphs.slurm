#!/bin/bash

#SBATCH --job-name=mbg_w
#SBATCH --account=ln0002k
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=00:10:00
#SBATCH --ntasks=1

# get environment and helpers
source /cluster/bin/jobsetup

# Load modules
module load bedtools

declare -a fnames=(N1_CHIP N1_INPUT T1_CHIP T1_INPUT)

# Convert sorted bam to bedGraph files
genomeCoverageBed -ibam ${fnames[$TASK_ID]}_sorted.bam -bg -g data/sacisizes.txt > ${fnames[$TASK_ID]}.bedGraph
