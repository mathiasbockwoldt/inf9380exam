#!/bin/bash

#SBATCH --job-name=excs2_09_w
#SBATCH --account=ln0002k
#SBATCH --mem-per-cpu=3500M
#SBATCH --time=01:00:00			########################################
#SBATCH --ntasks=1

# # # # # # # # # # # # # # # # # # # #
# This SLURM script starts the intersection step (step 9) of the ChIP-seq
# variant calling workflow. This is the WORKER script.
# # # # # # # # # # # # # # # # # # # #

# get environment and helpers
source /cluster/bin/jobsetup

# Load modules
module load bedtools

declare -a fnames=(N1_CHIP N1_INPUT T1_CHIP T1_INPUT)

# Run intersectBed to find all reads that intersect with the peaks
intersectBed -bed -wb -abam ${fnames[$TASK_ID]}_sorted.bam -b N1_T1_merged_peaks.bed > ${fnames[$TASK_ID]}_intersect.txt
