#!/bin/bash

#SBATCH --job-name=ibed_w
#SBATCH --account=ln0002k
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=00:10:00
#SBATCH --ntasks=1

# get environment and helpers
source /cluster/bin/jobsetup

# Load modules
module load bedtools

declare -a fnames=(N1_CHIP N1_INPUT T1_CHIP T1_INPUT)

# Run intersectBed to find all reads that intersect with the peaks
intersectBed -bed -wb -abam ${fnames[$TASK_ID]}_sorted.bam -b N1_T1_merged_peaks.bed > ${fnames[$TASK_ID]}_intersect.txt
