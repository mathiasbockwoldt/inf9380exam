#!/bin/bash

#SBATCH --job-name=ibed
#SBATCH --account=ln0002k
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=00:05:00
#SBATCH --ntasks=1

# get environment and helpers
source /cluster/bin/jobsetup

# Load modules
module load bedtools

# Run intersectBed to find all reads that intersect with the peaks
intersectBed -bed -wb -abam N1_CHIP_sorted.bam -b N1_T1_merged_peaks.bed > N1_CHIP_intersect.txt
intersectBed -bed -wb -abam T1_CHIP_sorted.bam -b N1_T1_merged_peaks.bed > T1_CHIP_intersect.txt
intersectBed -bed -wb -abam N1_INPUT_sorted.bam -b N1_T1_merged_peaks.bed > N1_INPUT_intersect.txt
intersectBed -bed -wb -abam T1_INPUT_sorted.bam -b N1_T1_merged_peaks.bed > T1_INPUT_intersect.txt
