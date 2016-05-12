#!/bin/bash

#SBATCH --job-name=csmerge
#SBATCH --account=ln0002k
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=00:10:00
#SBATCH --ntasks=1

# get environment and helpers
source /cluster/bin/jobsetup

# Load modules
module load bedtools

# Get one file with all peaks and sort it
cat N1_peaks.bed T1_peaks.bed > N1_T1_peaks.bed
sort -k1,1 -k2,2n N1_T1_peaks.bed > N1_T1_sorted_peaks.bed

# Run bedtools merge for peaks found in N1 and T1 to have the complete set for comparison
bedtools merge -i N1_T1_sorted_peaks.bed > N1_T1_merged_peaks.bed
