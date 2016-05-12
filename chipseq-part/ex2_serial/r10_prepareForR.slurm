#!/bin/bash

#SBATCH --job-name=excs2_10
#SBATCH --account=ln0002k
#SBATCH --mem-per-cpu=3500M
#SBATCH --time=01:00:00			########################################
#SBATCH --ntasks=1

# # # # # # # # # # # # # # # # # # # #
# This SLURM script starts the preparation step for R (step 10) of the ChIP-seq
# variant calling workflow.
# # # # # # # # # # # # # # # # # # # #

# get environment and helpers
source /cluster/bin/jobsetup

# Get the columns for chr, start, and end for each peak. Sort and count them to find the number of reads that overlap them
# Format to please R (no spaces)
cut -f13,14,15 N1_CHIP_intersect.txt | sort | uniq -c | awk '{print $2"_"$3"_"$4"\t"$1}' > N1_CHIP_countR.txt
cut -f13,14,15 T1_CHIP_intersect.txt | sort | uniq -c | awk '{print $2"_"$3"_"$4"\t"$1}' > T1_CHIP_countR.txt
cut -f13,14,15 T1_INPUT_intersect.txt | sort | uniq -c | awk '{print $2"_"$3"_"$4"\t"$1}' > T1_INPUT_countR.txt
cut -f13,14,15 N1_INPUT_intersect.txt | sort | uniq -c | awk '{print $2"_"$3"_"$4"\t"$1}' > N1_INPUT_countR.txt
