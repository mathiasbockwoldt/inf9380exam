#!/bin/bash

#SBATCH --job-name=prepr
#SBATCH --account=ln0002k
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=00:05:00
#SBATCH --ntasks=1

# get environment and helpers
source /cluster/bin/jobsetup

# Get the columns for chr, start, and end for each peak. Sort and count them to find the number of reads that overlap them
# Format to please R (no spaces)
cut -f13,14,15 N1_CHIP_intersect.txt | sort | uniq -c | awk '{print $2"_"$3"_"$4"\t"$1}' > N1_CHIP_countR.txt
cut -f13,14,15 T1_CHIP_intersect.txt | sort | uniq -c | awk '{print $2"_"$3"_"$4"\t"$1}' > T1_CHIP_countR.txt
cut -f13,14,15 T1_INPUT_intersect.txt | sort | uniq -c | awk '{print $2"_"$3"_"$4"\t"$1}' > T1_INPUT_countR.txt
cut -f13,14,15 N1_INPUT_intersect.txt | sort | uniq -c | awk '{print $2"_"$3"_"$4"\t"$1}' > N1_INPUT_countR.txt
