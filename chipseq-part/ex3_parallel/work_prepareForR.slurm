#!/bin/bash

#SBATCH --job-name=prepr_w
#SBATCH --account=ln0002k
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=00:10:00
#SBATCH --ntasks=1

# get environment and helpers
source /cluster/bin/jobsetup

declare -a fnames=(N1_CHIP N1_INPUT T1_CHIP T1_INPUT)

# Get the columns for chr, start, and end for each peak. Sort and count them to find the number of reads that overlap them
# Format to please R (no spaces)
cut -f13,14,15 ${fnames[$TASK_ID]}_intersect.txt | sort | uniq -c | awk '{print $2"_"$3"_"$4"\t"$1}' > ${fnames[$TASK_ID]}_countR.txt
