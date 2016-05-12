#!/bin/bash

#SBATCH --job-name=sbg
#SBATCH --account=ln0002k
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=00:05:00
#SBATCH --ntasks=1

# get environment and helpers
source /cluster/bin/jobsetup

# Sort bedgraph files
sort -k1,1 -k2,2n N1_INPUT.bedGraph > N1_INPUT_sorted.bedGraph
sort -k1,1 -k2,2n T1_INPUT.bedGraph > T1_INPUT_sorted.bedGraph
sort -k1,1 -k2,2n T1_CHIP.bedGraph > T1_CHIP_sorted.bedGraph
sort -k1,1 -k2,2n N1_CHIP.bedGraph > N1_CHIP_sorted.bedGraph
