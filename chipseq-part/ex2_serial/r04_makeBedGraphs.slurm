#!/bin/bash

#SBATCH --job-name=mbg
#SBATCH --account=ln0002k
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=00:05:00
#SBATCH --ntasks=1

# get environment and helpers
source /cluster/bin/jobsetup

# Load modules
module load bedtools

# Convert sorted bam to bedGraph files
genomeCoverageBed -ibam N1_CHIP_sorted.bam -bg -g data/sacisizes.txt > N1_CHIP.bedGraph
genomeCoverageBed -ibam T1_CHIP_sorted.bam -bg -g data/sacisizes.txt > T1_CHIP.bedGraph
genomeCoverageBed -ibam T1_INPUT_sorted.bam -bg -g data/sacisizes.txt > T1_INPUT.bedGraph
genomeCoverageBed -ibam N1_INPUT_sorted.bam -bg -g data/sacisizes.txt > N1_INPUT.bedGraph
