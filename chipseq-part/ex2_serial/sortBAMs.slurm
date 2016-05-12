#!/bin/bash

#SBATCH --job-name=sbam
#SBATCH --account=ln0002k
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=00:05:00
#SBATCH --ntasks=1

# get environment and helpers
source /cluster/bin/jobsetup

# Load modules
module load samtools

# Sort bam files
samtools sort N1_INPUT.bam N1_INPUT_sorted
samtools sort T1_INPUT.bam T1_INPUT_sorted
samtools sort T1_CHIP.bam T1_CHIP_sorted
samtools sort N1_CHIP.bam N1_CHIP_sorted
