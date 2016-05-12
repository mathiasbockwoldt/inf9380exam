#!/bin/bash

#SBATCH --job-name=mbam
#SBATCH --account=ln0002k
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=00:05:00
#SBATCH --ntasks=1

# get environment and helpers
source /cluster/bin/jobsetup

# Load modules
module load samtools

# Convert sam to bam files
samtools view N1_CHIP.sam -bS > N1_CHIP.bam
samtools view T1_CHIP.sam -bS > T1_CHIP.bam
samtools view N1_INPUT.sam -bS > N1_INPUT.bam
samtools view T1_INPUT.sam -bS > T1_INPUT.bam
