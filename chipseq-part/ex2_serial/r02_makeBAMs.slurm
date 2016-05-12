#!/bin/bash

#SBATCH --job-name=excs2_02
#SBATCH --account=ln0002k
#SBATCH --mem-per-cpu=3500M
#SBATCH --time=01:00:00			########################################
#SBATCH --ntasks=1

# # # # # # # # # # # # # # # # # # # #
# This SLURM script starts the conversion step (step 2) of the ChIP-seq variant
# calling workflow.
# # # # # # # # # # # # # # # # # # # #

# get environment and helpers
source /cluster/bin/jobsetup

# Load modules
module load samtools

# Convert sam to bam files
samtools view N1_CHIP.sam -bS > N1_CHIP.bam
samtools view T1_CHIP.sam -bS > T1_CHIP.bam
samtools view N1_INPUT.sam -bS > N1_INPUT.bam
samtools view T1_INPUT.sam -bS > T1_INPUT.bam
