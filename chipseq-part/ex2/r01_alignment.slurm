#!/bin/bash

#SBATCH --job-name=excs2_01
#SBATCH --account=ln0002k
#SBATCH --mem-per-cpu=3500M
#SBATCH --time=01:00:00			########################################
#SBATCH --ntasks=1

# # # # # # # # # # # # # # # # # # # #
# This SLURM script starts the alignment step (step 1) of the ChIP-seq variant
# calling workflow.
# # # # # # # # # # # # # # # # # # # #

# get environment and helpers
source /cluster/bin/jobsetup

# Load modules
module load bowtie

# Aligning reads with Bowtie
bowtie -S -v 2 -m 1 data/saci data/N1_CHIP.fastq > N1_CHIP.sam
bowtie -S -v 2 -m 1 data/saci data/N1_INPUT.fastq > N1_INPUT.sam
bowtie -S -v 2 -m 1 data/saci data/T1_CHIP.fastq > T1_CHIP.sam
bowtie -S -v 2 -m 1 data/saci data/T1_INPUT.fastq > T1_INPUT.sam
