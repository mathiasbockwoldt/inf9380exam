#!/bin/bash

#SBATCH --job-name=excs2_00
#SBATCH --account=ln0002k
#SBATCH --mem-per-cpu=3500M
#SBATCH --time=01:00:00			########################################
#SBATCH --ntasks=1

# # # # # # # # # # # # # # # # # # # #
# This SLURM script starts the preparation step (step 0) of the ChIP-seq variant
# calling workflow.
# # # # # # # # # # # # # # # # # # # #

# get environment and helpers
source /cluster/bin/jobsetup

# Load modules
module load bowtie

# Copy important files and programs
cp -r /work/projects/norbis/workflows/chip/data .
cp -r /work/projects/norbis/workflows/chip/programs .

# Add execution rights to programs
chmod a+rx programs/*

# Indexing the reference genome
bowtie-build -f data/saci.fasta data/saci
