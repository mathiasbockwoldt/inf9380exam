#!/bin/bash

#SBATCH --job-name=preparation
#SBATCH --account=ln0002k
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=00:05:00
#SBATCH --ntasks=1

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
