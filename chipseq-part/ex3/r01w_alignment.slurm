#!/bin/bash

#SBATCH --job-name=excs2_01_w
#SBATCH --account=ln0002k
#SBATCH --mem-per-cpu=3500M
#SBATCH --time=01:00:00			########################################
#SBATCH --ntasks=1

# # # # # # # # # # # # # # # # # # # #
# This SLURM script starts the alignment step (step 1) of the ChIP-seq variant
# calling workflow. This is the WORKER script.
# # # # # # # # # # # # # # # # # # # #

# get environment and helpers
source /cluster/bin/jobsetup

# Load modules
module load bowtie

declare -a fnames=(N1_CHIP N1_INPUT T1_CHIP T1_INPUT)

# Aligning reads with Bowtie
bowtie -S -v 2 -m 1 data/saci data/${fnames[$TASK_ID]}.fastq > ${fnames[$TASK_ID]}.sam
