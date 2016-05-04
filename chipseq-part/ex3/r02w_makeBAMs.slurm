#!/bin/bash

#SBATCH --job-name=excs2_02_w
#SBATCH --account=ln0002k
#SBATCH --mem-per-cpu=3500M
#SBATCH --time=01:00:00			########################################
#SBATCH --ntasks=1

# # # # # # # # # # # # # # # # # # # #
# This SLURM script starts the conversion step (step 2) of the ChIP-seq variant
# calling workflow. This is the WORKER script.
# # # # # # # # # # # # # # # # # # # #

# get environment and helpers
source /cluster/bin/jobsetup

# Load modules
module load samtools

declare -a fnames=(N1_CHIP N1_INPUT T1_CHIP T1_INPUT)

# Convert sam to bam files
samtools view ${fnames[$TASK_ID]}.sam -bS > ${fnames[$TASK_ID]}.bam
