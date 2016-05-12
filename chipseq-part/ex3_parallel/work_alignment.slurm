#!/bin/bash

#SBATCH --job-name=align_w
#SBATCH --account=ln0002k
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=00:10:00
#SBATCH --ntasks=1

# get environment and helpers
source /cluster/bin/jobsetup

# Load modules
module load bowtie

declare -a fnames=(N1_CHIP N1_INPUT T1_CHIP T1_INPUT)

# Aligning reads with Bowtie
bowtie -S -v 2 -m 1 data/saci data/${fnames[$TASK_ID]}.fastq > ${fnames[$TASK_ID]}.sam
