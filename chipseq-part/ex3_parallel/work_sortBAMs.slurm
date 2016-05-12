#!/bin/bash

#SBATCH --job-name=sbam_w
#SBATCH --account=ln0002k
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=00:10:00
#SBATCH --ntasks=1

# get environment and helpers
source /cluster/bin/jobsetup

# Load modules
module load samtools

declare -a fnames=(N1_CHIP N1_INPUT T1_CHIP T1_INPUT)

# Sort bam files
samtools sort ${fnames[$TASK_ID]}.bam ${fnames[$TASK_ID]}_sorted
