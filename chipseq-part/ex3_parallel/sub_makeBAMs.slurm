#!/bin/bash

#SBATCH --job-name=mbam
#SBATCH --account=ln0002k
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=00:10:00
#SBATCH --ntasks=1

# get environment and helpers
source /cluster/bin/jobsetup

# Convert sam to bam files
arrayrun 0-3 work_makeBAMs.slurm
