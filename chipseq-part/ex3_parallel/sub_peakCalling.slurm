#!/bin/bash

#SBATCH --job-name=peakc
#SBATCH --account=ln0002k
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=00:10:00
#SBATCH --ntasks=1

# get environment and helpers
source /cluster/bin/jobsetup

# Call peaks using macs
arrayrun 0-1 work_peakCalling.slurm
