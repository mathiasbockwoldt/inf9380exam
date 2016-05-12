#!/bin/bash

#SBATCH --job-name=sbg
#SBATCH --account=ln0002k
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=01:00:00
#SBATCH --ntasks=1

# get environment and helpers
source /cluster/bin/jobsetup

# Sort bedgraph files
arrayrun 0-3 work_sortBedGraphs.slurm
