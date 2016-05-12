#!/bin/bash

#SBATCH --job-name=ibed
#SBATCH --account=ln0002k
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=00:10:00
#SBATCH --ntasks=1

# get environment and helpers
source /cluster/bin/jobsetup

# Run intersectBed to find all reads that intersect with the peaks
arrayrun 0-3 work_intersectBed.slurm
