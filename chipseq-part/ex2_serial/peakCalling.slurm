#!/bin/bash

#SBATCH --job-name=peakc
#SBATCH --account=ln0002k
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=00:05:00
#SBATCH --ntasks=1

# get environment and helpers
source /cluster/bin/jobsetup

# Load modules
module load macs

# Call peaks using macs
macs -t N1_CHIP.sam -c N1_INPUT.sam -n N1 -g 2.23e6 -f SAM	######################
macs -t T1_CHIP.sam -c T1_INPUT.sam -n T1 -g 2.23e6 -f SAM  ### The number must be dynamic
