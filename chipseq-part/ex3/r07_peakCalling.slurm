#!/bin/bash

#SBATCH --job-name=excs2_07
#SBATCH --account=ln0002k
#SBATCH --mem-per-cpu=3500M
#SBATCH --time=01:00:00			########################################
#SBATCH --ntasks=1

# # # # # # # # # # # # # # # # # # # #
# This SLURM script starts the peak calling step (step 7) of the ChIP-seq
# variant calling workflow.
# # # # # # # # # # # # # # # # # # # #

# get environment and helpers
source /cluster/bin/jobsetup

# Load modules
module load macs

# Call peaks using macs
macs -t N1_CHIP.sam -c N1_INPUT.sam -n N1 -g 2.23e6 -f SAM	######################
macs -t T1_CHIP.sam -c T1_INPUT.sam -n T1 -g 2.23e6 -f SAM  ### The number must be dynamic
