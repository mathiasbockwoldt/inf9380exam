#!/bin/bash

#SBATCH --job-name=excs2_07_w
#SBATCH --account=ln0002k
#SBATCH --mem-per-cpu=3500M
#SBATCH --time=01:00:00			########################################
#SBATCH --ntasks=1

# # # # # # # # # # # # # # # # # # # #
# This SLURM script starts the peak calling step (step 7) of the ChIP-seq
# variant calling workflow. This is the WORKER script.
# # # # # # # # # # # # # # # # # # # #

# get environment and helpers
source /cluster/bin/jobsetup

# Load modules
module load macs

declare -a in1=(N1_CHIP T1_CHIP)
declare -a in2=(N1_INPUT T1_INPUT)
declare -a out=(N1 T1)

# Call peaks using macs
macs -t ${in1[$TASK_ID]}.sam -c ${in2[$TASK_ID]}.sam -n ${out[$TASK_ID]} -g 2.23e6 -f SAM
########## The number must be dynamic
