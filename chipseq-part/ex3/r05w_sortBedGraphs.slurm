#!/bin/bash

#SBATCH --job-name=excs2_05_w
#SBATCH --account=ln0002k
#SBATCH --mem-per-cpu=3500M
#SBATCH --time=01:00:00			########################################
#SBATCH --ntasks=1

# # # # # # # # # # # # # # # # # # # #
# This SLURM script starts the bedgraph sorting step (step 5) of the ChIP-seq
# variant calling workflow. This is the WORKER script.
# # # # # # # # # # # # # # # # # # # #

# get environment and helpers
source /cluster/bin/jobsetup

declare -a fnames=(N1_CHIP N1_INPUT T1_CHIP T1_INPUT)

# Sort bedgraph files
sort -k1,1 -k2,2n ${fnames[$TASK_ID]}.bedGraph > ${fnames[$TASK_ID]}_sorted.bedGraph
