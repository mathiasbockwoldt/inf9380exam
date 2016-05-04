#!/bin/bash

#SBATCH --job-name=excs2_06_w
#SBATCH --account=ln0002k
#SBATCH --mem-per-cpu=3500M
#SBATCH --time=01:00:00			########################################
#SBATCH --ntasks=1

# # # # # # # # # # # # # # # # # # # #
# This SLURM script starts the conversion step (step 6) of the ChIP-seq variant
# calling workflow. This is the WORKER script.
# # # # # # # # # # # # # # # # # # # #

# get environment and helpers
source /cluster/bin/jobsetup

declare -a fnames=(N1_CHIP N1_INPUT T1_CHIP T1_INPUT)

# Convert bedGraph to bwig files
programs/bedGraphToBigWig -unc ${fnames[$TASK_ID]}_sorted.bedGraph data/sacisizes.txt ${fnames[$TASK_ID]}_bigWig.bwig
