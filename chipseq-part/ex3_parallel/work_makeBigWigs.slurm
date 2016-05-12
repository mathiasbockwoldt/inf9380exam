#!/bin/bash

#SBATCH --job-name=mbw_w
#SBATCH --account=ln0002k
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=00:10:00
#SBATCH --ntasks=1

# get environment and helpers
source /cluster/bin/jobsetup

declare -a fnames=(N1_CHIP N1_INPUT T1_CHIP T1_INPUT)

# Convert bedGraph to bwig files
programs/bedGraphToBigWig -unc ${fnames[$TASK_ID]}_sorted.bedGraph data/sacisizes.txt ${fnames[$TASK_ID]}_bigWig.bwig
