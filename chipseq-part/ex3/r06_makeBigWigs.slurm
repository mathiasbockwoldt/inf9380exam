#!/bin/bash

#SBATCH --job-name=excs2_06
#SBATCH --account=ln0002k
#SBATCH --mem-per-cpu=3500M
#SBATCH --time=01:00:00			########################################
#SBATCH --ntasks=1

# # # # # # # # # # # # # # # # # # # #
# This SLURM script starts the conversion step (step 6) of the ChIP-seq variant
# calling workflow.
# # # # # # # # # # # # # # # # # # # #

# get environment and helpers
source /cluster/bin/jobsetup

# Convert bedGraph to bwig files
programs/bedGraphToBigWig -unc N1_INPUT_sorted.bedGraph data/sacisizes.txt N1_INPUT_bigWig.bwig
programs/bedGraphToBigWig -unc T1_INPUT_sorted.bedGraph data/sacisizes.txt T1_INPUT_bigWig.bwig
programs/bedGraphToBigWig -unc T1_CHIP_sorted.bedGraph data/sacisizes.txt T1_CHIP_bigWig.bwig
programs/bedGraphToBigWig -unc N1_CHIP_sorted.bedGraph data/sacisizes.txt N1_CHIP_bigWig.bwig
