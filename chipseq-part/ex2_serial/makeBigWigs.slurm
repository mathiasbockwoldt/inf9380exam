#!/bin/bash

#SBATCH --job-name=mbw
#SBATCH --account=ln0002k
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=00:05:00
#SBATCH --ntasks=1

# get environment and helpers
source /cluster/bin/jobsetup

# Convert bedGraph to bwig files
programs/bedGraphToBigWig -unc N1_INPUT_sorted.bedGraph data/sacisizes.txt N1_INPUT_bigWig.bwig
programs/bedGraphToBigWig -unc T1_INPUT_sorted.bedGraph data/sacisizes.txt T1_INPUT_bigWig.bwig
programs/bedGraphToBigWig -unc T1_CHIP_sorted.bedGraph data/sacisizes.txt T1_CHIP_bigWig.bwig
programs/bedGraphToBigWig -unc N1_CHIP_sorted.bedGraph data/sacisizes.txt N1_CHIP_bigWig.bwig
