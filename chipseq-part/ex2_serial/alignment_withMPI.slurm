#!/bin/bash

#SBATCH --job-name=alignmpi
#SBATCH --account=ln0002k
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=00:10:00
#SBATCH --ntasks=1

# get environment and helpers
source /cluster/bin/jobsetup

# Load modules
module load bowtie
module load pmap

# Aligning reads with Bowtie
./bowtie_mpi.sh N1_CHIP
./bowtie_mpi.sh N1_INPUT
./bowtie_mpi.sh T1_CHIP
./bowtie_mpi.sh T1_INPUT
