#!/bin/bash

#SBATCH --job-name=ex_cs_1
#SBATCH --account=ln0002k
#SBATCH --mem-per-cpu=256M
#SBATCH --time=00:12:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16

# # # # # # # # # # # # # # # # # # # #
# This SLURM script starts bowtie in serial and in parallel using MPI multiple
# times to get an idea of how good the speed-up is.
# I'm aware of that the design of this script is not best practice, but it get's
# its work done.
# # # # # # # # # # # # # # # # # # # #

# get environment and helpers
source /cluster/bin/jobsetup

# Load modules
module load bowtie
module load pmap

# First, we run bowtie five times in serial mode
echo "Plain bowtie"
time bowtie -S -v 2 -m 1 data/saci data/N1_CHIP.fastq > N1_CHIP_nompi.sam
time bowtie -S -v 2 -m 1 data/saci data/N1_CHIP.fastq > N1_CHIP_nompi.sam
time bowtie -S -v 2 -m 1 data/saci data/N1_CHIP.fastq > N1_CHIP_nompi.sam
time bowtie -S -v 2 -m 1 data/saci data/N1_CHIP.fastq > N1_CHIP_nompi.sam
time bowtie -S -v 2 -m 1 data/saci data/N1_CHIP.fastq > N1_CHIP_nompi.sam

# Second, we run bowtie in parallel using MPI with 1 to 16 cores five times each
for i in {1..16}
do
	echo "Num Processors: $i"
	time ./mpirun_script.sh $i
	time ./mpirun_script.sh $i
	time ./mpirun_script.sh $i
	time ./mpirun_script.sh $i
	time ./mpirun_script.sh $i
done
