#!/bin/bash

mpirun -np $1 pmap_dist . . data/N1_CHIP.fastq -i data/ saci
mpirun -np $1 pmap -i data/ saci . . bowtie -S -v 2 -m 1
mv out.txt N1_CHIP.sam
