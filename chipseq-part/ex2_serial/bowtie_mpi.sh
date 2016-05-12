#!/bin/bash

mpirun -np 8 pmap_dist . . data/${1}.fastq -i data/ saci
mpirun -np 8 pmap -i data/ saci . . bowtie -S -v 2 -m 1
awk 'NR <= 3 || !/^@/' < out.txt > ${1}.sam
rm out.txt
