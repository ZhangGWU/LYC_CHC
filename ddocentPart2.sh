#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=gompert-kp
#SBATCH --partition=gompert-kp
#SBATCH --job-name=Aqlddo2
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=linyi.zhang@usu.edu

mawk -v x=4 '$1 >= x' uniqCperindv > uniq.k.4.c.4.seqs

cut -f2 uniq.k.4.c.4.seqs > totaluniqseq

mawk '{c= c + 1; print ">scaffold_" c "\n" $1}' totaluniqseq > uniq.fasta

cd-hit-est -i uniq.fasta -o reference.fasta -M 0 -T 0 -c 0.8

