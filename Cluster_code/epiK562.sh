#!/bin/bash

#PBS -A open
#PBS -l nodes=1:ppn=1
#PBS -l walltime=48:00:00
#PBS -l pmem=2gb
#PBS -N epigen_K562

cd $PBS_O_WORKDIR

module load r/3.4

Rscript epigen_preprocessing_K562_clust.R
