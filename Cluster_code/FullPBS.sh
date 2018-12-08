#!/bin/bash

#PBS -A open
#PBS -l nodes=1:ppn=1
#PBS -l walltime=48:00:00
#PBS -l pmem=100gb
#PBS -N fullPenalized

cd $PBS_O_WORKDIR

module load r/3.4

Rscript chr1_lasso_overall_analysis_v8_clust.R
