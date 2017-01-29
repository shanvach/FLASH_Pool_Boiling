#!/bin/bash
#PBS -l nodes=64:nehalem
#PBS -j oe
#PBS -o output_IBsphere_INTEL_64h.txt
#PBS -l walltime=1:00:00
source /home/balaras/marcos/My_Setup.sh
cd $PBS_O_WORKDIR

export OPENBLAS_NUM_THREADS=1
export OPENMP_NUM_THREADS=1

mpiexec -n 64 ./flash4 > case_IBsphere_INTEL_64h.scr 
