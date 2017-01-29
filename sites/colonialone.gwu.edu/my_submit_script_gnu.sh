#!/bin/bash

#SBATCH -o ib_rotcyl_black_ug%j.out
#SBATCH -e ib_rotcyl_black_ug%j.err
#SBATCH -p defq -n 16
# SBATCH -N 1
# SBATCH -n 8
#SBATCH -D /lustre/groups/balarasgrp/mvanella/IB_ROTCYL_UG_BLACKDBG
#SBATCH -J Fstflash_ts
#SBATCH --export=NONE
#SBATCH -t 0:30:00
# SBATCH --mem-per-cpu=1000000
# SBATCH --array=1-16
#SBATCH --nice=100

module load openmpi/gcc/64/1.8

mpirun ./flash4

