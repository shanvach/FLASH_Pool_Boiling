#!/bin/bash

#SBATCH -o fs_%j.out
#SBATCH -e fs_%j.err
#SBATCH -p short -n 8
# SBATCH -p debug -n 8
# SBATCH -p defq -n 16
# SBATCH -N 1
# SBATCH -n 8
#SBATCH -D ./
#SBATCH -J IT64
# SBATCH --export=NONE
#SBATCH -t 1:00:00
# SBATCH --mem-per-cpu=1000000
# SBATCH --array=1-16
#SBATCH --nice=100

module load intel/2013.0.028
#source /c1/apps/intel-cluster-studio/2013.0.028/composer_xe_2013/bin/compilervars.sh intel64
#source /c1/apps/intel-cluster-studio/2013.0.028/impi/4.1.0.024/intel64/bin/mpivars.sh
module load openmpi/intel/64/1.7.4

mpirun ./flash4

