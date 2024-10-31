#!/bin/bash
#SBATCH -n 10
#SBATCH -p general
#SBATCH -o output.txt
#SBATCH -e error.txt

module load gnu12 mpich openblas 

mpirun -np 10 ./nrlmol_exe > print.log
