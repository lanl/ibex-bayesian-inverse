#!/bin/bash

#SBATCH -N 1 
#SBATCH --ntasks-per-node=66
#SBATCH -t 96:00:00
#SBATCH -p normal_q
#SBATCH -A ibex

module reset
module load R/4.1.0-foss-2021a

export OMP_NUM_THREADS=2

echo "Started script"
R CMD BATCH "--args -v --infile=inputs.txt" ibex_bayes_test.R
echo "Finished calibration script"
