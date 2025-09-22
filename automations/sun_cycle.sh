#!/bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=64
#SBATCH -t 4:00:00
#SBATCH -p normal_q
#SBATCH -A ibex

module reset
module load parallel
module load R/4.4.2-gfbf-2024a

export OMP_NUM_THREADS=2

ncopies=64
nparallel=32

seq 1 $ncopies | parallel --slf -j$nparallel --wd $PWD --env OMP_NUM_THREADS "module reset; module load R/4.4.2-gfbf-2024a; R CMD BATCH \"--args index={}\" sun_cycle_calib.R sun_cycle_{}.Rout"
