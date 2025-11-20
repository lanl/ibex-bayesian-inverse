#!/bin/bash

#SBATCH -N 10
#SBATCH --ntasks-per-node=8
#SBATCH -t 72:00:00
#SBATCH -p normal_q
#SBATCH -A ibex

module reset
module load parallel
module load R/4.4.2-gfbf-2024a

export OMP_NUM_THREADS=2

map=$1 # map year on which to base simulated field data
ncopies=10
nf="nodes_$(( RANDOM % 9000 + 1000)).txt"
scontrol show hostname $SLURM_NODELIST > $nf

seq 1 $ncopies | parallel --slf $nf -j$ncopies --wd $PWD --env OMP_NUM_THREADS "module reset; module load R/4.4.2-gfbf-2024a; R CMD BATCH \"--args index={} seed={} map=$map\" scale_disc_test.R scale_disc_test_{}.Rout"

# # Collect all the results
R CMD BATCH scale_disc_test_collect.R

# Remove node file
if [[ -f $nf ]];
then
  rm $nf
fi
