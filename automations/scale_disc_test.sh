#!/bin/bash

#SBATCH -N 15
#SBATCH --ntasks-per-node=8
#SBATCH -t 72:00:00
#SBATCH -p normal_q
#SBATCH -A ibex

module reset
module load parallel
module load R/4.4.2-gfbf-2024a

export OMP_NUM_THREADS=2

map=$1 # map year on which to base simulated field data
ncv=$2 # number of folds in cross validation
ncopies=$(($ncv*10))
nparallel=15
nf="nodes_$(( RANDOM % 9000 + 1000)).txt"
scontrol show hostname $SLURM_NODELIST > $nf

seq 1 $ncopies | parallel --slf $nf -j$nparallel --wd $PWD --env OMP_NUM_THREADS "module reset; module load R/4.4.2-gfbf-2024a; R CMD BATCH \"--args index={} seed={} ncvs=$ncv map=$map\" scale_disc_test.R scale_disc_test_{}_$seed.Rout"

# # Collect all the results
# R CMD BATCH "--args seed=$seed" real_data_cv_collect.R

# Remove node file
if [[ -f $nf ]];
then
  rm $nf
fi
