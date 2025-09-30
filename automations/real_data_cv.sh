#!/bin/bash

#SBATCH -N 5
#SBATCH --ntasks-per-node=20
#SBATCH -t 24:00:00
#SBATCH -p normal_q
#SBATCH -A ibex

module reset
module load parallel
module load R/4.4.2-gfbf-2024a

export OMP_NUM_THREADS=2

ncopies=10
nparallel=10
seed=$((RANDOM))
nf="nodes_$(( RANDOM % 9000 + 1000)).txt"
scontrol show hostname $SLURM_NODELIST > $nf

seq 1 $ncopies | parallel --slf $nf -j$nparallel --wd $PWD --env OMP_NUM_THREADS "module reset; module load R/4.4.2-gfbf-2024a; R CMD BATCH \"--args seed=$seed fold_seed=$((RANDOM)) ncvs=10 fid={} nmcmcs=10000\" real_data_cv.R real_data_cv_{}_$seed.Rout"

# Collect all the results
R CMD BATCH "--args seed=$seed" real_data_cv_collect.R

# Remove node file
if [[ -f $nf ]];
then
  rm $nf
fi
