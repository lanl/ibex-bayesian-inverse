#!/bin/bash

#SBATCH -N 16
#SBATCH --ntasks-per-node=20
#SBATCH -t 1:00:00
#SBATCH -p normal_q
#SBATCH -A ascclass

module reset
module load parallel
module load R/4.2.2-foss-2022b

export OMP_NUM_THREADS=1

ncopies=14154
nparallel=1
nf="nodes_$(( RANDOM % 9000 + 1000)).txt"

# Write file for list of nodes
scontrol show hostname $SLURM_NODELIST > $nf

seq 1 $ncopies | parallel --slf $nf -j$nparallel --wd $PWD --env OMP_NUM_THREADS "module reset; module load R/4.2.2-foss-2022b; R CMD BATCH \"--args row={}\" ibex_sim_map.R ibex_sim_map_{}.Rout"

# Collect all the results
R CMD BATCH ibex_sim_map_collect.R
