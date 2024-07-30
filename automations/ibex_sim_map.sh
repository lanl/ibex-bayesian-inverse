#!/bin/bash

#SBATCH -N 32 
#SBATCH --ntasks-per-node=30
#SBATCH -t 8:00:00
#SBATCH -p normal_q
#SBATCH -A ascclass

module reset
module load parallel
module load R/4.2.2-foss-2022b

export OMP_NUM_THREADS=1

if ! [[ $1 =~ ^20[0-2][0-9]A$ ]];
then
  echo "Argument 1 must be a year between 2009-2022"
  echo "Indicates what year the real data should come from."
  echo "Usage: ibex_sim_map.sh [year] [numlocs]"
  exit 1
fi

if ! [[ $2 =~ ^[0-9]+$ ]];
then
  echo "Argument 2 must be an integer"
  echo "Indicates number of unique locations in year file"
  echo "Usage: ibex_sim_map.sh [year] [numlocs]"
  exit 1
fi

ncopies=$2
nparallel=1
nf="nodes_$(( RANDOM % 9000 + 1000)).txt"

# Write file for list of nodes
scontrol show hostname $SLURM_NODELIST > $nf

seq 1 $ncopies | parallel --slf $nf -j$nparallel --wd $PWD --env OMP_NUM_THREADS "module reset; module load R/4.2.2-foss-2022b; R CMD BATCH \"--args row={} map=$1\" ibex_sim_map.R ibex_sim_map_{}.Rout"

# Collect all the results
R CMD BATCH "--args map=$1" ibex_sim_map_collect.R

# Remove node file
if [[ -f $nf ]];
then
  rm $nf
fi

