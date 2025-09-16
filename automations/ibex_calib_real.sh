#!/bin/bash

#SBATCH -N 1 
#SBATCH --ntasks-per-node=8
#SBATCH -t 96:00:00
#SBATCH -p normal_q
#SBATCH -A ascclass

module reset
module load R/4.1.0-foss-2021a

export OMP_NUM_THREADS=2

if ! [[ $1 =~ ^[0-9]+$ ]];
then
  echo "Argument 1 must be integer. Indicates number of Monte Carlo repetitions."
  echo "Usage: ibex_calib_real.sh [MC reps] [year]"
  exit 1
fi

if ! [[ $2 =~ ^20[0-2][0-9]A$ ]];
then
  echo "Argument 2 must be a year between 2009-2022."
  echo "Indicates what year the real data should come from."
  echo "Usage: ibex_calib_real.sh [MC reps] [year]"
  exit 1
fi

echo "Started script"
for (( i=1; i<=$1; i++ ))
do
  echo "Beginning MC iteration $i."
  R CMD BATCH "--args -r -v --procs=14 --tls=0 --fyear=$2" calib.R
  echo "Finished MC iteration $i."
done
echo "Finished calibration script"
