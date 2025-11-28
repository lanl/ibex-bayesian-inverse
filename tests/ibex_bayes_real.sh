#!/bin/bash

if ! [[ $1 =~ ^[0-9]+$ ]];
then
  echo "Argument 1 must be integer. Indicates number of Monte Carlo repetitions."
  echo "Usage: ibex_calib_real.sh [MC reps] [year]"
  exit 1
fi

if ! [[ $2 =~ ^(20[0-2][0-9]|all|mod_align)$ ]];
then
  echo "Argument 2 must be a year between 2009-2022, the string 'all', or the string 'mod_align'."
  echo "Indicates where the real data should come from."
  echo "Usage: ibex_calib_real.sh [MC reps] [year]"
  exit 1
fi

echo "Started script"
for (( i=1; i<=$1; i++ ))
do
  echo "Beginning MC iteration $i."
  R CMD BATCH "--args -v -r --y=$2" ibex_bayes_test.R
  echo "Finished MC iteration $i."
done
echo "Finished calibration script"

# Collect all the results
R CMD BATCH ibex_bayes_test_collect.R
