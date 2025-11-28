#!/bin/bash

###############################################################################
#### Shell script for testing recovery of a multiplicative scale discrepancy.
#### More detail in scale_disc_test.R
###############################################################################

if ! [[ $1 =~ ^[0-9]+$ ]];
then
  echo "Argument 1 must be an integer. Indicates the year of satellite data to use."
  echo "Usage: scale_disc_test.sh [year] [nprocs]"
  exit 1
fi

if ! [[ $2 =~ ^[0-9]+$ ]];
then
  echo "Argument 2 must be an integer. Indicates the number of parallel processes to run."
  echo "Usage: scale_disc_test.sh [year] [nprocs]"
  exit 1
fi

map=$1 # map year on which to base simulated field data
ncopies=$2 # number of processes to run in parallel

seq 1 10 | parallel -j $ncopies --wd $PWD "R CMD BATCH \"--args index={} seed={} map=$map\" scale_disc_test.R scale_disc_test_{}.Rout"

# Collect all the results
R CMD BATCH scale_disc_test_collect.R
