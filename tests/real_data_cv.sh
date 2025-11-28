#!/bin/bash

###############################################################################
#### Bash script for launch cross-validation code for solving a Bayesian
#### inverse problem. Details in real_data_cv.R
###############################################################################

if ! [[ $1 =~ ^[0-9]+$ ]];
then
  echo "Argument 1 must be an integer. Indicates number of MCMC iterations."
  echo "Usage: real_data_cv.sh [number of mcmcs] [number of folds] [nprocs]"
  exit 1
fi

if ! [[ $2 =~ ^[0-9]+$ ]];
then
  echo "Argument 2 must be an integer. Indicates number of cross-validation folds."
  echo "Usage: real_data_cv.sh [number of mcmcs] [number of folds] [nprocs]"
  exit 1
fi

if ! [[ $3 =~ ^[0-9]+$ ]];
then
  echo "Argument 3 must be an integer. Indicates the number of parallel processes to run."
  echo "Usage: real_data_cv.sh [number of mcmcs] [number of folds] [nprocs]"
  exit 1
fi

nmcmc=$1 # number of mcmcs to run
ncv=$2 # number of folds in cross validation
nparallel=$3 # number of parallel instances

seed=$((RANDOM))

seq 1 $ncv | parallel -j $nparallel --wd $PWD "R CMD BATCH \"--args seed=$seed fold_seed={} ncvs=$ncv fid={} nmcmcs=$nmcmc\" real_data_cv.R real_data_cv_{}_$seed.Rout"

# Collect all the results
R CMD BATCH "--args seed=$seed" real_data_cv_collect.R
