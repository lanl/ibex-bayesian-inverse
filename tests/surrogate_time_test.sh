#!/bin/bash

###############################################################################
#### Testing speed of surrogate fits and predictions on IBEX simulator data.
#### More details in surrogate_time_test.R
###############################################################################

if ! [[ $1 =~ ^[01]$ ]];
then
  echo "Argument 1 must be either 0 or 1."
  echo "Indicates if the test should varying dimesion of the response (1) or number of simulator runs (0)."
  echo "Usage: surrogate_time_test.sh [flag for varying response dimension]"
  exit 1
fi

inc_out=$1 # flag indicating if test should be on increasing dimension of response

R CMD BATCH "--args seed=$((RANDOM)) inc_out=$inc_out" surrogate_time_test.R surrogate\_time\_test.Rout
