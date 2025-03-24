#!/bin/bash

export OMP_NUM_THREADS=1
R CMD BATCH "--args seed=$((RANDOM))" surrogate_time_test.R surrogate\_time\_test.Rout
