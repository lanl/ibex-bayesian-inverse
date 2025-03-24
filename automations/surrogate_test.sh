#!/bin/bash

export OMP_NUM_THREADS=1
R CMD BATCH "--args seed=$((RANDOM)) method='all' start=1" surrogate_test.R surrogate\_test.Rout
