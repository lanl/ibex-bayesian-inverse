#!/bin/bash

###############################################################################
#### Hold one out testing of different surrogate methods. Details in
#### surrogate_test.R
###############################################################################

R CMD BATCH "--args seed=$((RANDOM)) method='all' start=1" surrogate_test.R surrogate\_test.Rout
