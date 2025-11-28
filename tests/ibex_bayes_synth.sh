#!/bin/bash

infile=$1 # input file with simulator paramters

echo "Started script"
R CMD BATCH "--args -v --if=$infile" ibex_bayes_test.R
echo "Finished calibration script"

# Collect all the results
R CMD BATCH ibex_bayes_test_collect.R
