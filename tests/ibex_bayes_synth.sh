#!/bin/bash

if [ -z "$1" ];
then
  echo "Error: missing path to input file"
  echo "Usage: ibex_bayes_synth.sh [input file path]"
  exit 1
fi

infile=$1 # input file with simulator paramters

if [ ! -f "$infile" ];
then
  echo "Error: file does not exist"
  echo "Argument 1 must be a valid file path."
  echo "Refer to ibex_bayes_test.R for expected file format."
  exit 1
fi

echo "Started script"
R CMD BATCH "--args -v --if=$infile" ibex_bayes_test.R
echo "Finished calibration script"

# Collect all the results
R CMD BATCH ibex_bayes_test_collect.R
