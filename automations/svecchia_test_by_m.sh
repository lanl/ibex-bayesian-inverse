#!/bin/bash

R CMD BATCH "--args seed=$((RANDOM)) start=1" svecchia_test_by_m.R svecchia\_test\_by\_m.Rout
