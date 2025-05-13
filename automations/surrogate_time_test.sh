#!/bin/bash

R CMD BATCH "--args seed=$((RANDOM)) large_n=1 inc_out=0" surrogate_time_test.R surrogate\_time\_test.Rout
