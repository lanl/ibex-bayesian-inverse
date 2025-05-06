# python3 -m venv ~/py_envs
# source ~/py_envs/bin/activate

import csv
import math
import numpy as np
import os
from random import sample, choices
import sys
import time
from sepia.SepiaModel import SepiaModel
from sepia.SepiaData import SepiaData
from sepia.SepiaPredict import SepiaEmulatorPrediction

exp_pows = range(7, 45, 1)
large_ns = range(20000, 75001, 5000)
ns = []
for i in range(len(exp_pows)):
    ns.append(int(round(10 + math.pow(1.25, exp_pows[i]), 0)))

for i in range(len(large_ns)):
    ns.append(large_ns[i])

num_ns = len(ns)
mcs = 10
fit_times = np.zeros((mcs, num_ns))
pred_times = np.zeros((mcs, num_ns))
n_samp = 10000

for i in range(num_ns):
    n = ns[i]
    if i >= len(exp_pows):
        print("n = ", str(n))
    else:
        print("n = 1.25^" + str(exp_pows[i]) + " = " + str(n))
    xtest_fn = '../automations/xtest_n' + str(n) + '.csv'
    xtest = np.genfromtxt(xtest_fn, delimiter=',', skip_header=1)
    xtest = np.reshape(xtest[0,2:6], (1,4))

    for j in range(mcs):
        xtrain_fn = '../automations/xtrain_n' + str(n) + '_mc' + str(j+1) + '.csv'
        ytrain_fn = '../automations/ytrain_n' + str(n) + '_mc' + str(j+1) + '.csv'
        tic = time.time()
        xtrain = np.genfromtxt(xtrain_fn, delimiter=',', skip_header=1)
        ytrain = np.genfromtxt(ytrain_fn, delimiter=',', skip_header=1)
        yx = xtrain[range(n),0:2]
        x = xtrain[range(0, xtrain.shape[0], n),2:xtrain.shape[1]+1]
        y = np.zeros((x.shape[0], n))
        beg = 0
        end = n
        for k in range(x.shape[0]):
            y[k,] = ytrain[range(beg, end)]
            beg += n
            end += n
        data = SepiaData(x_sim=x, y_sim=y, y_ind_sim=yx)
        data.transform_xt()
        data.standardize_y()
        ## create basis
        data.create_K_basis(n_pc=0.99)
        model = SepiaModel(data)
        model.verbose = False
        ## fit model
        model.tune_step_sizes(50, 20, verbose=False)
        model.do_mcmc(n_samp)
        pred_samples = model.get_samples(nburn=int(.1*n_samp),effectivesamples=True)
        toc = time.time()
        fit_times[j,i] = toc - tic
        tic = time.time()
        n_pred=xtest.shape[0]
        ## Prediction
        pred = SepiaEmulatorPrediction(x_pred=xtest, samples=pred_samples, model=model,
           storeMuSigma=True)
        toc = time.time()
        pred_times[j,i] = toc - tic
        print("Finished Monte Carlo iteration " + str(j+1) + "/" + str(mcs) + " of n=" + str(n))

with open('sepia_fit_times.csv', 'w', newline='') as file:
    # Create a CSV writer object
    writer = csv.writer(file)
    # Write each row of data to the CSV file
    writer.writerows(fit_times)

with open('sepia_pred_times.csv', 'w', newline='') as file:
    # Create a CSV writer object
    writer = csv.writer(file)
    # Write each row of data to the CSV file
    writer.writerows(pred_times)
