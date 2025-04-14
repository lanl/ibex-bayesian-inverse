# python3 -m venv ~/py_envs
# source ~/py_envs/bin/activate
# ensure fullMatrices=True is set

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

## Read in the data
model_data = np.genfromtxt('../data/sims.csv', delimiter=',', names=True)

## Format the data
x = model_data[['parallel_mean_free_path', 'ratio']]
x = np.unique(x, axis=0)
x = np.sort(x, order=['parallel_mean_free_path', 'ratio'])
x = np.array(x.tolist())
nruns = x.shape[0]

y_np = np.sort(model_data, order=['parallel_mean_free_path', 'ratio', 'lat', 'lon'])
y_np = y_np['blurred_ena_rate']
nresponses = int(len(y_np) / nruns)
y = np.zeros((nruns, nresponses))
start = 0
end = nresponses
for i in range(66):
    y[i,:] = y_np[start:end]
    start += nresponses
    end += nresponses

latlon = np.sort(np.unique(model_data[['lat', 'lon']], axis=0), order=['lat', 'lon'])
yx = np.zeros((nresponses, 3))
for i in range(len(yx)):
    yx[i,0] = math.cos(math.pi*latlon['lon'][i]/180)*math.cos(math.pi*latlon['lat'][i]/180)
    yx[i,1] = math.sin(math.pi*latlon['lon'][i]/180)*math.cos(math.pi*latlon['lat'][i]/180)
    yx[i,2] = math.sin(math.pi*latlon['lat'][i]/180)

large_n = int(sys.argv[1])
exp_pows = range(7, 45, 1)
ns = range(20000, 75001, 5000)
num_ns = len(ns) if large_n else len(exp_pows)
mcs = 5
fit_times = np.zeros((mcs, len(exp_pows)))
pred_times = np.zeros((mcs, len(exp_pows)))
n_samp = 10000

for i in range(num_ns):
    ## select training data size
    if large_n:
        n = ns[i]
        print("n = ", str(n))
    else:
        n = int(round(10 + math.pow(1.25, exp_pows[i]), 0))
        print("n = 1.25^" + str(exp_pows[i]) + " = " + str(n))

    for j in range(mcs):

        sim_num = sample(range(66), 1)
        non_sim_nums = np.delete(range(66), sim_num)
        x_write = x[nonsim_nums,]
        inds = sample(range(16200), n) if n <= 16200 else choices(range(16200), k=n)
        y_write = y[non_sim_nums,]
        y_write = y_write[:,inds]
        yx_write = yx[inds,]
        with open('x_iter.csv', 'w', newline='') as file:
            # Create a CSV writer object
            writer = csv.writer(file)
            # Write each row of data to the CSV file
            writer.writerows(x_write)
        with open('y_iter.csv', 'w', newline='') as file:
            # Create a CSV writer object
            writer = csv.writer(file)
            # Write each row of data to the CSV file
            writer.writerows(y_write)
        with open('yx_iter.csv', 'w', newline='') as file:
            # Create a CSV writer object
            writer = csv.writer(file)
            # Write each row of data to the CSV file
            writer.writerows(yx_write)

        tic = time.time()
        ## read in data
        x_iter = np.genfromtxt('x_iter.csv', delimiter=',')
        y_iter = np.genfromtxt('y_iter.csv', delimiter=',')
        yx_iter = np.genfromtxt('yx_iter.csv', delimiter=',')

        ## create object
        data = SepiaData(x_sim=x_iter, y_sim=y_iter, y_ind_sim=yx_iter)
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
        xx = x[sim_num, :]
        n_pred=xx.shape[0]
        ## Prediction
        pred = SepiaEmulatorPrediction(x_pred=xx, samples=pred_samples, model=model,
           storeMuSigma=True)
        toc = time.time()
        pred_times[j,i] = toc - tic
        print("Finished Monte Carlo iteration " + str(j+1) + "/5 of n=" + str(n))

        ## Delete old files
        os.remove('x_iter.csv')
        os.remove('y_iter.csv')
        os.remove('yx_iter.csv')

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
