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

## Flag for incrementing output by dimension, or by number of runs
inc_out = False 

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

if inc_out:
    exp_pows = range(7, 45, 1)
    large_ns = range(20000, 75001, 5000)
    ns = []
    for i in range(len(exp_pows)):
        ns.append(int(round(10 + math.pow(1.25, exp_pows[i]), 0)))
    for i in range(len(large_ns)):
        ns.append(large_ns[i])
else:
    ns = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 500, 1000, 1500, 2000, 2500]

num_ns = len(ns)
mcs = 5
n_samp = 10000
pcs = [3,4,5,6]

for j in range(len(pcs)):
    fit_times = np.zeros((mcs, num_ns))
    pred_times = np.zeros((mcs, num_ns))
    iter_pc = pcs[j]
    of_fit_name = 'sepia_fit_times_' + str(iter_pc)
    of_fit_name = of_fit_name + '_dim' if inc_out else of_fit_name + '_ns'
    of_fit_name = of_fit_name + '.csv'
    of_pred_name = 'sepia_pred_times_' + str(iter_pc)
    of_pred_name = of_pred_name + '_dim' if inc_out else of_pred_name + '_ns'
    of_pred_name = of_pred_name + '.csv'
    for i in range(num_ns):
        n = ns[i]
        if inc_out:
            if i >= len(exp_pows):
                print("n = ", str(n))
            else:
                print("n = 1.25^" + str(exp_pows[i]) + " = " + str(n))
        else:
            print("n = ", str(n))
        for k in range(mcs):
            sim_num = sample(range(66), 1)
            non_sim_nums = np.delete(range(66), sim_num)
            if not inc_out:
                non_sim_nums = non_sim_nums[sample(range(len(non_sim_nums)), n)] if n <= len(non_sim_nums) else non_sim_nums[choices(range(len(non_sim_nums)), k=n)]
            x_write = x[non_sim_nums,:]
            if inc_out:
                inds = sample(range(nresponses), n) if n <= nresponses else choices(range(nresponses), k=n)
            else:
                inds = sample(range(16200), 10000)
            y_write = y[non_sim_nums,]
            y_write = y_write[:,inds]
            yx_write = yx[inds,:]
            if not inc_out:
                for l in range(y_write.shape[0]):
                    sd_k = np.std(y_write[l,:])
                    y_write[l,:] += np.random.normal(0, sd_l, size=y_write.shape[1])

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
            data.create_K_basis(n_pc=iter_pc)
            model = SepiaModel(data)
            model.verbose = False
            ## fit model
            model.tune_step_sizes(50, 20, verbose=False)
            model.do_mcmc(n_samp)
            pred_samples = model.get_samples(nburn=int(.1*n_samp),effectivesamples=True)
            toc = time.time()
            fit_times[k,i] = toc - tic
            print("Finished fitting in Monte Carlo iteration " + str(k+1) + "/5 of n=" + str(n) + "with nps=" + str(iter_pc) + " in " + str(toc-tic) " seconds.")
            tic = time.time()
            xx = x[sim_num, :]
            n_pred=xx.shape[0]
            ## Prediction
            pred = SepiaEmulatorPrediction(x_pred=xx, samples=pred_samples, model=model,
               storeMuSigma=True)
            toc = time.time()
            pred_times[k,i] = toc - tic
            print("Finished predicting in Monte Carlo iteration " + str(k+1) + "/5 of n=" + str(n) + "with nps=" + str(iter_pc) + " in " + str(toc-tic) " seconds.")
            ## Delete old files
            os.remove('x_iter.csv')
            os.remove('y_iter.csv')
            os.remove('yx_iter.csv')

    with open(of_fit_name, 'w', newline='') as file:
        # Create a CSV writer object
        writer = csv.writer(file)
        # Write each row of data to the CSV file
        writer.writerows(fit_times)

    with open(of_pred_name, 'w', newline='') as file:
        # Create a CSV writer object
        writer = csv.writer(file)
        # Write each row of data to the CSV file
        writer.writerows(pred_times)
