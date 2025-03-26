import csv
import math
import numpy as np
from random import sample, choices
import time
from sepia.SepiaModel import SepiaModel
from sepia.SepiaData import SepiaData
from sepia.SepiaPredict import SepiaEmulatorPrediction

x = np.genfromtxt('calib_params_ibex.csv', delimiter=',') ## pmfp, ratio
y = np.transpose(np.genfromtxt('ibex_responses.csv', delimiter=',')) ## 66 vectors of 16200
yx = np.genfromtxt('lat_lon_grid.csv', delimiter=',') ## lat lon grid

large_n = sys.argv[1]
exp_pows = range(7, 45, 1)
ns = range(20000, 75000, 5000)
num_ns = len(ns) if large_n else len(exp_pows)

mcs = 5
fit_times = np.zeros((mcs, len(exp_pows)))
pred_times = np.zeros((mcs, len(exp_pows)))
n_samp = 10000

for i in range(num_ns):
    ## select training data size
    n = ns[i] if large_n else int(round(10 + math.pow(1.25, exp_pows[i]), 0))
    if large_n:
        n = ns[i]
        print("n = ", str(n))
    else:
        n = int(round(10 + math.pow(1.25, exp_pows[i]), 0))
        print("n = 1.25^" + str(exp_pows[i]) + " = " + str(n))
    for j in range(mcs):
        sim_num = sample(range(66), 1)
        non_sim_nums = np.delete(range(66), sim_num)
        inds = sample(range(16200), n) if n <= 16200 else choices(range(16200), k=n)
        y_iter = y[non_sim_nums,]
        y_iter = y_iter[:,inds]
        data = SepiaData(x_sim = x[non_sim_nums,], y_sim = y_iter, y_ind_sim = yx[inds,])
        data.transform_xt()
        data.standardize_y()
        data.create_K_basis(n_pc=0.99)
        model = SepiaModel(data)
        model.verbose = False
        tic = time.time()
        model.tune_step_sizes(50, 20, verbose=False)
        model.do_mcmc(n_samp)
        toc = time.time()
        fit_times[j,i] = toc - tic
        xx = x[sim_num, :]
        n_pred=xx.shape[0]
        pred_samples = model.get_samples(nburn=int(.1*n_samp),effectivesamples=True)
        tic = time.time()
        pred = SepiaEmulatorPrediction(x_pred=xx, samples=pred_samples, model=model,
           storeMuSigma=True)
        toc = time.time()
        pred_times[j,i] = toc - tic
        print("Finished Monte Carlo iteration " + str(j+1) + "/5 of n=" + str(n))

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
