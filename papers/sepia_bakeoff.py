# python3 -m venv ~/py_envs
# source ~/py_envs/bin/activate

import csv
import math
import numpy as np
import time
from sepia.SepiaModel import SepiaModel
from sepia.SepiaData import SepiaData
from sepia.SepiaPredict import SepiaEmulatorPrediction

x = np.genfromtxt('calib_params_ibex.csv', delimiter=',', skip_header=1) ## pmfp, ratio
y = np.transpose(np.genfromtxt('ibex_responses.csv', delimiter=',', skip_header=1)) ## 66 vectors of 16200
yx = np.genfromtxt('xyz_grid.csv', delimiter=',', skip_header=1) ## xyz grid

n_samp = 10000
preds = np.zeros((len(x), yx.shape[0]))
s2s = np.zeros((len(x), yx.shape[0]))
rmses = np.zeros(len(x))
crps = np.zeros(len(x))

for i in range(len(x)):
    yy_true = y[[i],:][0]
    inds = np.delete(range(len(x)), i)
    data = SepiaData(x_sim = x[inds,], y_sim = y[inds,], y_ind_sim = yx)
    data.transform_xt()
    data.standardize_y()
    data.create_K_basis(n_pc=0.99)
    model = SepiaModel(data)
    model.verbose = False
    model.tune_step_sizes(50, 20, verbose=False)
    print("\n")
    model.do_mcmc(n_samp)
    xx = x[[i], :]
    n_pred=xx.shape[0]
    pred_samples = model.get_samples(nburn=int(.1*n_samp),effectivesamples=True)
    pred = SepiaEmulatorPrediction(x_pred=xx, samples=pred_samples, model=model,
       storeMuSigma=True)
    predy = pred.get_y()
    meany = np.mean(predy,0)[0]
    preds[i,] = np.mean(predy,0)
    s2s[i,] = np.var(predy,0)
    sumsq = 0
    for j in range(meany.shape[0]):
        sumsq += pow(meany[j] - yy_true[j], 2)
    rmses[i] = math.pow(sumsq/meany.shape[0], 0.5)
    print("\n")
    print("Finished iteration " + str(i) + "\n")

with open('sepia_preds.csv', 'w', newline='') as file:
    # Create a CSV writer object
    writer = csv.writer(file)  
    # Write each row of data to the CSV file
    writer.writerows(preds)

with open('sepia_rmses.csv', 'w', newline='') as file:
    # Create a CSV writer object
    writer = csv.writer(file)  
    # Write each row of data to the CSV file
    writer.writerow(rmses)
