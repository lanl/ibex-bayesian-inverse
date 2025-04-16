# python3 -m venv ~/py_envs
# source ~/py_envs/bin/activate

import csv
import math
import numpy as np
import random
import time
from scipy.stats import norm
from sepia.SepiaModel import SepiaModel
from sepia.SepiaData import SepiaData
from sepia.SepiaPredict import SepiaEmulatorPrediction

seed = random.randint(1000000, 9999999)
fn_rmse = 'sepia_basis_rmses_' + str(seed) + '.csv'
fn_crps = 'sepia_basis_crps_' + str(seed) + '.csv'

model_data = np.genfromtxt('../data/sims.csv', delimiter=',', names=True)
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

n_samp = 10000

nbases = []
nbases_small = range(1, 11, 1)
nbases_large = range(10, 65, 2)
for i in range(len(nbases_small)):
    nbases.append(nbases_small[i])

for i in range(len(nbases_large)):
    nbases.append(nbases_large[i])

metrics = np.zeros((len(nbases), nruns, 2))

for i in range(len(nbases)):
    nbase = nbases[i]
    for j in range(nruns):
        yy_true = y[[j],:][0]
        inds = np.delete(range(nruns), j)
        data = SepiaData(x_sim=x[inds,], y_sim=y[inds,], y_ind_sim=yx)
        data.transform_xt()
        data.standardize_y()
        data.create_K_basis(n_pc=nbase, fm=False)
        model = SepiaModel(data)
        model.verbose = False
        model.tune_step_sizes(50, 20, verbose=False)
        print("\n")
        model.do_mcmc(n_samp)
        xx = x[[j], :]
        n_pred=xx.shape[0]
        pred_samples = model.get_samples(nburn=int(.1*n_samp),effectivesamples=True)
        pred = SepiaEmulatorPrediction(x_pred=xx, samples=pred_samples, model=model,
           storeMuSigma=True)
        predy = pred.get_y()
        meany = np.mean(predy,0)[0]
        s2s = np.var(predy,0)[0]
        sumsq = 0
        for k in range(meany.shape[0]):
            sumsq += pow(meany[k] - yy_true[k], 2)
        metrics[i,j,0] = math.pow(sumsq/meany.shape[0], 0.5)
        scores = np.zeros(yx.shape[0])
        for k in range(meany.shape[0]):
            sigma = pow(s2s[k], 0.5)
            if (sigma == 0):
                scores[k] = np.nan
                continue
            z = (yy_true[k] - meany[k])/sigma
            scores[k] = sigma*(-1/(math.pow(math.pi, 0.5)) + 2*norm.pdf(z) + z*(2*norm.cdf(z)-1))
        metrics[i,j,1] = np.nanmean(scores)
        print("\n")
        print("Finished holdout iteration " + str(j) + "\n")
        with open(fn_rmse, 'w', newline='') as file:
            # Create a CSV writer object
            writer = csv.writer(file)
            # Write each row of data to the CSV file
            writer.writerows(metrics[:,:,0])
        with open(fn_crps, 'w', newline='') as file:
            # Create a CSV writer object
            writer = csv.writer(file)
            # Write each row of data to the CSV file
            writer.writerows(metrics[:,:,1])

    print("Finished iteration with " + str(nbase) + " basis functions\n")

with open(fn_rmse, 'w', newline='') as file:
    # Create a CSV writer object
    writer = csv.writer(file)
    # Write each row of data to the CSV file
    writer.writerows(metrics[:,:,0])

with open(fn_crps, 'w', newline='') as file:
    # Create a CSV writer object
    writer = csv.writer(file)
    # Write each row of data to the CSV file
    writer.writerows(metrics[:,:,1])
