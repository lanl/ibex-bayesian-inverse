# python3 -m venv ~/py_envs
# source ~/py_envs/bin/activate

###############################################################################
###############################################################################
## Hold one out test for SEPIA on IBEX simulator output. Used to generate
## Figure 6. SEPIA is run once for each unique combination of simulator model
## parameters. Predictions are made at a held out combination of model
## parameters. RMSE and CRPS are calculated. Varying number of bases
## (principal components) are used
## DATA NEEDED: sims.csv
###############################################################################
###############################################################################

import csv
import math
import numpy as np
import time
from scipy.stats import norm
from scipy.spatial.distance import pdist
from sepia.SepiaModel import SepiaModel
from sepia.SepiaData import SepiaData
from sepia.SepiaPredict import SepiaEmulatorPrediction

def energy_score(draws, observed):
    draws = np.asarray(draws)
    observed = np.asarray(observed)
    n = draws.shape[0]
    dist_to_obs = np.linalg.norm(draws - observed,  axis=1)
    pairwise_distance_sum = pdist(draws).sum()
    return(dist_to_obs.mean()-pairwise_distance_sum / n**2)

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

pcs = [3,4,5,6]
n_samp = 10000

for j in range(len(pcs)):
    metrics = np.zeros((nruns, 4))
    resids = np.zeros((16200, nruns))
    iter_pc = pcs[j]
    met_ofname = '../papers/sepia_metrics_' + str(iter_pc) + '.csv'
    resids_ofname = '../papers/sepia_resids_' + str(iter_pc) + '.csv'
    for i in range(nruns):
        yy_true = y[[i],:][0]
        inds = np.delete(range(nruns), i)
        data = SepiaData(x_sim=x[inds,], y_sim=y[inds,], y_ind_sim=yx)
        data.transform_xt()
        data.standardize_y()
        data.create_K_basis(n_pc=iter_pc)
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
        s2s = np.var(predy,0)[0]
        sumsq = 0
        for j in range(meany.shape[0]):
            sumsq += pow(meany[j] - yy_true[j], 2)
        metrics[i,0] = math.pow(sumsq/meany.shape[0], 0.5)
        scores = np.zeros(yx.shape[0])
        for j in range(meany.shape[0]):
            sigma = pow(s2s[j], 0.5)
            if (sigma == 0):
                scores[j] = np.nan
                continue
            z = (yy_true[j] - meany[j])/sigma
            scores[j] = sigma*(-1/(math.pow(math.pi, 0.5)) + 2*norm.pdf(z) + z*(2*norm.cdf(z)-1))
        metrics[i,1] = np.nanmean(scores)
        metrics[i,2] = energy_score(draws=predy[8000:9000,0,:], observed=yy_true)
        metrics[i,3] = np.max(np.abs(meany - yy_true))
        resids[:,i] = meany - yy_true
        print("\n")
        print("Finished iteration " + str(i) + "with nps=" + str(iter_pc) + "\n")
    with open(met_ofname, 'w', newline='') as file:
        # Create a CSV writer object
        writer = csv.writer(file)
        # Write each row of data to the CSV file
        writer.writerows(metrics)
    with open(resids_ofname, 'w', newline='') as file:
        # Create a CSV writer object
        writer = csv.writer(file)
        # Write each row of data to the CSV file
        writer.writerows(resids)
