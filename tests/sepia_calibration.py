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

import os
import csv
import math
import numpy as np
import time

os.chdir("../../SEPIA/sepia")

from scipy.interpolate import interp2d
from scipy.interpolate import RegularGridInterpolator
from scipy.stats import norm
from scipy.spatial.distance import pdist
from sepia.SepiaModel import SepiaModel
from sepia.SepiaData import SepiaData
from sepia.SepiaPredict import SepiaEmulatorPrediction

os.chdir("../../ibex-bayesian-inverse/tests")

## SIMULATION DATA
model_data = np.genfromtxt('../data/sims.csv', delimiter=',', names=True)

nlon = len(np.unique(model_data['lon']))
nlat = len(np.unique(model_data['lat']))
nresponses = nlon*nlat
nruns = np.unique(model_data[['parallel_mean_free_path', 'ratio']]).shape[0]
param_names = ['ratio', 'pmfp']

model_dat_sort = np.sort(model_data, order=['ratio', 'parallel_mean_free_path', 'lon', 'lat'])
logy_np = np.log(model_dat_sort['blurred_ena_rate']+0.65)
logy_sim = np.zeros((nruns, nresponses))
start = 0
end = nresponses
for i in range(66):
    logy_sim[i,:] = logy_np[start:end]
    start += nresponses
    end += nresponses

logy_sim_ind_lon_lat = np.zeros((nresponses, 2))
lons = np.sort(np.unique(model_dat_sort['lon']))
lats = np.sort(np.unique(model_dat_sort['lat']))
logy_sim_ind_lon_lat[:,0] = np.repeat(lons, nlat)
logy_sim_ind_lon_lat[:,1] = np.tile(lats, nlon)
t_sim=np.stack(np.array(np.sort(np.unique(model_dat_sort[['ratio', 'parallel_mean_free_path']]), order=['ratio', 'parallel_mean_free_path']), dtype=object))

## read in real data
real_data = np.genfromtxt('../data/synth_sat_data.csv', delimiter=',', dtype=None, encoding="utf-8", names=True)
unique_combos = np.unique(real_data[['ratio', 'parallel_mean_free_path']])

results = {}
## for each unique setting, run calibration.
for i in range(len(unique_combos)):
    ## OBSERVED DATA
    iter_ratio = unique_combos[i]['ratio']
    iter_pmfp = unique_combos[i]['parallel_mean_free_path']
    print("Started iteration " + str(i) + " with ratio=" + str(iter_ratio) + " and pmfp=" + str(iter_pmfp))
    iter_data = real_data[(real_data['ratio']==iter_ratio) & (real_data['parallel_mean_free_path']==iter_pmfp)]
    x_obs = iter_data[['lon', 'lat']]
    x_obs = np.array(x_obs.tolist(), dtype=object)
    y_obs = iter_data['sim_counts'] / iter_data['time']
    y_obs[np.isnan(y_obs)] = 0
    y_obs = y_obs - iter_data['background']
    logy_obs = np.log(y_obs + 0.65)
    logy_obs = [logy_obs]
    logy_ind_obs = [x_obs]

    ## SEPIA DATA
    data = SepiaData(x_sim=None, t_sim=t_sim, y_sim=logy_sim, y_ind_sim=logy_sim_ind_lon_lat,
        x_obs=None, y_obs=logy_obs, y_ind_obs=logy_ind_obs)
    data.transform_xt()

    ## K BASIS
    logysimmean = np.mean(logy_sim,0)
    logysimsd = np.std(logy_sim)
    logysimStd = (logy_sim - np.tile(logysimmean,nruns).reshape(logy_sim.shape))/logysimsd
    data.sim_data.orig_y_mean = logysimmean
    data.sim_data.orig_y_sd = logysimsd
    data.sim_data.y_std = logysimStd

    U, s, V = np.linalg.svd(logysimStd.T, full_matrices=False)
    numPC = 3
    data.sim_data.K = U[:,0:numPC]*s[0:numPC]/np.sqrt(nruns)
    data.sim_data.K = data.sim_data.K.T

    # obs
    latmat = np.repeat(lats,nlon).reshape((nlon,nlat),order='F')
    lonmat = np.repeat(lons,nlat).reshape((nlon,nlat))
    # compute simulator mean values simdat.ymean interpolated to the data values...
    interp = RegularGridInterpolator((lons, lats), data.sim_data.orig_y_mean.reshape((nlon,nlat),order='F'), method='linear', bounds_error=False, fill_value=None)

    data.obs_data.orig_y_mean = []
    data.obs_data.orig_y_sd = []
    for k in range(1):
        points = np.column_stack([data.obs_data.y_ind[k][:,0], data.obs_data.y_ind[k][:,1]])
        ymk = interp(points)
        data.obs_data.orig_y_mean.append(ymk.flatten())
        data.obs_data.orig_y_sd.append(data.sim_data.orig_y_sd)

    # now compute the centered, scaled observed arrival times yStd
    data.obs_data.y_std = []
    for k in range(1):
        data.obs_data.y_std.append((data.obs_data.y[k] - data.obs_data.orig_y_mean[k])/data.sim_data.orig_y_sd)

    # for now, hack this in - if it used the inbuilt methods it would happen automatically
    tSigy_std=[]
    for i in range(len(data.obs_data.y)):
        tSigy_std.append(np.atleast_2d(np.diag(np.ones(data.obs_data.y[i].shape))))

    data.obs_data.Sigy_std = tSigy_std
    del tSigy_std

    data.obs_data.K = []
    for k in range(1):
        data.obs_data.K.append(np.zeros((data.obs_data.y_std[k].shape[0], numPC)))
        for j in range(numPC):
            f = RegularGridInterpolator((lons, lats), np.reshape(data.sim_data.K[j,:],(nlon, nlat),order='F'), method='linear', bounds_error=False, fill_value=None)
            fpoints = np.column_stack([data.obs_data.y_ind[k][:,0], data.obs_data.y_ind[k][:,1]])
            data.obs_data.K[k][:,j] = f(fpoints)

    for k in range(1):
        data.obs_data.K[k] = data.obs_data.K[k].T

    model = SepiaModel(data)
    model.tune_step_sizes(50, 20, update_vals=True)
    model.do_mcmc(10000)
    samples_dict = model.get_samples()
    samples_dict['theta'][:,0] = samples_dict['theta'][:,0]*(0.1-0.001)+0.001
    samples_dict['theta'][:,1] = samples_dict['theta'][:,1]*(3000-500)+500
    results[(iter_ratio, iter_pmfp)] = samples_dict['theta']

    with open("sepia_sim_calib_results.csv", "w", newline="") as f:
        writer = csv.writer(f)
        # Header — adjust based on array length
        n = len(next(iter(results.values())))
        header = ["ratio", "pmfp", "parameter"] + [f"val_{i}" for i in range(n)]
        writer.writerow(header)
        
        for (ratio, pmfp), arr in results.items():
            for col_idx in range(arr.shape[1]):          # loop over the 2 columns
                row = [ratio, pmfp, param_names[col_idx]] + list(arr[:, col_idx])
                writer.writerow(row)
    print("Finished iteration " + str(i) + " with ratio=" + str(iter_ratio) + " and pmfp=" + str(iter_pmfp))

