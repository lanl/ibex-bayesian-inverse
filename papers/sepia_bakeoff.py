import csv
import numpy as np
import time
from sepia.SepiaModel import SepiaModel
from sepia.SepiaData import SepiaData
from sepia.SepiaPredict import SepiaEmulatorPrediction

x = np.genfromtxt('calib_params_ibex.csv', delimiter=',', skip_header=1) ## pmfp, ratio
y = np.transpose(np.genfromtxt('ibex_responses.csv', delimiter=',', skip_header=1)) ## 66 vectors of 16200
yx = np.genfromtxt('xyz_grid.csv', delimiter=',', skip_header=1) ## xyz grid

n_samp = 10000
times = np.zeros((len(x), 4))
preds = np.zeros((len(x), yx.shape[0]))
stddevs = np.zeros((len(x), yx.shape[0]))
for i in range(len(x)):
    yy_true = y[[i],:]
    inds = np.delete(range(len(x)), i)
    data = SepiaData(x_sim = x[inds,], y_sim = y[inds,], y_ind_sim = yx)
    data.transform_xt()
    data.standardize_y()
    data.create_K_basis(n_pc=0.99)
    model = SepiaModel(data)
    model.verbose = False
    tic = time.time()
    model.tune_step_sizes(50, 20, verbose=False)
    print("\n")
    model.do_mcmc(n_samp)
    toc = time.time()
    times[i,0] = toc - tic
    xx = x[[i], :]
    n_pred=xx.shape[0]
    pred_samples = model.get_samples(nburn=int(.1*n_samp),effectivesamples=True)
    tic = time.time()
    pred = SepiaEmulatorPrediction(x_pred=xx, samples=pred_samples, model=model,
       storeMuSigma=True)
    toc = time.time()
    times[i,1] = toc - tic
    predy = pred.get_y()
    preds[i,] = np.mean(predy,0)
    stddevs[i,] = np.std(predy,0)
    print("\n")
    print("Finished iteration " + str(i) + "\n")

with open('sepia_times.csv', 'w', newline='') as file:
    # Create a CSV writer object
    writer = csv.writer(file)  
    # Write each row of data to the CSV file
    writer.writerows(times)

with open('sepia_preds.csv', 'w', newline='') as file:
    # Create a CSV writer object
    writer = csv.writer(file)  
    # Write each row of data to the CSV file
    writer.writerows(preds)
