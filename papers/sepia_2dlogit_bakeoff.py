# python3 -m venv ~/py_envs
# source ~/py_envs/bin/activate

import csv
import math
import numpy as np
import time
from scipy.stats import norm
from sepia.SepiaModel import SepiaModel
from sepia.SepiaData import SepiaData
from sepia.SepiaPredict import SepiaEmulatorPrediction

def f(x, mu, nu):
    t1 = mu[0]*np.exp(mu[0]*x[:,0]-7)/(nu[0]+np.exp(mu[0]*x[:,0]-7))
    t2 = mu[1]*np.exp(mu[1]*x[:,1]-3)/(nu[1]+np.exp(mu[1]*x[:,1]-3))
    return t1+t2

x1_range = [x / 89 for x in range(90)]
x2_range = [x / 89 for x in range(90)]
# x1_range = [x / 4 for x in range(5)]
# x2_range = [x / 4 for x in range(5)]
yx = np.array([(y,x) for x in x1_range for y in x2_range])

Utrue = np.reshape(np.array((11, 0.5, 8, 2.3)), (1,4))
Utrueunit = np.reshape(np.array((11, 0.5, 8, 2.3)), (1,4))
Utrueunit[0,0] = (Utrueunit[0,0] - mu_min) / (mu_max - mu_min)
Utrueunit[0,1] = (Utrueunit[0,1] - nu_min) / (nu_max - nu_min)
Utrueunit[0,2] = (Utrueunit[0,2] - mu_min) / (mu_max - mu_min)
Utrueunit[0,3] = (Utrueunit[0,3] - nu_min) / (nu_max - nu_min)

XX = np.append(yx, np.reshape(np.repeat(Utrueunit[0,0], yx.shape[0]), (yx.shape[0], 1)), 1)
XX = np.append(XX, np.reshape(np.repeat(Utrueunit[0,1], XX.shape[0]), (XX.shape[0], 1)), 1)
XX = np.append(XX, np.reshape(np.repeat(Utrueunit[0,2], XX.shape[0]), (XX.shape[0], 1)), 1)
XX = np.append(XX, np.reshape(np.repeat(Utrueunit[0,3], XX.shape[0]), (XX.shape[0], 1)), 1)
YY = f(x=yx, mu=Utrue[0,[0,2]], nu=Utrue[0,[1,3]])

mu_min = 5
mu_max = 15
nu_min = 0
nu_max = 3

mcs = 30
n_samp = 10000
metrics = np.zeros((mcs, 2))

for i in range(mcs):
    uunit = np.genfromtxt('../automations/surrogate_2dlogit_U' + str(i+1) + ".csv", delimiter=',')
    u = np.genfromtxt('../automations/surrogate_2dlogit_U' + str(i+1) + ".csv", delimiter=',')
    u[:,0] = u[:,0] * (mu_max - mu_min) + mu_min
    u[:,1] = u[:,1] * (nu_max - nu_min) + nu_min
    u[:,2] = u[:,2] * (mu_max - mu_min) + mu_min
    u[:,3] = u[:,3] * (nu_max - nu_min) + nu_min
    y = np.zeros((u.shape[0], yx.shape[0],))
    for j in range(u.shape[0]):
        y[j,] = f(x=yx, mu=u[j,[0,2]], nu=u[j,[1,3]])
    data = SepiaData(x_sim=u, y_sim=np.log(y), y_ind_sim=yx)
    data.transform_xt()
    data.standardize_y()
    data.create_K_basis(n_pc=0.99)
    print("Created basis")
    model = SepiaModel(data)
    model.verbose = False
    model.tune_step_sizes(50, 20, verbose=False)
    print("\n")
    model.do_mcmc(n_samp)
    n_pred=Utrue.shape[0]
    pred_samples = model.get_samples(nburn=int(.1*n_samp))
    pred = SepiaEmulatorPrediction(x_pred=Utrueunit, samples=pred_samples, model=model, storeMuSigma=True)
    predy = np.exp(pred.get_y())
    meany = np.mean(predy,0)[0]
    s2s = np.var(predy,0)[0]
    sumsq = 0
    for j in range(meany.shape[0]):
        sumsq += pow(meany[j] - YY[j], 2)
        metrics[i,0] = math.pow(sumsq/meany.shape[0], 0.5)
    scores = np.zeros(yx.shape[0])
    for j in range(meany.shape[0]):
        sigma = pow(s2s[j], 0.5)
        if (sigma == 0):
            scores[j] = np.nan
            continue
        z = (YY[j] - meany[j])/sigma
        scores[j] = sigma*(-1/(math.pow(math.pi, 0.5)) + 2*norm.pdf(z) + z*(2*norm.cdf(z)-1))
    metrics[i,1] = np.nanmean(scores)
    with open('sepia_2dlogit_metrics.csv', 'w', newline='') as file:
        # Create a CSV writer object
        writer = csv.writer(file)  
        # Write each row of data to the CSV file
        writer.writerows(metrics)
    print("\n")
    print("Finished iteration " + str(i) + "\n")
