# python3 -m venv ~/py_envs
# source ~/py_envs/bin/activate

import csv
import math
import matplotlib.pyplot as plt
import numpy as np
import random
import time
from sepia.SepiaModel import SepiaModel
from sepia.SepiaData import SepiaData
from sepia.SepiaPredict import SepiaEmulatorPrediction

x = np.genfromtxt('calib_params_ibex.csv', delimiter=',', skip_header=1) ## pmfp, ratio
y = np.transpose(np.genfromtxt('ibex_responses.csv', delimiter=',', skip_header=1)) ## 66 vectors of 16200
yx = np.genfromtxt('xyz_grid.csv', delimiter=',', skip_header=1) ## xyz grid
n_samp = 10000

# Testing the effect of adding more simulation runs (with fullMatrices=True when creating the basis)
times = np.zeros(len(x)-1)
for i in range(1, len(x)):
    inds = range(i+1)
    tic = time.time()
    data = SepiaData(x_sim = x[inds,], y_sim = y[inds,], y_ind_sim = yx)
    data.transform_xt()
    data.standardize_y()
    data.create_K_basis(n_pc=0.99)
    model = SepiaModel(data)
    toc = time.time()
    times[i-1] = toc - tic
    print("Finished iteration " + str(i) + "\n")

plt.plot(range(1, len(x)), times, marker='o', linestyle='-', color='b')
plt.xlabel('Number of simulation runs included')
plt.ylabel('Execution time to create basis')
plt.grid(True)
plt.show()

# Testing the effect of making simulation output higher and higher dimension (with fullMatrices=True when creating the basis)
ns = range(100, 16300, 100)
times = np.zeros(len(ns))
for i in range(len(ns)):
    inds = random.sample(range(16200), ns[i])
    tic = time.time()
    data = SepiaData(x_sim = x, y_sim = y[:,inds], y_ind_sim = yx[inds,])
    data.transform_xt()
    data.standardize_y()
    data.create_K_basis(n_pc=0.99)
    model = SepiaModel(data)
    toc = time.time()
    times[i] = toc - tic
    print("Finished iteration " + str(i) + "\n")

plt.plot(ns, times, marker='o', linestyle='-', color='b')
plt.xlabel('Dimension of simulator output')
plt.ylabel('Execution time to create basis')
plt.grid(True)
plt.show()

# Testing the effect of increasing one dimension of a matrix 
ns = range(100000, 1000000000, 100000)
create_times = np.zeros(len(ns))
sum_times = np.zeros(len(ns))
for i in range(len(ns)):
    tic = time.time()
    iter_mat = np.random.random((66, ns[i]))
    toc = time.time()
    create_times[i] = toc - tic
    tic = time.time()
    iter_sum = np.sum(iter_mat)
    toc = time.time()
    sum_times[i] = toc - tic
    print(ns[i])

plt.plot(ns[1:50], create_times[1:50], marker='o', linestyle='-', color='b',
 label="matrix creation")
plt.plot(ns[1:50], sum_times[1:50], marker='o', linestyle='-', color='r',
 label="matrix sum")
plt.xlabel('Length of simulator output')
plt.ylabel('Execution time')
plt.grid(True)
plt.legend(loc="upper left")
plt.show()

# Testing the effect of increasing both dimensions of a matrix 
ns = range(100, 1000000000, 100)
create_times = np.zeros(len(ns))
sum_times = np.zeros(len(ns))
for i in range(len(ns)):
    tic = time.time()
    iter_mat = np.random.random((66, ns[i]))
    toc = time.time()
    create_times[i] = toc - tic
    tic = time.time()
    for j in range(ns[i]):
        iter_sum = np.sum(iter_mat)
    toc = time.time()
    sum_times[i] = toc - tic
    print(ns[i])

plt.plot(ns[1:99], create_times[1:99], marker='o', linestyle='-', color='b',
 label="matrix creation")
plt.plot(ns[1:99], sum_times[1:99], marker='o', linestyle='-', color='r',
 label="matrix sum")
plt.xlabel('Dimension of square matrix')
plt.ylabel('Execution time')
plt.grid(True)
plt.legend(loc="upper left")
plt.show()

# Write files of different sizes to test loading them in
ns = range(100, 20100, 100)
for i in range(len(ns)):
    iter_mat = np.random.random((66, ns[i]))
    fn = 'python_time_test_n' + str(ns[i]) + ".csv"
    with open(fn, 'w', newline='') as file:
        # Create a CSV writer object
        writer = csv.writer(file)  
        # Write each row of data to the CSV file
        writer.writerows(iter_mat)
    print(ns[i])

readin_times = np.zeros(len(ns))
ns = range(100, 20100, 100)
for i in range(len(ns)):
    fn = 'python_time_test_n' + str(ns[i]) + ".csv"
    tic = time.time()
    dat = np.genfromtxt(fn, delimiter=',') ## 66 vectors of length n
    toc = time.time()
    readin_times[i] = toc - tic
    print(ns[i])

plt.plot(ns, readin_times, marker='o', linestyle='-', color='b',
 label="file input time")
plt.xlabel('Length of simulator output')
plt.ylabel('Execution time')
plt.grid(True)
plt.legend(loc="upper left")
plt.show()