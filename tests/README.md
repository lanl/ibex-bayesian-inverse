# Test Scripts

This directory contains several scripts that test the performance and
functionality of the code within this respository. Below are descriptions,
required data files, and the corresponding figures for each set of scripts.
Each set can include bash scripts and/or R and Python code files.

### ibex_bayes*

- Description: Solves the inverse problem for ENA counts collected by the IBEX
  satellite in a Bayesian framework for Poisson observations. The file
  `ibex_bayes_test.R` contains code that calls the main function in this
  repository and returns posterior samples of computer model parameter. Two
  shell scripts, `ibex_bayes_real.sh` and `ibex_bayes_synth.sh`, execute the
  framework on real and synthetic satellite data, respectively.
- Required data: `sims.csv`, `ibex_real.csv`, `synth_sat_data.csv`
- Figures: 9, 10, 11, 13, 14, 15, 16

### real_data_cv*

- Description: Performs 10-fold cross-validation for IBEX satellite data
  collected between the years 2009-2011. Satellite data is split into 10
  separate folds. For each iteration, nine folds are used to estimate model
  parameters in our Poisson Bayesian inverse problem framework. Predictions of
  ENA rates are made at locations in the held out fold. CRPS is calculated to
  measure performance at different model parameter values.
- Required data: `sims.csv`, `ibex_real.csv`
- Figures: 12

### scale_disc_test*

- Description: Conducts a simulation to determine if our Poisson Bayesian
  inverse framework is able to recover a multiplicative scale discrepancy
  between the computer model and satellite observations. Ten different scale
  terms are considered. In each instance, synthetic satellite data (with a
  scale discrepancy) is generated within the script.
- Required data: `sims.csv`, `ibex_real.csv`
- Figures: 17

### surrogate_*

- Description: Code testing the ability of various surrogates to effectively
  model the IBEX simulator response. Surrogate models include R packages `laGP`
  and `deepgp` (with one layer), and our preferred method, Scaled Vecchia.
  Predictive performance is measured through RMSE and proper uncertainty
  quantification is assessed via CRPS. Computational thriftiness is also
  evaluated for each method. Varying dimensions of simulator response and
  different number of training runs available are both tested for their effect
  on execution time.
- Required data: `sims.csv`
- Figures: 6, 7, 8

### sepia*

- Description: Python code testing the ability of the SEPIA package to
  effectively model the IBEX simulator response. As above, predictive
  performance is measured through RMSE, proper uncertainty quantification is
  assessed via CRPS, and computational thriftiness is evaluated by varying
  simulator response dimension and training set size. Documentation on SEPIA's
  implementation and user guidance can be found here:
  https://sepia-lanl.readthedocs.io/en/latest/
- Required data: `sims.csv`
- Figures: 6, 7, 8
