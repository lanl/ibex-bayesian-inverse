# The IBEX Bayesian Inverse Problem

Code solving the Bayesian inverse problem for satellite data collected by the
Interstellar Boundary Explorer (IBEX) and output from a computer model
attempting to represent the rate at which energetic neutral atoms (ENA) are
generated throughout the heliosphere. The goal is to estimate the posterior
distributions of the input parameter to the computer simulation, given observed
data from the IBEX satellite.

Descriptions of directories and files contained in the top-level directory are
below.

### Directories

- **data**: contains IBEX satellite data, computer model simulation output, and
  synthetic satellite data for testing purposes provided in this dataset
- **papers**: includes any content necessary to compile the manuscript along
  with `R` scripts for toy examples and figure generation
- **presentations**: code used to produce figures and documents for conference
  presentations
- **tests**: various bash, `Python`, and `R` scripts conducting tests to
  evaluate the performance and functionality of our proposed Poisson Bayesian
  inverse problem framework

## Files

- `pois_bayes_inv.R`: contains the function that implements the proposed
  Poisson Bayesian inverse problem framework. Requires computer model
  input/output and field observations
- `mcmc.R`: support functions for our Markov chain Monte Carlo (MCMC),
  including proposals and likelihood evaluations
- `helper.R`: contains functions that aid in the processing of IBEX satellite
  and simulator data, particulary in preparation for running our Bayesian
  inverse problem method
- `check.R`: input validation functions
- `vecchia_scaled.R`: original implementation of the Scaled Vecchia GP approximation for use as a surrogate model of the IBEX simulator. Original code found here: https://github.com/katzfuss-group/scaledVecchia

---

Copyright 2023 for **O4646**

This program is Open-Source under the BSD-3 License.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software
without specific prior written permission. THIS SOFTWARE IS PROVIDED BY THE
COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
