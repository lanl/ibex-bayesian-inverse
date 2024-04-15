library(GPvecchia)
library(GpGp)
library(laGP)
library(tidyverse)

## Loads in scaled Vecchia approximation code
Sys.setenv('https_proxy'='http://proxyout.lanl.gov:8080')
source('https://raw.githubusercontent.com/katzfuss-group/scaledVecchia/master/vecchia_scaled.R')

###############################################################################
# Runs Markov chain Monte Carlo (McMC) for computer model calibration. Uses
# scaled Vecchia approximation or locally approximated Gaussian processes to
# fit surrogate model from the computer model output.
#
# @param Xm matrix of field inputs for the computer model
# @param Um matrix of calibration parameter combinations for the computer model
# @param Zm vector of outputs from each computer model run
# @param Xf matrix of inputs for runs of the field experiment
# @param Zf vector of outputs from each field experiment
# @param end neighborhood size for each GP if using laGP
# @param gpmeth method for GP fit (e.g. svecchia or if using laGP, nn)
# @param nmcmcs number of McMC iterations
# @param step step size for random walk proposal distribution
# @param thrds if using laGP, number of threads to pass for OpenMP to use
# @param vb flag indicating verbose output
# @param debug flag indicating if more output should be saved for debugging
#
# @return list with posterior samples of calibration parameters, along with
# likelihoods, proposals, acceptance rates, and covariances.
###############################################################################
mcmc <- function(Xm, Um, Zm, Xf, Zf, Of, end=NA, gpmeth="nn", nmcmcs=10000,
                 step=NA, thrds=2, vb=FALSE, debug=FALSE) {
  ## create objects to hold posterior samples and other metrics
  u <- matrix(data=NA, nrow=nmcmcs, ncol=ncol(Um))
  props <- matrix(data=NA, nrow=nmcmcs, ncol=ncol(Um))
  rates <- rep(NA, nmcmcs-min(nmcmcs, 10))
  covars <- list()
  colnames(u) <- colnames(Um)
  lls <- rep(NA, nmcmcs)
  XUm <- cbind(Xm, Um)
  if (gpmeth=="svecchia") {
    fit <- fit_scaled(y=Zm, inputs=as.matrix(cbind(Xm, Um)), nug=1e-4, ms=25)
  }
  ## initialize chains
  u[1,] <- props[1,] <- apply(Um, 2, mean)
  lhatc <- runif(length(Zf))
  lls[1] <- sum(Zf*log(lhatc) - lhatc)
  ## establish covariance for proposals
  pvar <- ifelse(is.na(step), 0.1, step)
  pcovar <- covars[[1]] <- matrix(c(pvar, 0, 0, pvar), nrow=2, byrow=TRUE)
  pmin <- apply(Um, 2, min)
  pmax <- apply(Um, 2, max)
  tic <- proc.time()[3]
  for (t in 2:nmcmcs) {
    if (vb) print(paste("Started iteration:", t))
    ## Propose u_prime and calculate proposal ratio
    up <- proposal(curr=u[t-1,], method="tmvnorm", pmin=pmin, pmax=pmax,
      pcovar=pcovar)
    props[t,] <- up$prop
    if (vb) print(paste("Iteration proposal:", up$prop[1], up$prop[2]))
    if (gpmeth=="svecchia") {
      ## surrogate model predictions using scaled vecchia fit
      XX <- cbind(Xf, up$prop)
      colnames(XX) <- c("x", "y", "z", "pmfp", "ratio")
      ## TODO: value of m should be user specified
      lhatp <- predictions_scaled(fit, as.matrix(XX), m=25, joint=FALSE,
        predvar=FALSE)
    } else {
      ## surrogate model predictions using a locally approximated GP
      d <- darg(NULL, cbind(Xm, Um))
      if (vb) print(paste0("Lengthscale prior information: ", "MLE: ", d$mle,
        "; start: ", d$start, "; max: ", d$max, "; min: ", d$min, "; a: ",
        d$ab[1], "; b: ", d$ab[2]))
      sfit <- aGPsep(X=XUm, Z=Zm, XX=cbind(Xf, up$prop), method=gpmeth, end=end,
        omp.threads=thrds, d=d, verb=0)
      lhatp <- sfit$mean
    }
    ## calculate MH ratio under truncated multivariate normal proposal
    ## TODO: is my proposal ratio correct?
    llp <- sum(Zf*log(Of$time*(lhatp + Of$bg)) - Of$time*(lhatp + Of$bg))
    lmh <- llp - lls[t-1] + up$pr
    ## accept or reject
    if (lmh > log(runif(n=1))) {
      u[t,] <- up$prop
      lls[t] <- llp
    } else {
      u[t,] <- u[t-1,]
      lls[t] <- lls[t-1]
    }
    ## update proposal covariance
    if (debug && t > 10) {
      winsize <- ifelse(t > 100, 100, t)
      accepted_vals <- u[(t-winsize+1):t,] %>%
        as.data.frame() %>% tidyr::drop_na() %>% dplyr::distinct()
      pcovar <- covars[[t-10]] <- as.matrix(cov(accepted_vals))
      rates[t-10] <- nrow(accepted_vals)/winsize
    }
    toc <- proc.time()[3]
    if (vb && t %% 1 == 0) {
      print(paste("Finished iteration", t))
      print(paste("Time elapsed:", toc-tic))
      print(paste("Iteration sample:", drop(u[t,1]), drop(u[t,2])))
    }
    if (t %% 100 == 0) {
      temp_res <- list(u=u, lls=lls, props=props, rates=rates, covars=covars)
      saveRDS(temp_res, file="../results/temp_mcmc.rds")
    }
  }
  return(list(u=u, lls=lls, props=props, rates=rates, covars=covars,
    time=proc.time()[3]-tic))
}
