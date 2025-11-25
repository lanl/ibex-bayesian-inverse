library(GPvecchia)
library(GpGp)
library(laGP)

source("check.R")
source("mcmc.R")
source("vecchia_scaled.R")

###############################################################################
# Runs Markov chain Monte Carlo (McMC) for a Poisson Bayesian inverse problem.
# This uses the scaled Vecchia approximation to fit a surrogate model over the
# output of a simulator
#
# @param xm matrix of field inputs for the simulator
# @param um matrix of model parameter combinations for the simulator
# @param ym vector of outputs from each simulator run
# @param xf matrix of inputs for runs of the field experiment
# @param yf vector of outputs from all field experiments
# @param e vector of exposure times from all field experiment
# @param lam0 vector of constants to add to lambda (for IBEX, background rates)
# @param T number of McMC iterations
# @param init list containing initial values (u)
# @param settings list containing settings (step size for random walk, beta
# parameter for prior on u, conditioning set size for Scaled Vecchia, and a
# flag indicating if a multiplicative scale discrepancy should be sampled)
# @param vb flag indicating verbose output
#
# @return list with posterior samples of model parameters and corresponding
# likelihoods
###############################################################################
pois_bayes_inv <- function(xm, um, ym, xf, yf, e, lam0, T,
  init=NULL, settings=NULL, vb=TRUE) {

  px <- ncol(xm)
  pu <- ncol(um)
  p <- px+pu
  nf <- nrow(xf)
  nm <- nrow(xm)

  ## create objects to hold posterior samples and other metrics
  us <- matrix(data=NA, nrow=T, ncol=pu)
  logscls <- rep(0, T)
  lls <- rep(NA, T)

  settings <- check_settings(settings)
  init <- check_init(init, um)
  ucov <- diag(settings$step, nrow=pu)

  xum <- as.matrix(cbind(xm, um))
  fit <- fit_scaled(y=ym, inputs=xum, nug=1e-4, ms=settings$m)
  us[1,] <- init$u
  xuf <- as.matrix(cbind(xf, matrix(us[1,], nrow=nf, ncol=pu, byrow=TRUE)))

  lam_curr <- predictions_scaled(fit, xuf, m=settings$m, joint=FALSE)
  lls[1] <- ll_pois(lambda=lam_curr+lam0, obs=yf, reps=e)

  for (t in 2:T) {
    # sample model parameter u
    ## propose u*
    ustar <- propose_u(ucurr=us[t-1,], method="tmvnorm", ucov=ucov)
    ## predict lambda*
    xuf[,(px+1):p] <- matrix(rep(ustar$prop, nf), nrow=nf, byrow=TRUE)
    lam_star <- predictions_scaled(fit, xuf, m=settings$m, joint=FALSE)

    ### calculate proposed likelihood
    ll_star <- ll_pois(lambda=lam_star*exp(logscls[t-1])+lam0, obs=yf, reps=e)

    ### calculate prior on u (model parameter)
    lp_star <- sum(dbeta(ustar$prop, shape1=settings$beta,
      shape2=settings$beta, log=TRUE))
    lp_curr <- sum(dbeta(us[t-1,], shape1=settings$beta,
      shape2=settings$beta, log=TRUE))

    ### calculate Metropolis-Hastings ratio
    lmh <- ll_star - lls[t-1] + lp_star - lp_curr + ustar$pr

    ## accept or reject
    if (lmh > log(runif(n=1))) {
      us[t,] <- ustar$prop
      lls[t] <- ll_star
      lam_curr <- lam_star
    } else {
      us[t,] <- us[t-1,]
      lls[t] <- lls[t-1]
    }

    # sample multiplicative discrepancy
    if (settings$sample_scl) {

      logscl_star <- propose_logscl(lscl_curr=logscls[t-1], sd=0.1)

      ### calculate proposed likelihood
      ll_star <- ll_pois(lambda=lam_star*exp(logscl_star$prop)+lam0,
        obs=yf, reps=e)

      ### calculate prior on the log scale
      lp_star <- dnorm(x=logscl_star$prop, log=TRUE)
      lp_curr <- dnorm(x=logscls[t-1], log=TRUE)

      ### calculate Metropolis-Hastings ratio
      lmh <- ll_star - lls[t-1] + lp_star - lp_curr + logscl_star$pr

      ## accept or reject
      if (lmh > log(runif(n=1))) {
        logscls[t] <- logscl_star$prop
        lls[t] <- ll_star
      } else {
        logscls[t] <- logscls[t-1]
        lls[t] <- lls[t-1]
      }
    }

    if (vb && t %% 10 == 0) {
      print(paste("Finished iteration", t))
    }
  }
  return(list(u=us, logscl=logscls, lls=lls))
}
