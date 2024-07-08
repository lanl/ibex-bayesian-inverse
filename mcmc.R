library(GPvecchia)
library(GpGp)
library(laGP)
library(tidyverse)

## Loads in scaled Vecchia approximation code
source('../vecchia_scaled.R')

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
  step=NA, thrds=2, vb=FALSE, debug=FALSE, true_u=NA, true_logscl=NA) {

  id <- sample(100000:999999, 1)
  ## create objects to hold posterior samples and other metrics
  u <- uprops <- matrix(data=NA, nrow=nmcmcs, ncol=ncol(Um))
  logscl <- logsclprops <- matrix(data=NA, nrow=nmcmcs, ncol=1)

  if (debug) {
    urates <- logsclrates <- rep(NA, nmcmcs-1)
  }
  covars <- list()
  colnames(u) <- colnames(Um)
  lls <- rep(NA, nmcmcs)
  XUm <- cbind(Xm, Um)
  if (gpmeth=="svecchia") {
    fit <- fit_scaled(y=Zm, inputs=as.matrix(cbind(Xm, Um)), nug=1e-4, ms=25)
  } else {
    fit <- NA
  }
  ## initialize chains
  if (any(is.na(true_u))) {
    u[1,] <- uprops[1,] <- apply(Um, 2, mean)
  } else {
    u <- uprops <- matrix(rep(true_u, nrow(u)), ncol=length(true_u),
     byrow=TRUE)
  }
  if (is.na(true_logscl)) {
    logscl[1,] <- logsclprops[1,] <- 0
  } else {
    logscl <- rep(true_logscl, nrow(logscl))
  }

  XX <- cbind(Xf, u[1,,drop=FALSE])
  colnames(XX) <- c("x", "y", "z", "pmfp", "ratio")
  ## TODO: value of m should be user specified
  lhat_curr <- predictions_scaled(fit, as.matrix(XX), m=25, joint=FALSE,
    predvar=FALSE)
  lls[1] <- sum(Zf*log(lhat_curr) - lhat_curr)
  ## establish covariance for proposals
  pvar <- ifelse(is.na(step), 0.1, step)
  pcovar <- covars[[1]] <- matrix(c(pvar, 0, 0, pvar), nrow=2, byrow=TRUE)
  pmin <- apply(Um, 2, min)
  pmax <- apply(Um, 2, max)
  tic <- proc.time()[3]

  for (t in 2:nmcmcs) {
    if (vb) print(paste("Started iteration:", t))

    ###########################################################################
    ## SAMPLE CALIBRATION PARAMETERS U
    ### Propose u_prime and calculate proposal ratio
    if (any(is.na(true_u))) {
      up <- propose_u(curr=u[t-1,], method="tmvnorm", pmin=pmin, pmax=pmax,
        pcovar=pcovar)
      uprops[t,] <- up$prop
      if (vb) {
        print(paste("Iteration proposal (calib params):", up$prop[1],
         up$prop[2]))
      }
      ### Predict simulator output at u_prime using fitted surrogate
      if (gpmeth=="svecchia") {
        ## Use scaled Vecchia GP
        XX <- cbind(Xf, up$prop)
        colnames(XX) <- c("x", "y", "z", "pmfp", "ratio")
        ## TODO: value of m should be user specified
        lhatp <- predictions_scaled(fit, as.matrix(XX), m=25, joint=FALSE,
          predvar=FALSE)
      } else {
        ## Use laGP
        d <- darg(NULL, cbind(Xm, Um))
        if (vb) print(paste0("Lengthscale prior information: ", "MLE: ", d$mle,
          "; start: ", d$start, "; max: ", d$max, "; min: ", d$min, "; a: ",
          d$ab[1], "; b: ", d$ab[2]))
        sfit <- aGPsep(X=XUm, Z=Zm, XX=cbind(Xf, up$prop), method=gpmeth, end=end,
          omp.threads=thrds, d=d, verb=0)
        lhatp <- sfit$mean
      }
      ### Calculate proposed likelihood
      llp <- sum(Zf*log(Of$time*(lhatp*exp(logscl[t-1]) + Of$bg)) -
       Of$time*(lhatp*exp(logscl[t-1]) + Of$bg))
      ### Calculate prior on u (calibration parameters)
      lpp <- dbeta(up$prop[1], shape1=1.1, shape2=1.1, log=TRUE) +
       dbeta(up$prop[2], shape1=1.1, shape2=1.1, log=TRUE)
      lp_curr <- dbeta(u[t-1,1], shape1=1.1, shape2=1.1, log=TRUE) +
       dbeta(u[t-1,2], shape1=1.1, shape2=1.1, log=TRUE)
      ### Calculate Metropolis-Hastings ratio
      ### { L(xp|Y)*p(xp)*g(xt|xp) } / { L(xt|Y)*p(xt)*g(xp|xt) }
      lmh <- llp - lls[t-1] + lpp - lp_curr + up$pr
      ## accept or reject
      if (lmh > log(runif(n=1))) {
        u[t,] <- up$prop
        lls[t] <- llp
        lhat_curr <- lhatp
      } else {
        u[t,] <- u[t-1,]
        lls[t] <- lls[t-1]
      }
    } else {
      lls[t] <- sum(Zf*log(Of$time*(lhat_curr*exp(logscl[t-1]) + Of$bg)) -
       Of$time*(lhat_curr*exp(logscl[t-1]) + Of$bg))
    }

    ###########################################################################

    ## update u proposal covariance
    if (any(is.na(true_u)) && t %% 100 == 0) {
      accepted_vals <- u[(t-winsize+1):t,] %>%
        as.data.frame() %>% tidyr::drop_na() %>% dplyr::distinct()
      pcovar <- covars[[floor(t / 100)+1]] <- as.matrix(cov(accepted_vals))
    }

    ###########################################################################
    ## SAMPLE SCALE PARAMETER
    if (is.na(true_logscl)) {
      ## TODO: adaptive estimation for sd
      logsclp <- propose_logscl(curr=logscl[t-1], sd=0.1)
      logsclprops[t] <- logsclp$prop
      if (vb) print(paste("Iteration proposal (scale):", exp(logsclp$prop)))
      ### Calculate proposed likelihood
      llp <- sum(Zf*log(Of$time*(lhat_curr*exp(logsclp$prop) + Of$bg)) -
       Of$time*(lhat_curr*exp(logsclp$prop) + Of$bg))
      ### Calculate prior on log scale parameter
      lpp <- dnorm(x=logsclp$prop, log=TRUE)
      lp_curr <- dnorm(x=logscl[t-1], log=TRUE)
      ### Calculate Metropolis-Hastings ratio
      ### { L(xp|Y)*p(xp)*g(xt|xp) } / { L(xt|Y)*p(xt)*g(xp|xt) }
      lmh <- llp - lls[t-1] + lpp - lp_curr + logsclp$pr
      ## accept or reject
      if (lmh > log(runif(n=1))) {
        logscl[t,] <- logsclp$prop
        lls[t] <- llp
      } else {
        logscl[t,] <- logscl[t-1]
        lls[t] <- lls[t-1]
      }
    } else {
      lls[t] <- sum(Zf*log(Of$time*(lhat_curr*exp(true_logscl) + Of$bg)) -
       Of$time*(lhat_curr*exp(true_logscl) + Of$bg))
    }
    ###########################################################################

    toc <- proc.time()[3]
    ## calculate acceptance rates
    if (debug) {
      winsize <- min(t, 100)
      accepted_u <- u[(t-winsize+1):t,] %>%
        as.data.frame() %>% tidyr::drop_na() %>% dplyr::distinct()
      urates[t-10] <- nrow(accepted_u)/winsize
      accepted_logscl <- unique(logscl[(t-winsize+1):t])
      logsclrates[t-10] <- length(accepted_logscl)/winsize
    }

    if (vb && t %% 1 == 0) {
      print(paste("Finished iteration", t))
      print(paste("Time elapsed:", toc-tic))
      print(paste("Iteration sample (calib params):", drop(u[t,1]),
        drop(u[t,2])))
      print(paste("Iteration sample (scl):", exp(logscl[t])))
    }
    if (debug && t %% 100 == 0) {
      temp_res <- list(u=u, logscl=logscl, lls=lls,
        props=list(u=uprops, logscl=logsclprops), covars=covars)
      if (debug) temp_res$rates <- list(u=urates, logscl=logsclrates)
      saveRDS(temp_res, file=paste0("../results/temp_mcmc_", id, ".rds"))
    }
  }
  ret <- list(u=u, logscl=logscl, lls=lls, props=list(u=uprops, logscl=logsclprops),
    covars=covars, time=proc.time()[3]-tic)
  if (debug) ret$rates <- list(u=urates, logscl=logsclrates)
  return(ret)
}
