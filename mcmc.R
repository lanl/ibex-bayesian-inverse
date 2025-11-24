###############################################################################
# Proposes value of model parameter(s) for next iteration of an MCMC. Proposals
# can come from a number of different distributions (depending on the value of
# method)
#
# @param ucurr the current sample for the model parameter in the McMC
# @param method method by which to propose a new value (e.g. unif, norm)
# @param umin minimum value the sampled parameter is allowed to be
# @param umax maximum value the sampled parameter is allowed to be
# @param ucov covariance matrix of parameters (for norm or truncated
# norm method)
#
# @return list with proposal and proposal ratio
###############################################################################
propose_u <- function(ucurr, method, umin=0, umax=1,
  ucov=diag(1, length(ucurr))) {

  } if (method=="unif") {
    ## Propose u_prime from a unit hypercube (centered around current value)
    mins <- apply(matrix(c(ucurr, 1-ucurr), ncol=2, byrow=FALSE), 1, min)
    maxs <- apply(matrix(c(ucurr, 1-ucurr), ncol=2, byrow=FALSE), 1, max)
    return(list(prop=t(mapply(runif, n=1, min=mins, max=maxs)), pr=0))
  } else if (method=="tmvnorm") {
    umin <- rep(umin, length(ucurr))
    umax <- rep(umax, length(ucurr))
    prop <- tmvtnorm::rtmvnorm(n=1, mean=ucurr, sigma=ucov, lower=umin,
      upper=umax, algorithm="rejection")
    pr <- tmvtnorm::dtmvnorm(ucurr, mean=drop(prop), sigma=ucov, lower=umin,
      upper=umax, log=TRUE) -
      tmvtnorm::dtmvnorm(drop(prop), mean=ucurr, sigma=ucov, lower=umin,
        upper=umax, log=TRUE)
    return(list(prop=prop, pr=pr))
  } else if (method=="norm") {
    prop <- rmvnorm(n=1, mean=ucurr, sigma=ucov)
    pr <- dmvnorm(ucurr, mean=drop(prop), sigma=ucov, log=TRUE) -
      dmvnorm(drop(prop), mean=ucurr, sigma=ucov, log=TRUE)
    return(list(prop=prop, pr=0))
  } else {
    stop("specified method not implemented")
  }
}

###############################################################################
# Proposes value of the log multiplicative scale discrepancy for the next
# iteration of an MCMC. Proposals come from a random walk
#
# @param lscl_curr the current sample in the MCMC
# @param sd step size for the random walk (i.e. standard deviation of proposal)
#
# @return list with proposal and proposal ratio
###############################################################################
propose_logscl <- function(lscl_curr, sd) {
  prop <- rnorm(n=1, mean=lscl_curr, sd=sd)
  pr <- dnorm(lscl_curr, mean=prop, sd=sd, log=TRUE) -
    dnorm(prop, mean=lscl_curr, sd=sd, log=TRUE)
  return(list(prop=prop, pr=0))
}

###############################################################################
# Calculates the log likelihood of observed Poisson counts under variable
# exposure, given an underlying mean function
#
# @param lambda given mean under which to evaluate the likelihood of counts
# @param obs observed counts
# @param reps exposure time under which the counts were observed
#
# @return the log likelihood of observed counts
###############################################################################
ll_pois <- function(lambda, obs, reps) {
  sum(obs*log(lambda) + obs*log(reps) - lambda*reps)
}
