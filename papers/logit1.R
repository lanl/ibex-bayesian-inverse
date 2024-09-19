
f <- function(x, mu, nu) {
  mu*exp(mu*x-5)/(nu+exp(mu*x-5))
}

true_mu <- 10
true_nu <- 1

xf <- seq(0, 1, length=8)
xf <- rbind(xf, xf, xf, xf)
lam <- f(x=xf, mu=true_mu, nu=true_nu)
yf <- rpois(lam, lam)

xm <- seq(0, 1, length=50)
lam_m <- f(x=xm, mu=true_mu, nu=true_nu)
mus <- seq(5, 15, length=5)
nus <- seq(0.25, 1.75, length=5)
calib_params <- as.matrix(expand.grid(mus, nus))
colnames(calib_params) <- c("mu", "nu")

ym <- matrix(NA, ncol=nrow(calib_params), nrow=length(xm))
for (i in 1:nrow(calib_params)) {
  mu <- calib_params[i,1]
  nu <- calib_params[i,2]
  ym[,i] <- f(x=xm, mu=mu, nu=nu)
}

ylims <- range(yf, ym)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 0.2, 0.2))
pdf("logit1_examp.pdf", width=5, height=5)
plot(x=xf, y=yf, xlab="X", ylab=expression(lambda), ylim=ylims)
matplot(x=xm, y=ym, type="l", col="lightgrey", lty=1, add=TRUE)
lines(x=xm, y=lam_m, type="l", col=2, lty=2, lwd=2)
legend("topleft", c("counts", "simulation", "true mean"),
  col=c(1, "lightgrey", 2), pch=c(1, NA, NA), lty=c(NA, 1, 2), lwd=c(1, 1, 2),
  cex=0.85)
dev.off()

library(GPvecchia)
library(GpGp)
library(laGP)
library(tidyverse)

## Loads in scaled Vecchia approximation code
source('../vecchia_scaled.R')

nmcmcs <- 10000
u <- uprops <- matrix(data=NA, nrow=nmcmcs, ncol=ncol(calib_params))

colnames(u) <- colnames(calib_params)
lls <- rep(NA, nmcmcs)

for (i in 1:ncol(calib_params)) {
  Umin <- min(calib_params[,i])
  Umax <- max(calib_params[,i])
  Urange <- diff(range(calib_params[,i]))
  calib_params[,i] <- (calib_params[,i] - Umin)/Urange
}

Xm <- matrix(rep(xm, nrow(calib_params)), ncol=1)
Xm <- cbind(Xm, rep(NA, 1250), rep(NA, 1250))
row <- 1
for (i in 1:nrow(calib_params)) {
  for (j in 1:length(xm)) {
    Xm[row,2:3] <- calib_params[i,] 
    row <- row + 1 
  }
}
XUm <- Xm

fit <- fit_scaled(y=ym, inputs=XUm, nug=1e-4, ms=25)

## initialize chains
u[1,] <- uprops[1,] <- c(0.5, 0.5)

XX <- cbind(xf, rep(NA, nrow(xf)), rep(NA, nrow(xf)))
for (i in 1:nrow(XX)) {
  XX[i,2:3] <- u[1,]
}
colnames(XX) <- c("x", "mu", "nu")

pmin <- apply(calib_params, 2, min)
pmax <- apply(calib_params, 2, max)
## TODO: value of m should be user specified
lhat_curr <- predictions_scaled(fit, as.matrix(XX), m=25, joint=FALSE,
  predvar=FALSE)
lhat_curr[lhat_curr <= 0] <- 0.1
lls[1] <- sum(yf*log(lhat_curr) - lhat_curr)

accept <- 1
lmhs <- rep(NA, nmcmcs)
for (t in 2:nmcmcs) {

  ###########################################################################
  ## SAMPLE CALIBRATION PARAMETERS U
  ### Propose u_prime and calculate proposal ratio
  up <- propose_u(curr=u[t-1,], method="tmvnorm", pmin=pmin, pmax=pmax,
    pcovar=matrix(c(1, 0, 0, 1), byrow=TRUE, ncol=2))
  uprops[t,] <- up$prop
  # print(paste("Iteration proposal (calib params):", up$prop[1], up$prop[2]))
  ### Predict simulator output at u_prime using fitted surrogate
  ## Use scaled Vecchia GP
  XX <- cbind(xf, rep(NA, nrow(xf)), rep(NA, nrow(xf)))
  for (i in 1:nrow(XX)) {
    XX[i,2:3] <- up$prop
  }
  colnames(XX) <- c("x", "mu", "nu")
  ## TODO: value of m should be user specified
  lhatp <- predictions_scaled(fit, as.matrix(XX), m=25, joint=FALSE,
    predvar=FALSE)
  lhatp[lhatp < 0] <- 0.1

  ### Calculate proposed likelihood
  llp <- sum(yf*log(lhatp) - lhatp)
  ### Calculate prior on u (calibration parameters)
  lpp <- dbeta(up$prop[1], shape1=1.1, shape2=1.1, log=TRUE) +
   dbeta(up$prop[2], shape1=1.1, shape2=1.1, log=TRUE)
  lp_curr <- dbeta(u[t-1,1], shape1=1.1, shape2=1.1, log=TRUE) +
   dbeta(u[t-1,2], shape1=1.1, shape2=1.1, log=TRUE)
  ### Calculate Metropolis-Hastings ratio
  ### { L(xp|Y)*p(xp)*g(xt|xp) } / { L(xt|Y)*p(xt)*g(xp|xt) }
  lmh <- llp - lls[t-1] + lpp - lp_curr + up$pr
  lmhs[t] <- lmh

  ## accept or reject
  if (lmh > log(runif(n=1))) {
    u[t,] <- up$prop
    lls[t] <- llp
    lhat_curr <- lhatp
    accept <- accept + 1
  } else {
    u[t,] <- u[t-1,]
    lls[t] <- lls[t-1]
  }

  print(paste("Finished iteration", t))
}