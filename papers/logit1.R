####################################################################
## FIGURES 2/3: Toy 1D example for standard calibration
####################################################################

library(lhs)
library(plgp)

source("../helper.R")

f <- function(x, mu, nu) {
  if (is.null(nrow(x))) {
    x <- matrix(x, ncol=1)
  }
  mu*exp(mu*x-5)/(nu+exp(mu*x-5))
}

true_mu <- 10
true_nu <- 1
nf <- 8
repsf <- 4
nm <- 20

## Set up field data
xf <- rep(seq(0, 1, length=nf), repsf)
lam <- f(x=xf, mu=true_mu, nu=true_nu)
set.seed(9812987)
yf <- rnorm(n=length(lam), mean=lam)

## Set up computer model data
xm <- seq(0, 1, length=nm)
lam_m <- f(x=xm, mu=true_mu, nu=true_nu)
calib_params <- randomLHS(n=60, k=2)
mu_range <- 10
mu_min <- 5
mu_max <- 15
nu_range <- 2
nu_min <- 0
nu_max <- 2
calib_params[,1] <- calib_params[,1]*mu_range + mu_min
calib_params[,2] <- calib_params[,2]*nu_range + nu_min
colnames(calib_params) <- c("mu", "nu")
ym <- matrix(NA, ncol=nrow(calib_params), nrow=length(xm))
for (i in 1:nrow(calib_params)) {
  mu <- calib_params[i,1]
  nu <- calib_params[i,2]
  ym[,i] <- f(x=xm, mu=mu, nu=nu)
}

ylims <- range(yf, ym)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 0.2, 0.2))
pdf("logit1_obs.pdf", width=5, height=5)
matplot(x=xm, y=ym[,sample(1:ncol(ym), 10)], type="l", col="lightgrey", lty=1,
  lwd=1.5, xlab="X", ylab="Y", ylim=ylims)
points(x=as.vector(t(xf)), y=yf, col=2, pch=8)
lines(x=xm, y=lam_m, lwd=2, lty=2)
legend("topleft", c("observations", "sample computer model runs", "truth"),
  col=c(2, "lightgrey", 1), pch=c(8, NA, NA), lty=c(NA, 1, 2), lwd=c(1, 1.5, 2),
  bg="white", cex=0.75)
dev.off()

nmcmcs <- 20000
u <- uprops <- matrix(data=NA, nrow=nmcmcs, ncol=ncol(calib_params))

colnames(u) <- colnames(uprops) <- colnames(calib_params)
lls <- rep(NA, nmcmcs)
umins <- umaxs <- uranges <- rep(NA, ncol(calib_params))

calib_params_unit <- calib_params
for (i in 1:ncol(calib_params)) {
  umins[i] <- min(calib_params[,i])
  umaxs[i] <- max(calib_params[,i])
  uranges[i] <- diff(range(calib_params[,i]))
  calib_params_unit[,i] <- (calib_params[,i] - umins[i])/uranges[i]
}

## initialize chains
u[1,] <- uprops[1,] <- c(0.5, 0.5)

XX <- matrix(xf, ncol=1)
DX <- distance(XX)
Sigma <- exp(-DX) + diag(1, nrow(DX))
SigmaInv <- solve(Sigma)
SigmaDet <- determinant(Sigma)$modulus
attributes(SigmaDet) <- NULL
tau2hat <- drop(t(yf) %*% SigmaInv %*% yf / length(yf))

lhat_curr <- f(x=XX, mu=u[1,1]*uranges[1]+umins[1], nu=u[1,2]*uranges[2]+umins[2])
lls[1] <- -0.5*nrow(XX)*tau2hat - 0.5*determinant(Sigma)$modulus - drop(0.5*t(yf - lhat_curr) %*% SigmaInv %*% (yf-lhat_curr))/tau2hat

accept <- 1
lmhs <- rep(NA, nmcmcs)
mod_lhatps <- matrix(NA, nrow=length(xm), ncol=nmcmcs)
mod_lhatps[,1] <- f(x=xm, mu=uprops[1,1]*uranges[1]+umins[1],
 nu=uprops[1,2]*uranges[2]+umins[2])
for (t in 2:nmcmcs) {

  ###########################################################################
  ## SAMPLE CALIBRATION PARAMETERS U
  ### Propose u_prime and calculate proposal ratio
  up <- propose_u(curr=u[t-1,], method="tmvnorm", pmin=rep(0, 2), pmax=rep(1, 2),
    pcovar=matrix(c(0.15, 0, 0, 0.15), byrow=TRUE, ncol=2))
  uprops[t,] <- up$prop
  ## Evaluate simulator at u_prime
  lhatp <- f(x=XX, mu=uprops[t,1]*uranges[1]+umins[1],
    nu=uprops[t,2]*uranges[2]+umins[2])
  mod_lhatps[,t] <- f(x=xm, mu=uprops[t,1]*uranges[1]+umins[1],
    nu=uprops[t,2]*uranges[2]+umins[2])

  ### Calculate proposed likelihood
  llp <- -0.5*nrow(XX)*tau2hat - 0.5*SigmaDet - drop(0.5*t(yf - lhatp) %*% SigmaInv %*% (yf-lhatp))/tau2hat

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
  if (t %% 100 == 0) {
    print(paste("Finished iteration", t))
  }
}

## Visualize model evaluations:
par(mfrow=c(1,1), mar=c(5.1, 4.1, 0.2, 0.2))
pdf("logit1_est.pdf", width=5, height=5)
ylims <- range(c(mod_lhatps, yf, lam_m))
matplot(x=xm, y=mod_lhatps[,seq(15001, 20000, by=10)], type="l", lty=1,
  col="lightgrey", xlab="X", ylab="Y", ylim=ylims)
points(x=xf, y=yf, col=2, pch=8)
u_postmean <- apply(u[seq(15001, 20000, by=10),], 2, mean)
lines(x=xm, y=f(xm, u_postmean[1]*uranges[1]+umins[1],
  u_postmean[2]*uranges[2]+umins[2]), col=4, lwd=2, lty=4)
lines(x=xm, y=lam_m, lty=2, lwd=2)
legend("topleft", c("observations", expression("model runs at u"^(t)),
  expression("model at " * bar(u)["post"]), "truth"),
  col=c(2, "lightgrey", 4, 1), lty=c(NA, 1, 4, 2), lwd=2, pch=c(8, rep(NA, 3)), bg="white",
  y.intersp=1.3, cex=0.75)
dev.off()

## Visualize posterior draws of u
par(mfrow=c(1,1), mar=c(5.1, 4.1, 0.2, 0.2))
pdf("logit1_post_draws.pdf", width=5, height=5)
plot(x=u[seq(10001, 20000, by=10),1]*uranges[1]+umins[1],
 y=u[seq(10001, 20000, by=10),2]*uranges[2]+umins[2], xlab=expression(mu),
 ylab=expression(nu), col="lightgrey")
points(x=true_mu, y=true_nu, col=3, pch=8, lwd=2, cex=1.5)
points(x=u_postmean[1]*uranges[1]+umins[1],
  y=u_postmean[2]*uranges[2]+umins[2], col=4, pch=9, lwd=2, cex=1.5)
legend("topleft", c("posterior draws of u", "posterior mean", "truth"),
  col=c("lightgrey", 4, 3), lty=NA, lwd=2, pch=c(1, 9, 8), cex=0.75, bg="white")
dev.off()


####################################################################
## FIGURES 2/3: Toy 1D example for Poisson calibration
####################################################################

library(lhs)

source("../helper.R")

f <- function(x, mu, nu) {
  if (is.null(nrow(x))) {
    x <- matrix(x, ncol=1)
  }
  mu*exp(mu*x-5)/(nu+exp(mu*x-5))
}

true_mu <- 10
true_nu <- 1
nf <- 8
repsf <- 4
nm <- 20

## Set up field data
xf <- rep(seq(0, 1, length=nf), repsf)
lam <- f(x=xf, mu=true_mu, nu=true_nu)
yf <- rpois(lam, lam)

## Set up computer model data
xm <- seq(0, 1, length=nm)
lam_m <- f(x=xm, mu=true_mu, nu=true_nu)
calib_params <- randomLHS(n=60, k=2)
mu_range <- 10
mu_min <- 5
mu_max <- 15
nu_range <- 2
nu_min <- 0
nu_max <- 2
calib_params[,1] <- calib_params[,1]*mu_range + mu_min
calib_params[,2] <- calib_params[,2]*nu_range + nu_min
colnames(calib_params) <- c("mu", "nu")
ym <- matrix(NA, ncol=nrow(calib_params), nrow=length(xm))
for (i in 1:nrow(calib_params)) {
  mu <- calib_params[i,1]
  nu <- calib_params[i,2]
  ym[,i] <- f(x=xm, mu=mu, nu=nu)
}

ylims <- range(yf, ym)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 0.2, 0.2))
pdf("logit1_obs.pdf", width=5, height=5)
matplot(x=xm, y=ym, type="l", col="lightgrey", lty=1, lwd=1.5, xlab="X",
  ylab=expression(lambda), ylim=ylims)
points(x=as.vector(t(xf)), y=yf)
legend("topleft", c("observed counts", "computer model output"),
  col=c(1, "lightgrey"), pch=c(1, NA), lty=c(NA, 1), lwd=c(1, 1.5), cex=1.25)
dev.off()

nmcmcs <- 20000
u <- uprops <- matrix(data=NA, nrow=nmcmcs, ncol=ncol(calib_params))

colnames(u) <- colnames(uprops) <- colnames(calib_params)
lls <- rep(NA, nmcmcs)
umins <- umaxs <- uranges <- rep(NA, ncol(calib_params))

calib_params_unit <- calib_params
for (i in 1:ncol(calib_params)) {
  umins[i] <- min(calib_params[,i])
  umaxs[i] <- max(calib_params[,i])
  uranges[i] <- diff(range(calib_params[,i]))
  calib_params_unit[,i] <- (calib_params[,i] - umins[i])/uranges[i]
}

## initialize chains
u[1,] <- uprops[1,] <- c(0.5, 0.5)

XX <- matrix(xf, ncol=1)

lhat_curr <- f(x=XX, mu=u[1,1]*uranges[1]+umins[1], nu=u[1,2]*uranges[2]+umins[2])
lls[1] <- sum(yf*log(lhat_curr) - lhat_curr)

accept <- 1
lmhs <- rep(NA, nmcmcs)
mod_lhatps <- matrix(NA, nrow=length(xm), ncol=nmcmcs)
mod_lhatps[,1] <- f(x=xm, mu=uprops[1,1]*uranges[1]+umins[1],
 nu=uprops[1,2]*uranges[2]+umins[2])
for (t in 2:nmcmcs) {

  ###########################################################################
  ## SAMPLE CALIBRATION PARAMETERS U
  ### Propose u_prime and calculate proposal ratio
  up <- propose_u(curr=u[t-1,], method="tmvnorm", pmin=rep(0, 2), pmax=rep(1, 2),
    pcovar=matrix(c(0.15, 0, 0, 0.15), byrow=TRUE, ncol=2))
  uprops[t,] <- up$prop
  ## Evaluate simulator at u_prime
  lhatp <- f(x=XX, mu=uprops[t,1]*uranges[1]+umins[1],
    nu=uprops[t,2]*uranges[2]+umins[2])
  mod_lhatps[,t] <- f(x=xm, mu=uprops[t,1]*uranges[1]+umins[1],
    nu=uprops[t,2]*uranges[2]+umins[2])

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
  if (t %% 100 == 0) {
    print(paste("Finished iteration", t))
  }
}

## Visualize model evaluations:
par(mfrow=c(1,1), mar=c(5.1, 4.1, 0.2, 0.2))
pdf("logit1_est.pdf", width=5, height=5)
ylims <- range(c(mod_lhatps, yf, lam_m))
matplot(x=xm, y=mod_lhatps[,seq(15001, 20000, by=10)], type="l", lty=1,
  col="lightgrey", xlab="X", ylab="lambda", ylim=ylims)
points(x=xf, y=yf, lwd=2)
u_postmean <- apply(u[seq(15001, 20000, by=10),], 2, mean)
lines(x=xm, y=f(xm, u_postmean[1]*uranges[1]+umins[1],
  u_postmean[2]*uranges[2]+umins[2]), col=2, lwd=2, lty=2)
lines(x=xm, y=lam_m, col=3, lty=3, lwd=2)
legend("topleft", c("field counts", "posterior draws of f(x)",
  "posterior mean estimate", "true lambda"),
  col=c(1, "lightgrey", 2, 3), lty=c(NA, 1:3), lwd=2, pch=c(1, rep(NA, 3)))
dev.off()

## Visualize posterior draws of u
par(mfrow=c(1,1), mar=c(5.1, 4.1, 0.2, 0.2))
pdf("logit1_post_draws.pdf", width=5, height=5)
plot(x=u[seq(10001, 20000, by=10),1]*uranges[1]+umins[1],
 y=u[seq(10001, 20000, by=10),2]*uranges[2]+umins[2], xlab=expression(mu),
 ylab=expression(nu))
points(x=true_mu, y=true_nu, col=2, pch=8, lwd=2, cex=1.5)
points(x=u_postmean[1]*uranges[1]+umins[1],
  y=u_postmean[2]*uranges[2]+umins[2], col=3, pch=8, lwd=2, cex=1.5)
legend("topleft", c("posterior draws of u", "posterior mean", "truth"),
  col=c(1:3), lty=NA, lwd=2, pch=c(1, 8, 8))
dev.off()

####################################################################
## FIGURE 5: 1D example for Poisson calibration w/ limited data
####################################################################

library(laGP)

## Uses same set up as above

nmcmcs <- 20000
u <- uprops <- matrix(data=NA, nrow=nmcmcs, ncol=ncol(calib_params))

colnames(u) <- colnames(uprops) <- colnames(calib_params)
lls <- rep(NA, nmcmcs)
umins <- umaxs <- uranges <- rep(NA, ncol(calib_params))

calib_params_unit <- calib_params
for (i in 1:ncol(calib_params)) {
  umins[i] <- min(calib_params[,i])
  umaxs[i] <- max(calib_params[,i])
  uranges[i] <- diff(range(calib_params[,i]))
  calib_params_unit[,i] <- (calib_params[,i] - umins[i])/uranges[i]
}

## initialize chains
u[1,] <- uprops[1,] <- c(0.5, 0.5)

Xm <- matrix(rep(xm, ncol(ym)), ncol=1)
Um <- matrix(NA, nrow=length(xm)*ncol(ym), ncol=ncol(calib_params_unit))
for (i in 1:nrow(calib_params_unit)) {
  for (j in 1:length(xm)) {
    Um[((i-1)*length(xm)+j),] <- calib_params_unit[i,]
  }
}

Xsurr <- cbind(Xm, Um)
tic <- proc.time()[3]
gpfit <- newGP(X=Xsurr, Z=c(ym), d=0.1, g=1e-6, dK=TRUE)
mle <- mleGP(gpfit, param="d")
toc <- proc.time()[3]
toc -tic

XXf <- matrix(NA, nrow=length(xf), ncol=3)
XXf[,1] <- xf
XXf[,2] <- (10 - umins[1])/uranges[1]
XXf[,3] <- (1.1 - umins[2])/uranges[2]

lhat_curr <- predGP(gpfit, XX=XXf, lite=TRUE)$mean
lls[1] <- sum(yf*log(lhat_curr) - lhat_curr)

accept <- 1
lhatps <- lhat_truth <- matrix(NA, nrow=8, ncol=nmcmcs)
for (t in 2:nmcmcs) {

  ###########################################################################
  ## SAMPLE CALIBRATION PARAMETERS U
  ### Propose u_prime and calculate proposal ratio
  up <- propose_u(curr=u[t-1,], method="tmvnorm", pmin=rep(0, 2), pmax=rep(1, 2),
    pcovar=matrix(c(0.15, 0, 0, 0.15), byrow=TRUE, ncol=2))
  uprops[t,] <- up$prop
  ## Evaluate simulator at u_prime
  XXiter <- XXf
  XXiter[,2] <- uprops[t,1]
  XXiter[,3] <- uprops[t,2]
  lhatp <- predGP(gpfit, XX=XXiter[1:8,], lite=TRUE)$mean
  ## Set negative values to 0.0001
  lhatp[which(lhatp<=0)] <- 0.0001
  lhatps[,t] <- lhatp
  lhatp <- rep(lhatp, 4)
  lhat_truth[,t] <- f(x=XXiter[1:8,1], mu=XXiter[1,2]*uranges[1]+umins[1],
    nu=XXiter[1,3]*uranges[2]+umins[2])

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
  if (t %% 100 == 0) {
    print(paste("Finished iteration", t))
  }
}

## Visualize model evaluations:
par(mfrow=c(1,1), mar=c(5.1, 4.1, 0.2, 0.2))
pdf("logit1_est_wsurr.pdf", width=5, height=5)
ylims <- range(c(mod_lhatps, yf, lam_m))
matplot(x=xm, y=mod_lhatps[,seq(15001, 20000, by=10)], type="l", lty=1,
  col="lightgrey", xlab="X", ylab="lambda", ylim=ylims)
points(x=xf, y=yf, lwd=2)
u_postmean <- apply(u[seq(15001, 20000, by=10),], 2, mean)
lines(x=xm, y=f(xm, u_postmean[1]*uranges[1]+umins[1],
  u_postmean[2]*uranges[2]+umins[2]), col=2, lwd=3, lty=2)
lines(x=xm, y=lam_m, col=3, lty=3, lwd=4)
legend("topleft", c("field counts", "posterior draws of f(x)",
  "posterior mean estimate", "true lambda"),
  col=c(1, "lightgrey", 2, 3), lty=c(NA, 1:3), lwd=3, pch=c(1, rep(NA, 3)))
dev.off()

## Visualize posterior draws of u
par(mfrow=c(1,1), mar=c(5.1, 4.1, 0.2, 0.2))
pdf("logit1_post_draws_wsurr.pdf", width=5, height=5)
plot(x=u[seq(10001, 20000, by=10),1]*uranges[1]+umins[1],
 y=u[seq(10001, 20000, by=10),2]*uranges[2]+umins[2], xlab=expression(mu),
 ylab=expression(nu))
points(x=true_mu, y=true_nu, col=2, pch=8, lwd=2, cex=1.5)
points(x=u_postmean[1]*uranges[1]+umins[1],
  y=u_postmean[2]*uranges[2]+umins[2], col=3, pch=8, lwd=2, cex=1.5)
legend("topleft", c("posterior draws of u", "posterior mean", "truth"),
  col=c(1:3), lty=NA, lwd=2, pch=c(1, 8, 8))
dev.off()
