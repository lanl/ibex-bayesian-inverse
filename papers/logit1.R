###############################################################################
###############################################################################
## Figures for toy examples showing current set up of Bayesian inverse
## problems and our proposed framework for Poisson Bayesian inverse problems
###############################################################################
###############################################################################

###############################################################################
## FIGURE 2: Toy 1D example for standard Bayesian inverse problem. A logistic
## curve with two input parameters to vary the shape
## DATA NEEDED:
###############################################################################

library(laGP)
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
calib_params <- randomLHS(n=30, k=2)
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

Y <- as.vector(ym)
X <- cbind(matrix(rep(xm, ncol(ym)), ncol=1),
  matrix(rep(calib_params_unit[,1], each=length(xm)), ncol=1),
  matrix(rep(calib_params_unit[,2], each=length(xm)), ncol=1))
g <- garg(list(mle=TRUE, max=1), Y)
d <- darg(list(mle=TRUE, max=0.25), X)
gpi <- newGPsep(X, Y, d=d$start, g=g$start, dK=TRUE)
mle <- jmleGPsep(gpi, c(d$min, d$max), c(g$min, g$max), d$ab, g$ab)

## initialize chains
u[1,] <- uprops[1,] <- c(0.5, 0.5)

XX <- matrix(xf, ncol=1)
DX <- distance(XX)
Sigma <- exp(-DX) + diag(1, nrow(DX))
SigmaInv <- solve(Sigma)
SigmaDet <- determinant(Sigma)$modulus
attributes(SigmaDet) <- NULL
tau2hat <- drop(t(yf) %*% SigmaInv %*% yf / length(yf))

lhat_curr <- predGPsep(gpi, XX=cbind(XX, matrix(rep(uprops[1,], nrow(XX)),
  ncol=2, byrow=TRUE)), lite=TRUE, nonug=TRUE)$mean
lls[1] <- -0.5*nrow(XX)*tau2hat - 0.5*determinant(Sigma)$modulus -
  drop(0.5*t(yf - lhat_curr) %*% SigmaInv %*% (yf-lhat_curr))/tau2hat

accept <- 1
lmhs <- rep(NA, nmcmcs)
mod_lhatps <- mod_lhats_accept <- matrix(NA, nrow=length(xm), ncol=nmcmcs)
mod_lhatps[,1] <- mod_lhats_accept[,1] <- predGPsep(gpi,
  XX=cbind(matrix(xm, ncol=1), matrix(rep(uprops[1,], length(xm)),
  ncol=2, byrow=TRUE)), lite=TRUE, nonug=TRUE)$mean
for (t in 2:nmcmcs) {

  ###########################################################################
  ## SAMPLE CALIBRATION PARAMETERS U
  ### Propose u_prime and calculate proposal ratio
  up <- propose_u(curr=u[t-1,], method="tmvnorm", pmin=rep(0, 2), pmax=rep(1, 2),
    pcovar=matrix(c(0.05, 0, 0, 0.05), byrow=TRUE, ncol=2))
  uprops[t,] <- up$prop
  ## Evaluate surrogate at u_prime
  lhatp <- predGPsep(gpi, XX=cbind(XX, matrix(rep(uprops[t,], nrow(XX)),
    ncol=2, byrow=TRUE)), lite=TRUE, nonug=TRUE)$mean
  mod_lhatps[,t] <- predGPsep(gpi,
    XX=cbind(matrix(xm, ncol=1), matrix(rep(uprops[t,], length(xm)),
    ncol=2, byrow=TRUE)), lite=TRUE, nonug=TRUE)$mean

  ### Calculate proposed likelihood
  llp <- -0.5*nrow(XX)*tau2hat - 0.5*SigmaDet -
    drop(0.5*t(yf - lhatp) %*% SigmaInv %*% (yf-lhatp))/tau2hat

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
    mod_lhats_accept[,t] <- mod_lhatps[,t]
    accept <- accept + 1
  } else {
    u[t,] <- u[t-1,]
    lls[t] <- lls[t-1]
    mod_lhats_accept[,t] <- mod_lhats_accept[,t-1]
  }

  if (t %% 100 == 0) {
    print(paste("Finished iteration", t))
  }
}
deleteGPsep(gpi)

## Figure 2 (left panel)
ylims <- range(c(mod_lhats_accept[,seq(15001, 20000, by=10)], yf, lam_m))
ylims[1] <- ylims[1]-0.5
ylims[2] <- ylims[2]+0.5
pdf("logit1_obs.pdf", width=5, height=5)
matplot(x=xm, y=ym, type="l", col="lightgrey", lty=1,
  lwd=1.5, xlab="X", ylab="Y", ylim=ylims, mgp=c(2,0.75,0))
points(x=as.vector(t(xf)), y=yf, col=2, pch=8)
lines(x=xm, y=lam_m, lwd=2, lty=2)
legend("topleft", c("observations", "model runs", "truth"),
  col=c(2, "lightgrey", 1), pch=c(8, NA, NA), lty=c(NA, 1, 2), lwd=c(1, 1.5, 2),
  bg="white", cex=1.05)
dev.off()

## Figure 2 (middle panel)
## Visualize model evaluations
pdf("logit1_est.pdf", width=5, height=5)
matplot(x=xm, y=mod_lhats_accept[,seq(15001, 20000, by=10)], type="l", lty=1,
  col=adjustcolor("lightgrey", alpha.f=0.3), xlab="X", yaxt="n", ylim=ylims,
  mgp=c(2,0.75,0))
points(x=xf, y=yf, col=2, pch=8)
u_postmean <- apply(u[seq(15001, 20000, by=10),], 2, mean)
lines(x=xm, y=f(xm, u_postmean[1]*uranges[1]+umins[1],
  u_postmean[2]*uranges[2]+umins[2]), col=4, lwd=2, lty=4)
lines(x=xm, y=lam_m, lty=2, lwd=2)
legend("topleft", c(expression("model runs at u"^(t)),
  expression("model at " * bar(u)["post"])),
  col=c("lightgrey", 4), lty=c(1, 4), lwd=2, pch=c(rep(NA, 2)), bg="white",
  y.intersp=1.3, cex=1.05)
dev.off()

## Figure 2 (right panel)
## Visualize posterior draws of u
pdf("logit1_post_draws.pdf", width=5, height=5)
plot(x=u[seq(10001, 20000, by=10),1]*uranges[1]+umins[1],
 y=u[seq(10001, 20000, by=10),2]*uranges[2]+umins[2], xlab=expression(u[1]),
 ylab=expression(u[2]), col="lightgrey", mgp=c(2,0.75,0))
points(x=true_mu, y=true_nu, col=3, pch=8, lwd=2, cex=1.5)
points(x=u_postmean[1]*uranges[1]+umins[1],
  y=u_postmean[2]*uranges[2]+umins[2], col=4, pch=9, lwd=2, cex=1.5)
legend("topleft", c("posterior draws of u", "posterior mean", "truth"),
  col=c("lightgrey", 4, 3), lty=NA, lwd=2, pch=c(1, 9, 8), cex=1.05, bg="white")
dev.off()

###############################################################################
## FIGURE 3: Toy 1D example for standard Bayesian inverse problem. A logistic
## curve with two input parameters to vary the shape
## DATA NEEDED:
###############################################################################

library(laGP)
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
set.seed(51997)
repsf <- sample(3:7, 8, replace=TRUE)
nm <- 20

## Set up field data
xf <- rep(seq(0, 1, length=nf), repsf)
xf <- c()
for (i in 1:nf) {
  xf <- c(xf, rep(seq(0, 1, length=nf)[i], repsf[i]))
}
lam <- f(x=xf, mu=true_mu, nu=true_nu)
yf <- rpois(lam, lam)

## Set up computer model data
xm <- seq(0, 1, length=nm)
lam_m <- f(x=xm, mu=true_mu, nu=true_nu)
calib_params <- randomLHS(n=30, k=2)
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

Y <- as.vector(ym)
X <- cbind(matrix(rep(xm, ncol(ym)), ncol=1),
  matrix(rep(calib_params_unit[,1], each=length(xm)), ncol=1),
  matrix(rep(calib_params_unit[,2], each=length(xm)), ncol=1))
g <- garg(list(mle=TRUE, max=1), Y)
d <- darg(list(mle=TRUE, max=0.25), X)
gpi <- newGPsep(X, log(Y), d=d$start, g=g$start, dK=TRUE)
mle <- jmleGPsep(gpi, c(d$min, d$max), c(g$min, g$max), d$ab, g$ab)

## initialize chains
u[1,] <- uprops[1,] <- c(0.5, 0.5)

XX <- matrix(xf, ncol=1)

lhat_curr <- exp(predGPsep(gpi, XX=cbind(XX, matrix(rep(uprops[1,], nrow(XX)),
  ncol=2, byrow=TRUE)), lite=TRUE, nonug=TRUE)$mean)
lls[1] <- sum(yf*log(lhat_curr) - lhat_curr)

accept <- 1
lmhs <- rep(NA, nmcmcs)
mod_lhatps <- mod_lhats_accept <- matrix(NA, nrow=length(xm), ncol=nmcmcs)
mod_lhatps[,1] <- mod_lhats_accept[,1] <- predGPsep(gpi,
  XX=cbind(matrix(xm, ncol=1), matrix(rep(uprops[1,], length(xm)),
  ncol=2, byrow=TRUE)), lite=TRUE, nonug=TRUE)$mean
for (t in 2:nmcmcs) {

  ###########################################################################
  ## SAMPLE CALIBRATION PARAMETERS U
  ### Propose u_prime and calculate proposal ratio
  up <- propose_u(curr=u[t-1,], method="tmvnorm", pmin=rep(0, 2), pmax=rep(1, 2),
    pcovar=matrix(c(0.15, 0, 0, 0.15), byrow=TRUE, ncol=2))
  uprops[t,] <- up$prop
  ## Evaluate simulator at u_prime
  lhatp <- exp(predGPsep(gpi, XX=cbind(XX, matrix(rep(uprops[t,], nrow(XX)),
    ncol=2, byrow=TRUE)), lite=TRUE, nonug=TRUE)$mean)
  mod_lhatps[,t] <- exp(predGPsep(gpi,
    XX=cbind(matrix(xm, ncol=1), matrix(rep(uprops[t,], length(xm)),
    ncol=2, byrow=TRUE)), lite=TRUE, nonug=TRUE)$mean)

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
    mod_lhats_accept[,t] <- mod_lhatps[,t]
    accept <- accept + 1
  } else {
    u[t,] <- u[t-1,]
    lls[t] <- lls[t-1]
    mod_lhats_accept[,t] <- mod_lhats_accept[,t-1]
  }
  if (t %% 100 == 0) {
    print(paste("Finished iteration", t))
  }
}
deleteGPsep(gpi)

## Figure 3 (left panel)
ylims <- range(c(mod_lhats_accept[,seq(15001, 20000, by=10)], yf, lam_m))
ylims[1] <- ylims[1]-0.5
ylims[2] <- ylims[2]+0.5
par(mfrow=c(1,1), mar=c(5.1, 4.1, 0.05, 0.05))
pdf("logit1_pois_obs.pdf", width=5, height=5)
matplot(x=xm, y=ym, type="l", col="lightgrey", lty=1, lwd=1.5, xlab="X",
  ylab=expression(lambda), ylim=ylims, mgp=c(2,0.75,0))
points(x=as.vector(t(xf)), y=yf, col=2, pch=8)
est_means <- rep(NA, nf)
y_ind <- 1
for (i in 1:nf) {
  est_means[i] <- mean(yf[y_ind:(y_ind+repsf[i]-1)])
  y_ind <- y_ind + repsf[i]
}
points(x=seq(0, 1, length=nf), y=est_means, col=1, bg=2, pch=21)
lines(x=xm, y=lam_m, lwd=2, lty=2)
legend("topleft", c("observed counts", "counts/exposure", "model runs",
  "truth"), col=c(2, 1, "lightgrey", 1), pch=c(8, 21, NA, NA),
  lty=c(NA, NA, 1:2), lwd=c(1, 1, 1.5, 2), pt.bg=c(NA, 2, NA, NA),
  bg="white", cex=1.05)
dev.off()

## Figure 3 (middle panel)
## Visualize model evaluations
par(mfrow=c(1,1), mar=c(5.1, 4.1, 0.05, 0.05))
pdf("logit1_pois_est.pdf", width=5, height=5)
matplot(x=xm, y=mod_lhats_accept[,seq(15001, 20000, by=10)], type="l", lty=1,
  col="lightgrey", xlab="X", yaxt="n", ylim=ylims, mgp=c(2,0.75,0))
points(x=xf, y=yf, col=2, pch=8)
points(x=seq(0, 1, length=nf), y=est_means, col=1, bg=2, pch=21)
u_postmean <- apply(u[seq(15001, 20000, by=10),], 2, mean)
lines(x=xm, y=f(xm, u_postmean[1]*uranges[1]+umins[1],
  u_postmean[2]*uranges[2]+umins[2]), col=4, lwd=2, lty=4)
lines(x=xm, y=lam_m, lty=2, lwd=2)
legend("topleft", c(expression("model runs at u"^(t)),
  expression("model at " * bar(u)["post"])),
  col=c("lightgrey", 4), lty=c(1, 4), lwd=2, pch=c(rep(NA, 2)), bg="white",
  y.intersp=1.3, cex=1.05)
dev.off()

## Figure 3 (right panel)
## Visualize posterior draws of u
par(mfrow=c(1,1), mar=c(5.1, 4.1, 0.05, 0.05))
pdf("logit1_pois_post_draws.pdf", width=5, height=5)
plot(x=u[seq(10001, 20000, by=10),1]*uranges[1]+umins[1],
 y=u[seq(10001, 20000, by=10),2]*uranges[2]+umins[2], xlab=expression(u[1]),
 ylab=expression(u[2]), col="lightgrey", mgp=c(2,0.75,0))
points(x=true_mu, y=true_nu, col=3, pch=8, lwd=2, cex=1.5)
points(x=u_postmean[1]*uranges[1]+umins[1],
  y=u_postmean[2]*uranges[2]+umins[2], col=4, pch=9, lwd=2, cex=1.5)
legend("topleft", c("posterior draws of u", "posterior mean", "truth"),
  col=c("lightgrey", 4, 3), lty=NA, lwd=2, pch=c(1, 9, 8), cex=1.05, bg="white")
dev.off()
