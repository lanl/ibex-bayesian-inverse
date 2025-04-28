####################################################################
## FIGURES ?/?: Toy 2D example for Poisson calibration
####################################################################

library(lhs)

f <- function(x, mu, nu) {
  if (is.null(nrow(x)) || ncol(x) != 2) {
    stop("function requires an nx2 matrix of parameters")
  }
  if (length(mu) != 2 || length(nu) != 2) {
    stop("function requires two vectors of length two")
  }
  t1 <- mu[1]*exp(mu[1]*x[,1]-7)/(nu[1]+exp(mu[1]*x[,1]-7))
  t2 <- mu[2]*exp(mu[2]*x[,2]-3)/(nu[2]+exp(mu[2]*x[,2]-3))
  return(t1+t2)
}

true_mu <- c(11, 8)
true_nu <- c(0.5, 2.3)

set.seed(233807)
xf <- randomLHS(n=40, k=2)
xf <- rbind(xf, xf, xf)
lam <- as.vector(t(f(x=xf, mu=true_mu, nu=true_nu)))
yf <- rpois(lam, lam)

xm <- seq(0, 1, length=20)
Xm <- as.matrix(expand.grid(xm, xm))
lam_m <- f(x=Xm, mu=true_mu, nu=true_nu)
calib_params <- randomLHS(n=100, k=4)
mu_range <- 10
mu_min <- 5
nu_range <- 3
nu_min <- 0
calib_params[,1] <- calib_params[,1]*mu_range + mu_min
calib_params[,2] <- calib_params[,2]*nu_range + nu_min
calib_params[,3] <- calib_params[,3]*nu_range + nu_min
calib_params[,4] <- calib_params[,4]*nu_range + nu_min
colnames(calib_params) <- c("mu1", "nu1", "mu2", "nu2")
ym <- matrix(NA, ncol=nrow(calib_params), nrow=nrow(Xm))
for (i in 1:nrow(calib_params)) {
  mu <- calib_params[i,c(1,3)]
  nu <- calib_params[i,c(2,4)]
  ym[,i] <- f(x=Xm, mu=mu, nu=nu)
}

library(classInt)
ylims <- range(yf, ym)
cols <- heat.colors(128)
intervals <- classIntervals(c(yf, ym), 128, style="fixed",
 fixedBreaks=seq(ylims[1], ylims[2], length=129))
point_cols <- cols[findInterval(yf, intervals$brks, all.inside=T)]
leg_img <- as.raster(matrix(rev(cols),ncol=1))
pdf("logit2_obs.pdf", width=6, height=5)
par(mar=c(4.75,4.25,0.4,1), oma=c(2,1.75,1,2));
layout(matrix(c(1, 2), nrow=1, ncol=2, byrow=TRUE), width=c(1, 0.2))
plot(x=jitter(xf[,1],amount=0.02), y=jitter(xf[,2],amount=0.02),
  cex=1.5, pch=19, col=point_cols, xlab="x1", ylab="x2")
par(mar=c(4.75,0.4,1,1))
plot(c(0,2),c(0,1),type='n', axes=FALSE, xlab="", ylab="", main="Counts")
text(x=1.5, y=seq(0,1,length=4),
 labels=format(round(intervals$brks[seq(1, length(intervals$brks), length=4)], 0)),
 cex=1.25)
rasterImage(leg_img,0,0,1,1)
dev.off()

# mu1=11.2076745, nu1=2.8304068, mu2=0.4027237, nu2=1.8098630
rand_model_rns <- sample(1:ncol(ym), 2)
bks <- seq(ylims[1], ylims[2], length=129)
pdf("logit2_model1.pdf", width=6, height=5)
par(mar=c(4.75,4.25,0.4,1), oma=c(2,1.75,1,2));
layout(matrix(c(1, 2), nrow=1, ncol=2, byrow=TRUE), width=c(1, 0.2))
image(x=xm, y=xm, z=matrix(ym[,rand_model_rns[1]], ncol=length(xm)),
 col=cols, breaks=bks, xlab="x1", ylab="x2")
par(mar=c(4.75,0.4,1,1))
plot(c(0,2),c(0,1),type='n', axes=FALSE, xlab="", ylab="", main=expression(paste(lambda)))
text(x=1.5, y=seq(0,1,length=4),
 labels=format(round(intervals$brks[seq(1, length(intervals$brks), length=4)], 0)),
 cex=1.25)
rasterImage(leg_img,0,0,1,1)
dev.off()

# mu1=14.6179921, nu1=0.2186919, mu2=2.7904711, nu2=0.5410672
pdf("logit2_model2.pdf", width=6, height=5)
par(mar=c(4.75,4.25,0.4,1), oma=c(2,1.75,1,2));
layout(matrix(c(1, 2), nrow=1, ncol=2, byrow=TRUE), width=c(1, 0.2))
image(x=xm, y=xm, z=matrix(ym[,rand_model_rns[2]], ncol=length(xm)),
 col=cols, breaks=bks, xlab="x1", ylab="x2")
par(mar=c(4.75,0.4,1,1))
plot(c(0,2),c(0,1),type='n', axes=FALSE, xlab="", ylab="", main=expression(paste(lambda)))
text(x=1.5, y=seq(0,1,length=4),
 labels=format(round(intervals$brks[seq(1, length(intervals$brks), length=4)], 0)),
 cex=1.25)
rasterImage(leg_img,0,0,1,1)
dev.off()

library(GPvecchia)
library(GpGp)
library(laGP)
library(tidyverse)

## Loads in scaled Vecchia approximation code
source('../helper.R')
source('../vecchia_scaled.R')

nmcmcs <- 10000
u <- uprops <- matrix(data=NA, nrow=nmcmcs, ncol=ncol(calib_params))

colnames(u) <- colnames(uprops) <- colnames(calib_params)
lls <- rep(NA, nmcmcs)
umins <- umaxs <- uranges <- rep(NA, 2)

for (i in 1:ncol(calib_params)) {
  umins[i] <- min(calib_params[,i])
  umaxs[i] <- max(calib_params[,i])
  uranges[i] <- diff(range(calib_params[,i]))
  calib_params[,i] <- (calib_params[,i] - umins[i])/uranges[i]
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

fit <- fit_scaled(y=as.vector(ym), inputs=XUm, nug=1e-4, ms=25)

## initialize chains
u[1,] <- uprops[1,] <- c(0.5, 0.5)

XX <- matrix(as.vector(t(xf)), ncol=1)
XX <- cbind(XX, rep(NA, nrow(XX)), rep(NA, nrow(XX)))
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
    pcovar=matrix(c(0.15, 0, 0, 0.15), byrow=TRUE, ncol=2))
  uprops[t,] <- up$prop
  # print(paste("Iteration proposal (calib params):", up$prop[1], up$prop[2]))
  ### Predict simulator output at u_prime using fitted surrogate
  ## Use scaled Vecchia GP
  XX <- matrix(as.vector(t(xf)), ncol=1)
  XX <- cbind(XX, rep(NA, nrow(XX)), rep(NA, nrow(XX)))
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

  # print(paste("Current likelihood:", lls[t-1])
  # print(paste("Proposed likelihood:", llp)

  ## accept or reject
  if (lmh > log(runif(n=1))) {
    u[t,] <- up$prop
    lls[t] <- llp
    lhat_curr <- lhatp
    accept <- accept + 1
    # print("ACCEPTED!!")
  } else {
    u[t,] <- u[t-1,]
    lls[t] <- lls[t-1]
    # print("REJECTED!!")
  }
  # t <- t + 1
  print(paste("Finished iteration", t))
}

## Visualize this:
plot(x=xm, y=f(xm, u[t-1,1]*uranges[1]+umins[1], u[t-1,2]*uranges[2]+umins[2]),
 type="l", lwd=2, xlab="X", ylab="lambda", ylim=ylims)
lines(x=XX[1:8,1], y=lhat_curr[1:8], col=2, lwd=2, lty=2)
points(x=as.vector(t(xf)), y=yf)
lines(x=XX[1:8,1], y=lhatp[1:8], col=3, lwd=2, lty=3)
# lines(x=xm, y=f(xm, uprops[t,1]*uranges[1]+umins[1],
 # uprops[t,2]*uranges[2]+umins[2]), col=3, lty=3, lwd=2)
lines(x=xm, y=lam_m, col=4, lty=4, lwd=2)
legend("topleft", c("current u at model X", "current u at field X",
  "proposed u at field X", "true lambda", "field obs"),
  col=c(1:4, 1), lty=c(1:4, NA), lwd=2, pch=c(rep(NA, 4), 1))


par(mfrow=c(1,3))
plot(x=xm, y=f(xm, u[t-1,1]*uranges[1]+umins[1], u[t-1,2]*uranges[2]+umins[2]),
   type="n", lwd=2, xlab="X", ylab="lambda", ylim=ylims)
matplot(x=xm, y=draws, type="l", col="lightgrey", add=TRUE, lty=1)
matplot(x=XX[1:8,1], y=fdraws, col="lightpink", add=TRUE, pch=1)
points(x=as.vector(t(xf)), y=yf)
lines(x=xm, y=lam_m, col=4, lty=4, lwd=2)
lines(x=xm, y=apply(draws, 1, mean), col=1, lty=1, lwd=3)
legend("topleft", c("truth", "posterior draws at model X", "posterior draws at field X", "observed counts"), col=c(4, "lightgrey", "lightpink", 1), lty=c(4, 1, NA, NA), pch=c(NA, NA, 1, 1), cex=1.5, bg="white")
plot(x=1:900, y=u[seq(1001, 10000, by=10),1]*uranges[1]+umins[1], type="l", xlab="iteration", ylab="mu")
abline(h=true_mu, col=4, lty=4, lwd=2)
abline(h=mean(u[seq(1001, 10000, by=10),1]*uranges[1]+umins[1]), col=3, lwd=2)
abline(h=quantile(u[seq(1001, 10000, by=10),1]*uranges[1]+umins[1], prob=0.025), col=3, lty=3, lwd=2)
abline(h=quantile(u[seq(1001, 10000, by=10),1]*uranges[1]+umins[1], prob=0.975), col=3, lty=3, lwd=2)
plot(x=1:900, y=u[seq(1001, 10000, by=10),2]*uranges[2]+umins[2], type="l", xlab="iteration", ylab="nu")
abline(h=true_nu, col=4, lty=4, lwd=2)
abline(h=mean(u[seq(1001, 10000, by=10),2]*uranges[2]+umins[2]), col=3, lwd=2)
abline(h=quantile(u[seq(1001, 10000, by=10),2]*uranges[2]+umins[2], prob=0.025), col=3, lty=3, lwd=2)
abline(h=quantile(u[seq(1001, 10000, by=10),2]*uranges[2]+umins[2], prob=0.975), col=3, lty=3, lwd=2)