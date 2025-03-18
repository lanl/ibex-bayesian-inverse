#################

#################

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
nm <- 50

## Set up field data
xf <- rep(seq(0, 1, length=nf), repsf)
lam <- f(x=xf, mu=true_mu, nu=true_nu)
yf <- rpois(lam, lam)

## Set up computer model data
xm <- seq(0, 1, length=nm)
lam_m <- f(x=xm, mu=true_mu, nu=true_nu)
mus <- seq(5, 15, length=4)
nus <- seq(0.25, 1.75, length=3)
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

for (i in 1:ncol(calib_params)) {
  umins[i] <- min(calib_params[,i])
  umaxs[i] <- max(calib_params[,i])
  uranges[i] <- diff(range(calib_params[,i]))
  calib_params[,i] <- (calib_params[,i] - umins[i])/uranges[i]
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
