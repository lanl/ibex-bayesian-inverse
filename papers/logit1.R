#################

#################

f <- function(x, mu, nu) {
  if (is.null(nrow(x))) {
    x <- matrix(x, ncol=1)
  }
  mu*exp(mu*x-5)/(nu+exp(mu*x-5))
}

true_mu <- 10
true_nu <- 1
nf <- 8
nm <- 50
repsf <- 4

xf <- rep(seq(0, 1, length=nf), repsf)
lam <- f(x=xf, mu=true_mu, nu=true_nu)
yf <- rpois(lam, lam)

xm <- seq(0, 1, length=nm)
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
matplot(x=xm, y=ym, type="l", col="lightgrey", lty=1, lwd=1.5, xlab="X",
  ylab=expression(lambda), ylim=ylims)
points(x=as.vector(t(xf)), y=yf)
legend("topleft", c("observed counts", "computer model output"),
  col=c(1, "lightgrey"), pch=c(1, NA), lty=c(NA, 1), lwd=c(1, 1.5), cex=1.25)
dev.off()

nmcmcs <- 10000
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

pmin <- apply(calib_params, 2, min)
pmax <- apply(calib_params, 2, max)
lhat_curr <- f(x=XX, mu=u[1,1]*uranges[1]+umins[1], nu=u[1,2]*uranges[2]+umins[2])
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
  ## Evaluate simulator at u_prime
  lhatp <- f(x=XX, mu=uprops[t,1]*uranges[1]+umins[1], nu=uprops[t,2]*uranges[2]+umins[2])

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

## Visualize this:
plot(x=xm, y=f(xm, u[t-1,1]*uranges[1]+umins[1], u[t-1,2]*uranges[2]+umins[2]),
 type="l", lwd=2, xlab="X", ylab="lambda", ylim=ylims)
lines(x=xf[1:8], y=lhat_curr[1:8], col=2, lwd=2, lty=2)
points(x=xf, y=yf)
lines(x=XX[1:8], y=lhatp[1:8], col=3, lwd=2, lty=3)
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
