library(deepgp)
library(laGP)
library(lhs)

source("../helper.R")
source('../vecchia_scaled.R')

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

start <- 1
method <- "all"
seed <- 6756781
mcs <- 30

## read in the command line arguments
## run with: R CMD BATCH '--args seed=1 method=all start=1' surrogate_test.R
args <- commandArgs(TRUE)
if (length(args) > 0) {
  for(i in 1:length(args)) {
    eval(parse(text=args[[i]]))
  }
}
print(args)
settings <- list(seed=seed, method=method, start=start)
print(settings)

# x <- seq(0, 1, length=5)
x <- seq(0, 1, length=90)
X <- as.matrix(expand.grid(x, x))
colnames(X) <- NULL

Utrue <- matrix(c(11, 0.5, 8, 2.3), nrow=1)
YY <- f(x=X, mu=Utrue[c(1,3)], nu=Utrue[c(2,4)])
saveRDS(Utrue, paste0("surrogate_test_2dlogit_settings",
 format(Sys.time(), "%Y%m%d"), ".rds"))

mu_min <- 5
mu_max <- 15
nu_min <- 0
nu_max <- 3

Utrueunit <- Utrue
Utrueunit[1] <- (Utrue[1] - mu_min) / (mu_max - mu_min)
Utrueunit[2] <- (Utrue[2] - nu_min) / (nu_max - nu_min)
Utrueunit[3] <- (Utrue[3] - mu_min) / (mu_max - mu_min)
Utrueunit[4] <- (Utrue[4] - nu_min) / (nu_max - nu_min)
XX <- matrix(NA, nrow=nrow(X), ncol=ncol(X)+ncol(Utrueunit))
XX[,1:ncol(X)] <- X
XX[,(ncol(X)+1):ncol(XX)] <- rep(Utrueunit, each=nrow(XX))

## Calculating metrics
rmses <- crps <- fit_times <- pred_times <- matrix(NA, ncol=3, nrow=mcs)
colnames(rmses) <- colnames(crps) <- colnames(fit_times) <-
  colnames(pred_times) <- c("svecchia", "lagp", "deepgp")

for (i in 1:mcs) {

  U <- randomLHS(n=125, k=4)
  # U <- randomLHS(n=10, k=4)
  Uunit <- U
  U[,1] <- U[,1] * (mu_max - mu_min) + mu_min
  U[,2] <- U[,2] * (nu_max - nu_min) + nu_min
  U[,3] <- U[,3] * (mu_max - mu_min) + mu_min
  U[,4] <- U[,4] * (nu_max - nu_min) + nu_min
  write.table(Uunit, file=paste0("surrogate_2dlogit_U", i, ".csv"), sep=",",
    row.names=FALSE, col.names=FALSE)

  Xtrain <- matrix(NA, nrow=nrow(X)*nrow(U), ncol=ncol(X)+ncol(U))
  for (j in 1:ncol(X)) {
    Xtrain[,j] <- X[,j]
  }
  for (j in (ncol(X)+1):(ncol(X)+ncol(U))) {
    Xtrain[,j] <- rep(Uunit[,j-ncol(X)], each=nrow(X))
  }
  Ytrain <- rep(NA, nrow(Xtrain))
  for (j in 1:nrow(U)) {
    beg <- (j-1)*nrow(X)+1
    end <- j*nrow(X)
    Ytrain[beg:end] <- log(f(x=X, mu=U[j,c(1,3)], nu=U[j,c(2,4)]))
  }

  tic <- proc.time()[3]
  svecfit <- fit_scaled(y=Ytrain, inputs=Xtrain, nug=1e-4, ms=75)
  toc <- proc.time()[3]
  fit_times[i,1] <- toc-tic
  print("Finished SVecchia fit")

  tic <- proc.time()[3]
  svecpreds <- predictions_scaled(svecfit, XX, m=75, joint=FALSE, predvar=TRUE)
  toc <- proc.time()[3]
  pred_times[i,1] <- toc-tic
  rmses[i,1] <- sqrt(mean((svecpreds$means - YY)^2))
  crps[i,1] <- crps(y=YY, mu=svecpreds$means, s2=svecpreds$vars)
  print("Finished SVecchia predictions")

  ## laGP
  tic <- proc.time()[3]
  d <- darg(NULL, Xtrain)
  lagppreds <- aGPsep(X=Xtrain, Z=Ytrain, XX=XX, omp.threads=8, verb=0,
    end=25, method="nn", d=d)
  toc <- proc.time()[3]
  pred_times[i,2] <- toc-tic
  rmses[i,2] <- sqrt(mean((exp(lagppreds$mean) - YY)^2))
  crps[i,2] <- crps(y=YY, mu=lagppreds$mean, s2=lagppreds$var)
  print("Finished laGP fit and predictions")

  # ## Fully Bayesian GP
  # tic <- proc.time()[3]
  # dgp1fit <- fit_one_layer(x=Xtrain, y=Ytrain, nmcmc=1000, vecchia=TRUE, m=10)
  # toc <- proc.time()[3]
  # fit_times[i,3] <- toc-tic
  # print("Finished deep gp fit")

  # tic <- proc.time()[3]
  # dgp1preds <- predict(trim(dgp1fit, burn=100, thin=5), x_new=X)
  # toc <- proc.time()[3]
  # pred_times[i,3] <- toc-tic
  # rmses[i,3] <- sqrt(mean((dgp1preds$mean - YY)^2))
  # crps[i,3] <- crps(y=YY, mu=dgp1preds$mean, s2=dgp1preds$s2)
  # print("Finished deep gp predictions")

  res <- list(fit_times=fit_times, pred_times=pred_times, rmse=rmses, crps=crps)
  saveRDS(res, paste0("surrogate_test_2dlogit", format(Sys.time(), "%Y%m%d"), ".rds"))
  print(paste0("Finished holdout iteration ", i))
}

res <- list(fit_times=fit_times, pred_times=pred_times, rmse=rmses, crps=crps)
saveRDS(res, paste0("surrogate_test_2dlogit", format(Sys.time(), "%Y%m%d"), ".rds"))
