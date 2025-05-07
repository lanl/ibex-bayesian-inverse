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

seed <- 781691
inc_out <- TRUE

## read in the command line arguments
## run with: R CMD BATCH '--args seed=1 inc_out=1' surrogate_time_test_2dlogit.R
args <- commandArgs(TRUE)
if (length(args) > 0) {
  for(i in 1:length(args)) {
    eval(parse(text=args[[i]]))
  }
}

Utrue <- matrix(c(11, 0.5, 8, 2.3), nrow=1)

mu_min <- 5
mu_max <- 15
nu_min <- 0
nu_max <- 3

Utrueunit <- Utrue
Utrueunit[1] <- (Utrue[1] - mu_min) / (mu_max - mu_min)
Utrueunit[2] <- (Utrue[2] - nu_min) / (nu_max - nu_min)
Utrueunit[3] <- (Utrue[3] - mu_min) / (mu_max - mu_min)
Utrueunit[4] <- (Utrue[4] - nu_min) / (nu_max - nu_min)

if (inc_out) {
  U <- randomLHS(n=100, k=4)
  Uunit <- U
  U[,1] <- U[,1] * (mu_max - mu_min) + mu_min
  U[,2] <- U[,2] * (nu_max - nu_min) + nu_min
  U[,3] <- U[,3] * (mu_max - mu_min) + mu_min
  U[,4] <- U[,4] * (nu_max - nu_min) + nu_min

  ## Calculating metrics
  set.seed(seed)
  exp_pows <- 7:45

  ns <- c(round(10 + 1.25^exp_pows), seq(20000, 75000, by=5000))
  num_ns <- length(ns)
  outf <- "surrogate_time_test_2dlogit_"
  mcs <- 10

  ## Calculating metrics
  fit_times <- pred_times <- array(NA, dim=c(mcs, length(ns), 3))
  too_long <- c(FALSE, TRUE, TRUE)#rep(FALSE, 3)

  for (i in 1:num_ns) {
    ## select training data size
    n <- ns[i]
    if (i > length(exp_pows)) {
      cat("n = ", n)
    } else {
      cat("n = 1.25^", exp_pows[i], " = ", n, "\n", sep="")
    }

    N <- floor(sqrt(n))
    XX <- matrix(NA, nrow=N^2, ncol=6)
    xx <- seq(0, 1, length=N)
    XX[,1:2] <- as.matrix(expand.grid(xx, xx))
    XX[,3:6] <- rep(Utrueunit, each=nrow(XX))
    colnames(XX) <- NULL
    xtest_fn <- paste0("xtest_nout", n, ".csv")
    write.csv(XX, xtest_fn, row.names=FALSE)
    XX <- NULL

    for (m in 1:mcs) {
      xtrain_fn <- paste0("xtrain_nout", n, "_mc", m, ".csv")
      ytrain_fn <- paste0("ytrain_nout", n, "_mc", m, ".csv")

      X <- as.matrix(randomLHS(n=n, k=2))
      colnames(X) <- NULL

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
      write.csv(Xtrain, xtrain_fn, row.names=FALSE)
      write.csv(Ytrain, ytrain_fn, row.names=FALSE)
      Xtrain <- Ytrain <- NULL

      if (!too_long[1]) {
        tic <- proc.time()[3]
        Xtrain <- as.matrix(read.csv(xtrain_fn))
        colnames(Xtrain) <- NULL
        Ytrain <- drop(as.matrix(read.csv(ytrain_fn)))
        svecfit <- fit_scaled(y=Ytrain, inputs=Xtrain, nug=1e-4, ms=25)
        toc <- proc.time()[3]
        fit_times[m,i,1] <- toc-tic
        print("Finished SVecchia fit")

        tic <- proc.time()[3]
        XX <- as.matrix(read.csv(xtest_fn))
        colnames(XX) <- NULL
        svecpreds <- predictions_scaled(svecfit, XX, m=25, joint=FALSE, predvar=TRUE)
        toc <- proc.time()[3]
        pred_times[m,i,1] <- toc-tic
        print("Finished SVecchia predictions")
      }

      if (!too_long[2]) {
        tic <- proc.time()[3]
        Xtrain <- as.matrix(read.csv(xtrain_fn))
        colnames(Xtrain) <- NULL
        Ytrain <- drop(as.matrix(read.csv(ytrain_fn)))
        XX <- as.matrix(read.csv(xtest_fn))
        colnames(XX) <- NULL
        d <- darg(NULL, Xtrain)
        lagppreds <- aGPsep(X=Xtrain, Z=Ytrain, XX=XX, omp.threads=16, verb=0,
         end=25, method="nn", d=d)
        toc <- proc.time()[3]
        pred_times[m,i,2] <- toc-tic
        print("Finished laGP fit and predictions")
      }

      if (!too_long[3]) {
        tic <- proc.time()[3]
        Xtrain <- as.matrix(read.csv(xtrain_fn))
        colnames(Xtrain) <- NULL
        Ytrain <- drop(as.matrix(read.csv(ytrain_fn)))
        dgp1fit <- fit_one_layer(x=Xtrain, y=Ytrain, nmcmc=1000, vecchia=TRUE, m=10)
        toc <- proc.time()[3]
        fit_times[m,i,3] <- toc-tic
        print("Finished deep gp fit")

        tic <- proc.time()[3]
        XX <- as.matrix(read.csv(xtest_fn))
        colnames(XX) <- NULL
        dgp1preds <- predict(trim(dgp1fit, burn=100, thin=5), x_new=XX)
        toc <- proc.time()[3]
        pred_times[m,i,3] <- toc-tic
        print("Finished deep gp predictions")
      }
      print(paste0("Finished MC iteration ", m))
      res <- list(fit_times=fit_times, pred_times=pred_times)
      saveRDS(res, paste0(outf, format(Sys.time(), "%Y%m%d"), ".rds"))
    }
    too_long[1] <- ifelse(too_long[1], too_long[1], mean(pred_times[,i,1], na.rm=TRUE) >= 2700)
    too_long[2] <- ifelse(too_long[2], too_long[2], mean(pred_times[,i,2], na.rm=TRUE) >= 2700)
    too_long[3] <- ifelse(too_long[3], too_long[3], mean(pred_times[,i,3], na.rm=TRUE) >= 2700)
  }
} else {

  ns <- c(seq(10, 100, by=10), seq(500, 2500, by=500))
  mcs <- 10
  outf <- "surrogate_time_test_2dlogit_inc_runs_"

  ## DEFINE THE OUTPUT SIZE (10,000)
  nout <- 10000
  XX <- matrix(NA, nrow=nout, ncol=6)
  x <- seq(0, 1, length=sqrt(nout))
  X <- as.matrix(expand.grid(x, x))
  XX[,1:2] <- X
  XX[,3:6] <- rep(Utrueunit, each=nrow(XX))
  colnames(XX) <- NULL

  xtest_fn <- "xtest_nruns.csv"
  write.csv(XX, xtest_fn, row.names=FALSE)

  ## Calculating metrics
  fit_times <- pred_times <- array(NA, dim=c(mcs, length(ns), 3))
  too_long <- c(FALSE, TRUE, TRUE)#rep(FALSE, 3)

  ## LOOP THROUGH SIZES OF RUNS (10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 500, 1000)
  for (i in 1:length(ns)) {
    n <- ns[i]
    cat("n = ", n, "\n")

    #### LOOP THROUGH MONTE CARLOS
    for (m in 1:mcs) {

      xtrain_fn <- paste0("xtrain_nruns", n, "_mc", m, ".csv")
      ytrain_fn <- paste0("ytrain_nruns", n, "_mc", m, ".csv")

      ###### Create new U matrix
      U <- randomLHS(n=n, k=4)
      Uunit <- U
      U[,1] <- U[,1] * (mu_max - mu_min) + mu_min
      U[,2] <- U[,2] * (nu_max - nu_min) + nu_min
      U[,3] <- U[,3] * (mu_max - mu_min) + mu_min
      U[,4] <- U[,4] * (nu_max - nu_min) + nu_min

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

      write.csv(Xtrain, xtrain_fn, row.names=FALSE)
      write.csv(Ytrain, ytrain_fn, row.names=FALSE)
      Xtrain <- Ytrain <- NULL

      ###### Run fits
      if (!too_long[1]) {
        tic <- proc.time()[3]
        Xtrain <- as.matrix(read.csv(xtrain_fn))
        colnames(Xtrain) <- NULL
        Ytrain <- drop(as.matrix(read.csv(ytrain_fn)))
        svecfit <- fit_scaled(y=Ytrain, inputs=Xtrain, nug=1e-4, ms=25)
        toc <- proc.time()[3]
        fit_times[m,i,1] <- toc-tic
        print("Finished SVecchia fit")

        tic <- proc.time()[3]
        XX <- as.matrix(read.csv(xtest_fn))
        colnames(XX) <- NULL
        svecpreds <- predictions_scaled(svecfit, XX, m=25, joint=FALSE, predvar=TRUE)
        toc <- proc.time()[3]
        pred_times[m,i,1] <- toc-tic
        print("Finished SVecchia predictions")
      }

      if (!too_long[2]) {
        tic <- proc.time()[3]
        Xtrain <- as.matrix(read.csv(xtrain_fn))
        colnames(Xtrain) <- NULL
        Ytrain <- drop(as.matrix(read.csv(ytrain_fn)))
        XX <- as.matrix(read.csv(xtest_fn))
        colnames(XX) <- NULL
        d <- darg(NULL, Xtrain)
        lagppreds <- aGPsep(X=Xtrain, Z=Ytrain, XX=XX, omp.threads=16, verb=0,
         end=25, method="nn", d=d)
        toc <- proc.time()[3]
        pred_times[m,i,2] <- toc-tic
        print("Finished laGP fit and predictions")
      }

      if (!too_long[3]) {
        tic <- proc.time()[3]
        Xtrain <- as.matrix(read.csv(xtrain_fn))
        colnames(Xtrain) <- NULL
        Ytrain <- drop(as.matrix(read.csv(ytrain_fn)))
        dgp1fit <- fit_one_layer(x=Xtrain, y=Ytrain, nmcmc=1000, vecchia=TRUE, m=10)
        toc <- proc.time()[3]
        fit_times[m,i,3] <- toc-tic
        print("Finished deep gp fit")
        tic <- proc.time()[3]
        XX <- as.matrix(read.csv(xtest_fn))
        colnames(XX) <- NULL
        dgp1preds <- predict(trim(dgp1fit, burn=100, thin=5), x_new=XX)
        toc <- proc.time()[3]
        pred_times[m,i,3] <- toc-tic
        print("Finished deep gp predictions")
      }
      print(paste0("Finished MC iteration ", m))
      res <- list(fit_times=fit_times, pred_times=pred_times)
      saveRDS(res, paste0(outf, format(Sys.time(), "%Y%m%d"), ".rds"))
    }
    too_long[1] <- ifelse(too_long[1], too_long[1], mean(pred_times[,i,1], na.rm=TRUE) >= 2700)
    too_long[2] <- ifelse(too_long[2], too_long[2], mean(pred_times[,i,2], na.rm=TRUE) >= 2700)
    too_long[3] <- ifelse(too_long[3], too_long[3], mean(pred_times[,i,3], na.rm=TRUE) >= 2700)
  }
}
