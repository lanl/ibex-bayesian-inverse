library(deepgp)
library(dplyr)
library(laGP)

source("../helper.R")
source('../vecchia_scaled.R')

seed <- 781691

## read in the command line arguments
## run with: R CMD BATCH '--args seed=1' surrogate_time_test.R
args <- commandArgs(TRUE)
if (length(args) > 0) {
  for(i in 1:length(args)) {
    eval(parse(text=args[[i]]))
  }
}

## Load in computer model data from file
model_data <- read.csv(file="../data/sims.csv")
model_data[,c("x", "y", "z")] <- geo_to_spher_coords(lat=model_data$lat,
  lon=model_data$lon)
model_data <- model_data[,c("lat", "lon", "x", "y", "z", "parallel_mean_free_path",
 "ratio", "blurred_ena_rate")]

## Scale every column but the response to 0-1
for (i in 1:(ncol(model_data)-1)) {
  minval <- min(model_data[,i])
  maxval <- max(model_data[,i])
  valrange <- diff(range(model_data[,i]))
  model_data[,i] <- (model_data[,i] - minval)/(valrange)
}

pmfps <- unique(model_data$parallel_mean_free_path)
rats <- unique(model_data$ratio)
calib_grid <- expand.grid(pmfps, rats)
colnames(calib_grid) <- c("pmfp", "ratio")
nruns <- nrow(calib_grid)
nresponses <- nrow(model_data) / nruns

## Separate the output of each simulation run into an element in a list
sim_runs <- list()
for (i in 1:nrow(calib_grid)) {
  pmfp <- calib_grid[i,1]
  rat <- calib_grid[i,2]
  sim_runs[[i]] <- model_data[model_data$parallel_mean_free_path==pmfp &
    model_data$ratio==rat,]
}

## Setting data sizes
set.seed(seed)
exp_pows <- 7:45
ns <- c(round(10 + 1.25^exp_pows), seq(20000, 75000, by=5000))
num_ns <- length(ns)

## Defining some parameters of the bakeoff
outf <- "surrogate_time_test_"
mcs <- 5
fit_times <- pred_times <- array(NA, dim=c(5, length(exp_pows), 6))
too_long <- rep(FALSE, 6)

for (i in 1:num_ns) {
  ## select training data size
  n <- ns[i]
  if (i > length(exp_pows)) {
    cat("n = ", n)
  } else {
    cat("n = 1.25^", exp_pows[i], " = ", n, "\n", sep="")
  }

  for (m in 1:mcs) {

    ## Building test set, drawn from a random computer model run
    sim_run <- sample(1:nruns, 1)
    inds <- sample(1:nresponses, n, replace=(n > nresponses))
    Xtest <- sim_runs[[sim_run]][inds,c("parallel_mean_free_path", "ratio", "x", "y", "z")]
    Ytest <- sim_runs[[sim_run]][inds,c("blurred_ena_rate")]

    ## Building training set, taken from all non-test set runs. But at the same indices
    Xtrain_list <- Ytrain_list <- list()
    Xtrain <- data.frame(matrix(NA, nrow=0, ncol=5))
    Ytrain <- matrix(NA, nrow=0, ncol=1)
    train_runs <- (1:nruns)[-sim_run]
    for (j in 1:length(train_runs)) {
      iter_run <- train_runs[j]
      Xtrain_list[[j]] <- sim_runs[[iter_run]][inds,c("parallel_mean_free_path", "ratio", "x", "y", "z")]
      Ytrain_list[[j]] <- sim_runs[[iter_run]][inds,c("blurred_ena_rate"),drop=FALSE]
    }
    Xtrain <- do.call("rbind", Xtrain_list)
    Ytrain <- do.call("rbind", Ytrain_list)$blurred_ena_rate

    ## Save training sets so that each method needs to read in the file
    write.csv(Xtrain, "xtrain_iter.csv", row.names=FALSE)
    write.csv(Ytrain, "ytrain_iter.csv", row.names=FALSE)

    ms <- seq(25, 100, by=25)
    for (j in 1:length(ms)) {
      if (!too_long[j]) {
        tic <- proc.time()[3]
        Xtrain_iter <- read.csv("xtrain_iter.csv")
        Ytrain_iter <- read.csv("ytrain_iter.csv")[,1]
        svecfit <- fit_scaled(y=c(Ytrain_iter), inputs=as.matrix(Xtrain_iter),
          nug=1e-4, ms=ms[j])
        toc <- proc.time()[3]
        fit_times[m,i,j] <- toc-tic
        print(paste0("Finished SVecchia fit with m=", ms[j]))

        tic <- proc.time()[3]
        svecpreds <- predictions_scaled(svecfit, as.matrix(Xtest), m=ms[j], joint=FALSE,
          predvar=TRUE)
        toc <- proc.time()[3]
        pred_times[m,i,j] <- toc-tic
        print(paste0("Finished SVecchia predictions with m=", ms[j]))
      }
    }

    if (!too_long[5]) {
      tic <- proc.time()[3]
      Xtrain_iter <- read.csv("xtrain_iter.csv")
      Ytrain_iter <- read.csv("ytrain_iter.csv")[,1]
      d <- darg(NULL, Xtrain_iter)
      fit_times[m,i,5] <- 0
      lagppreds <- aGPsep(X=Xtrain_iter, Z=Ytrain_iter, XX=Xtest, omp.threads=16,
       verb=0, end=25, method="nn", d=d)
      toc <- proc.time()[3]
      pred_times[m,i,5] <- toc-tic
      print("Finished laGP fit and predictions")
    }

    if (!too_long[6]) {
      tic <- proc.time()[3]
      Xtrain_iter <- read.csv("xtrain_iter.csv")
      Ytrain_iter <- read.csv("ytrain_iter.csv")[,1]
      dgp1fit <- fit_one_layer(x=as.matrix(Xtrain_iter), y=Ytrain_iter, nmcmc=1000,
       vecchia=TRUE, m=10)
      toc <- proc.time()[3]
      fit_times[m,i,6] <- toc-tic
      print("Finished deep gp fit")

      tic <- proc.time()[3]
      dgp1preds <- predict(trim(dgp1fit, burn=100, thin=5), x_new=Xtest)
      toc <- proc.time()[3]
      pred_times[m,i,6] <- toc-tic
      print("Finished deep gp predictions")

      file.remove("xtrain_iter.csv")
      file.remove("ytrain_iter.csv")
    }
    res <- list(fit_times=fit_times, pred_times=pred_times)
    saveRDS(res, paste0(outf, format(Sys.time(), "%Y%m%d"), ".rds"))
  }
  too_long_fits <- apply(fit_times[,i,], 2, mean, na.rm=TRUE) >= 3600
  too_long_fits[is.na(too_long_fits)] <- TRUE
  too_long_preds <- apply(pred_times[,i,], 2, mean, na.rm=TRUE) >= 3600
  too_long_preds[is.na(too_long_preds)] <- TRUE
  too_long <- too_long | too_long_fits | too_long_preds
  print(paste0("Current status of each implementation after dimension size ", n, ":"))
  print(too_long)
}

res <- list(fit_times=fit_times, pred_times=pred_times)
saveRDS(res, paste0(outf, format(Sys.time(), "%Y%m%d"), ".rds"))
