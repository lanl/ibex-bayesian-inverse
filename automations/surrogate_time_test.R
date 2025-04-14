library(deepgp)
library(dplyr)
library(laGP)

source("../helper.R")
source('../vecchia_scaled.R')

seed <- 781691
large_n <- 0

## read in the command line arguments
## run with: R CMD BATCH '--args seed=1 large_n=0' surrogate_time_test.R
args <- commandArgs(TRUE)
if (length(args) > 0) {
  for(i in 1:length(args)) {
    eval(parse(text=args[[i]]))
  }
}

model_data <- read.csv(file="../data/sims.csv")
model_data <- model_data[model_data$ESA==4,]
model_data[,c("x", "y", "z")] <- geo_to_spher_coords(lat=model_data$lat,
  lon=model_data$lon)
model_data <- model_data[,c("lat", "lon", "x", "y", "z", "parallel_mean_free_path",
 "ratio", "blurred_ena_rate")]

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

sim_runs <- list()
for (i in 1:nrow(calib_grid)) {
  pmfp <- calib_grid[i,1]
  rat <- calib_grid[i,2]
  sim_runs[[i]] <- model_data[model_data$parallel_mean_free_path==pmfp &
    model_data$ratio==rat,]
}

## Calculating metrics
set.seed(seed)
exp_pows <- 7:45
ns <- seq(20000, 75000, by=5000)
num_ns <- ifelse(large_n, length(ns), length(exp_pows))
outf <- paste0("surrogate_time_test_", ifelse(large_n, "large_", ""))

mcs <- 5
fit_times <- pred_times <- array(NA, dim=c(5, length(exp_pows), 3))
too_long <- rep(FALSE, 3)

for (i in 1:length(exp_pows)) {

  ## select training data size
  if (large_n) {
    n <- ns[i]
    cat("n = ", n)
  } else {
    n <- round(10 + 1.25^exp_pows[i])
    cat("n = 1.25^", exp_pows[i], " = ", n, "\n", sep="")
  }

  for (m in 1:mcs) {

    sim_run <- sample(1:nruns, 1)
    inds <- sample(1:nresponses, n, replace=(n > nresponses))

    Xtest <- sim_runs[[i]][inds,c("parallel_mean_free_path", "ratio", "x", "y", "z")]
    Ytest <- sim_runs[[i]][inds,c("blurred_ena_rate")]

    Xtrain <- data.frame(matrix(NA, nrow=0, ncol=5))
    Ytrain <- matrix(NA, nrow=0, ncol=1)
    train_runs <- (1:nruns)[-sim_run]

    for (j in train_runs) {
      Xtrain <- rbind(Xtrain, sim_runs[[j]][inds,c("parallel_mean_free_path", "ratio", "x", "y", "z")])
      Ytrain <- rbind(Ytrain, sim_runs[[j]][inds,c("blurred_ena_rate"),drop=FALSE])
    }
    Ytrain <- Ytrain$blurred_ena_rate

    write.csv(Xtrain, "xtrain_iter.csv", row.names=FALSE)
    write.csv(Ytrain, "ytrain_iter.csv", row.names=FALSE)

    if (!too_long[1]) {
      tic <- proc.time()[3]
      Xtrain_iter <- read.csv("xtrain_iter.csv")
      Ytrain_iter <- read.csv("ytrain_iter.csv")[,1]
      svecfit <- fit_scaled(y=c(Ytrain_iter), inputs=as.matrix(Xtrain_iter),
        nug=1e-4, ms=25)
      toc <- proc.time()[3]
      fit_times[m,i,1] <- toc-tic
      print("Finished SVecchia fit")

      tic <- proc.time()[3]
      svecpreds <- predictions_scaled(svecfit, as.matrix(Xtest), m=25, joint=FALSE,
        predvar=TRUE)
      toc <- proc.time()[3]
      pred_times[m,i,1] <- toc-tic
      print("Finished SVecchia predictions")
    }

    if (!too_long[2]) {
      tic <- proc.time()[3]
      Xtrain_iter <- read.csv("xtrain_iter.csv")
      Ytrain_iter <- read.csv("ytrain_iter.csv")[,1]
      d <- darg(NULL, Xtrain_iter)
      lagppreds <- aGPsep(X=Xtrain_iter, Z=Ytrain_iter, XX=Xtest, omp.threads=16,
       verb=0, end=25, method="nn", d=d)
      toc <- proc.time()[3]
      pred_times[m,i,2] <- toc-tic
      print("Finished laGP fit and predictions")
    }

    if (!too_long[3]) {
      tic <- proc.time()[3]
      Xtrain_iter <- read.csv("xtrain_iter.csv")
      Ytrain_iter <- read.csv("ytrain_iter.csv")[,1]
      dgp1fit <- fit_one_layer(x=as.matrix(Xtrain_iter), y=Ytrain_iter, nmcmc=1000,
       vecchia=TRUE, m=10)
      toc <- proc.time()[3]
      fit_times[m,i,3] <- toc-tic
      print("Finished deep gp fit")

      tic <- proc.time()[3]
      dgp1preds <- predict(trim(dgp1fit, burn=100, thin=5), x_new=Xtest)
      toc <- proc.time()[3]
      pred_times[m,i,3] <- toc-tic
      print("Finished deep gp predictions")

      file.remove("xtrain_iter.csv")
      file.remove("ytrain_iter.csv")
    }
    res <- list(fit_times=fit_times, pred_times=pred_times)
    saveRDS(res, paste0(outf, format(Sys.time(), "%Y%m%d"), ".rds"))
  }
  too_long <- apply(pred_times[,i,], 2, mean) >= 3600
}

res <- list(fit_times=fit_times, pred_times=pred_times)
saveRDS(res, paste0(outf, format(Sys.time(), "%Y%m%d"), ".rds"))
