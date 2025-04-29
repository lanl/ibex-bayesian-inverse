library(deepgp)
library(laGP)

source("../helper.R")
source('../vecchia_scaled.R')

seed <- 679818781
start <- 1

## read in the command line arguments
## run with: R CMD BATCH '--args seed=1 start=1' surrogate_test.R
args <- commandArgs(TRUE)
if (length(args) > 0) {
  for(i in 1:length(args)) {
    eval(parse(text=args[[i]]))
  }
}
print(args)
settings <- list(seed=seed, start=start)
print(settings)

model_data <- read.csv(file="../data/sims.csv")
model_data <- model_data[order(model_data$parallel_mean_free_path, model_data$ratio,
  model_data$lat, model_data$lon),]
model_data[,c("x", "y", "z")] <- geo_to_spher_coords(lat=model_data$lat,
  lon=model_data$lon)
model_data <- model_data[,c("x", "y", "z", "parallel_mean_free_path", "ratio",
 "blurred_ena_rate")]

for (i in 1:(ncol(model_data)-1)) {
  minval <- min(model_data[,i])
  maxval <- max(model_data[,i])
  valrange <- diff(range(model_data[,i]))
  model_data[,i] <- (model_data[,i] - minval)/(valrange)
}

## Hold one out - 66 model runs
pmfps <- unique(model_data$parallel_mean_free_path)
ratios <- unique(model_data$ratio)
unique_runs <- expand.grid(pmfps, ratios)
mfs <- seq(5, 100, by=5)
mps <- seq(5, 100, by=5)

## Calculating metrics
rmses <- crps <- fit_times <- pred_times <-
  matrix(NA, ncol=length(mfs), nrow=nrow(unique_runs))
colnames(rmses) <- colnames(crps) <- colnames(fit_times) <-
  colnames(pred_times) <- paste0("m", mfs)

for (i in 1:length(mfs)) {
  mf <- mfs[i]
  mp <- mps[i]

  for (j in start:nrow(unique_runs)) {
    pmfp <- unique_runs[j,1]
    ratio <- unique_runs[j,2]
    Xtrain <- model_data[model_data$parallel_mean_free_path != pmfp |
      model_data$ratio != ratio,c("parallel_mean_free_path", "ratio", "x", "y", "z")]
    Ytrain <- model_data[model_data$parallel_mean_free_path != pmfp |
      model_data$ratio != ratio,c("blurred_ena_rate")]
    Xtest <- model_data[model_data$parallel_mean_free_path == pmfp &
      model_data$ratio == ratio,c("parallel_mean_free_path", "ratio", "x", "y", "z")]
    Ytest <- model_data[model_data$parallel_mean_free_path == pmfp &
      model_data$ratio == ratio,c("blurred_ena_rate")]

    tic <- proc.time()[3]
    svecfit <- fit_scaled(y=Ytrain, inputs=as.matrix(Xtrain), nug=1e-4, ms=mf)
    toc <- proc.time()[3]
    fit_times[j,i] <- toc-tic
    print(toc-tic)

    tic <- proc.time()[3]
    svecpreds <- predictions_scaled(svecfit, as.matrix(Xtest), m=mp, joint=FALSE,
      predvar=TRUE)
    toc <- proc.time()[3]
    print(toc-tic)
    pred_times[j,i] <- toc-tic
    rmses[j,i] <- sqrt(mean((svecpreds$means - Ytest)^2))
    crps[j,i] <- crps(y=Ytest, mu=svecpreds$means, s2=svecpreds$vars)
    print(paste0("Finished holdout iteration ", j))
    res <- list(fit_times=fit_times, pred_times=pred_times, rmse=rmses, crps=crps)
    saveRDS(res, paste0("svecchia_m_test_", format(Sys.time(), "%Y%m%d"), ".rds"))
  }
  print(paste0("Finished SVecchia iteration m = ", mf))
}

res <- list(fit_times=fit_times, pred_times=pred_times, rmse=rmses, crps=crps)
saveRDS(res, paste0("svecchia_m_test_", format(Sys.time(), "%Y%m%d"), ".rds"))
