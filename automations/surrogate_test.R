library(deepgp)
library(laGP)

source("../helper.R")
source('../vecchia_scaled.R')

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

model_data <- read.csv(file="../data/sims.csv")
model_data <- model_data[model_data$ESA==4,]
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

## Calculating metrics
rmses <- crps <- fit_times <- pred_times <-
  matrix(NA, ncol=3, nrow=nrow(unique_runs))
colnames(rmses) <- colnames(crps) <- colnames(fit_times) <-
  colnames(pred_times) <- c("svecchia", "lagp", "deepgp")

if (method=="svecchia" || method=="all") {
  for (i in start:nrow(unique_runs)) {
    pmfp <- unique_runs[i,1]
    ratio <- unique_runs[i,2]
    Xtrain <- model_data[model_data$parallel_mean_free_path != pmfp |
      model_data$ratio != ratio,c("parallel_mean_free_path", "ratio", "x", "y", "z")]
    Ytrain <- model_data[model_data$parallel_mean_free_path != pmfp |
      model_data$ratio != ratio,c("blurred_ena_rate")]
    Xtest <- model_data[model_data$parallel_mean_free_path == pmfp &
      model_data$ratio == ratio,c("parallel_mean_free_path", "ratio", "x", "y", "z")]
    Ytest <- model_data[model_data$parallel_mean_free_path == pmfp &
      model_data$ratio == ratio,c("blurred_ena_rate")]

    tic <- proc.time()[3]
    svecfit <- fit_scaled(y=Ytrain, inputs=as.matrix(Xtrain), nug=1e-4, ms=25)
    toc <- proc.time()[3]
    fit_times[i,1] <- toc-tic
    print("Finished SVecchia fit")

    tic <- proc.time()[3]
    svecpreds <- predictions_scaled(svecfit, as.matrix(Xtest), m=25, joint=FALSE,
      predvar=TRUE)
    toc <- proc.time()[3]
    pred_times[i,1] <- toc-tic
    rmses[i,1] <- sqrt(mean((svecpreds$means - Ytest)^2))
    crps[i,1] <- crps(y=Ytest, mu=svecpreds$means, s2=svecpreds$vars)
    print("Finished SVecchia predictions")
    print(paste0("Finished holdout iteration ", i))
    res <- list(fit_times=fit_times, pred_times=pred_times, rmse=rmses, crps=crps)
    saveRDS(res, paste0("surrogate_test_", format(Sys.time(), "_%Y%m%d"), ".rds"))
  }
}

if (method=="laGP" || method=="all") {
  for (i in start:nrow(unique_runs)) {
    pmfp <- unique_runs[i,1]
    ratio <- unique_runs[i,2]
    Xtrain <- model_data[model_data$parallel_mean_free_path != pmfp |
      model_data$ratio != ratio,c("parallel_mean_free_path", "ratio", "x", "y", "z")]
    Ytrain <- model_data[model_data$parallel_mean_free_path != pmfp |
      model_data$ratio != ratio,c("blurred_ena_rate")]
    Xtest <- model_data[model_data$parallel_mean_free_path == pmfp &
      model_data$ratio == ratio,c("parallel_mean_free_path", "ratio", "x", "y", "z")]
    Ytest <- model_data[model_data$parallel_mean_free_path == pmfp &
      model_data$ratio == ratio,c("blurred_ena_rate")]

    tic <- proc.time()[3]
    d <- darg(NULL, cbind(Xm, Um))
    lagppreds <- aGPsep(X=Xtrain, Z=Ytrain, XX=Xtest, omp.threads=1, verb=0,
      end=25, method="nn", d=d)
    toc <- proc.time()[3]
    pred_times[i,2] <- toc-tic
    rmses[i,2] <- sqrt(mean((lagppreds$mean - Ytest)^2))
    crps[i,2] <- crps(y=Ytest, mu=lagppreds$mean, s2=lagppreds$var)
    print("Finished laGP fit and predictions")
    print(paste0("Finished holdout iteration ", i))
    res <- list(fit_times=fit_times, pred_times=pred_times, rmse=rmses, crps=crps)
    saveRDS(res, paste0("surrogate_test_", format(Sys.time(), "_%Y%m%d"), ".rds"))
  }
}

if (method=="deepgp" || method=="all") {
  for (i in start:nrow(unique_runs)) {
    pmfp <- unique_runs[i,1]
    ratio <- unique_runs[i,2]
    Xtrain <- model_data[model_data$parallel_mean_free_path != pmfp |
      model_data$ratio != ratio,c("parallel_mean_free_path", "ratio", "x", "y", "z")]
    Ytrain <- model_data[model_data$parallel_mean_free_path != pmfp |
      model_data$ratio != ratio,c("blurred_ena_rate")]
    Xtest <- model_data[model_data$parallel_mean_free_path == pmfp &
      model_data$ratio == ratio,c("parallel_mean_free_path", "ratio", "x", "y", "z")]
    Ytest <- model_data[model_data$parallel_mean_free_path == pmfp &
      model_data$ratio == ratio,c("blurred_ena_rate")]

    tic <- proc.time()[3]
    dgp1fit <- fit_one_layer(x=as.matrix(Xtrain), y=Ytrain, nmcmc=1000, vecchia=TRUE, m=10)
    toc <- proc.time()[3]
    fit_times[i,3] <- toc-tic
    print("Finished deep gp fit")

    tic <- proc.time()[3]
    dgp1preds <- predict(trim(dgp1fit, burn=100, thin=5), x_new=Xtest)
    toc <- proc.time()[3]
    pred_times[i,3] <- toc-tic
    rmses[i,3] <- sqrt(mean((dgp1preds$mean - Ytest)^2))
    crps[i,3] <- crps(y=Ytest, mu=dgp1preds$mean, s2=dgp1preds$s2)
    print("Finished deep gp predictions")
    print(paste0("Finished holdout iteration ", i))
    res <- list(fit_times=fit_times, pred_times=pred_times, rmse=rmses, crps=crps)
    saveRDS(res, paste0("surrogate_test_", format(Sys.time(), "_%Y%m%d"), ".rds"))
  }
}

res <- list(fit_times=fit_times, pred_times=pred_times, rmse=rmses, crps=crps)
saveRDS(res, paste0("surrogate_test_", format(Sys.time(), "_%Y%m%d"), ".rds"))
