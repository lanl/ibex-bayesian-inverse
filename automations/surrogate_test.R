###############################################################################
#### Hold one out testing of different surrogate methods: Scaled Vecchia,
#### laGP, and deepgp with one layer. Testing data is from the IBEX
#### simulator. Each of 66 unique simulator runs is used as a holdout set once.
####
#### DATA NEEDED: sims.csv
###############################################################################

library(deepgp)
library(laGP)

setwd("..")
source("helper.R")
source("vecchia_scaled.R")
setwd("automations")

start <- 1
method <- "all"
seed <- 6756781

## read in the command line arguments
## run with: R CMD BATCH '--args seed=1 method=all start=1' surrogate_test.R
args <- commandArgs(TRUE)
if (length(args) > 0) {
  for(i in 1:length(args)) {
    eval(parse(text=args[[i]]))
  }
}
settings <- list(seed=seed, method=method, start=start)
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
unique_runs <- as.matrix(expand.grid(pmfps, ratios))
colnames(unique_runs) <- NULL

## Calculating metrics
rmses <- crps <- fit_times <- pred_times <-
  matrix(NA, ncol=6, nrow=nrow(unique_runs))
colnames(rmses) <- colnames(crps) <- colnames(fit_times) <-
  colnames(pred_times) <- c("svecchia25", "svecchia50", "svecchia75",
   "svecchia100", "lagp", "deepgp")

if (method=="svecchia" || method=="all") {
  ms <- seq(25, 100, by=25)
  for (i in 1:length(ms))) {
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
      svecfit <- fit_scaled(y=Ytrain, inputs=as.matrix(Xtrain), nug=1e-4, ms=ms[i])
      toc <- proc.time()[3]
      fit_times[j,i] <- toc-tic
      print(paste0("Finished SVecchia fit where m=", ms[i]))

      tic <- proc.time()[3]
      svecpreds <- predictions_scaled(svecfit, as.matrix(Xtest), m=ms[i], joint=FALSE,
        predvar=TRUE)
      toc <- proc.time()[3]
      pred_times[j,i] <- toc-tic
      rmses[j,i] <- sqrt(mean((svecpreds$means - Ytest)^2))
      crps[j,i] <- crps(y=Ytest, mu=svecpreds$means, s2=svecpreds$vars)
      print(paste0("Finished SVecchia predictions where m=", ms[i]))
      print(paste0("Finished holdout iteration ", j))
      res <- list(fit_times=fit_times, pred_times=pred_times, rmse=rmses, crps=crps)
      saveRDS(res, paste0("surrogate_test_", format(Sys.time(), "%Y%m%d"), ".rds"))
    }
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

    lagppreds <- NULL
    attempts <- 1
    while (is.null(lagppreds) && attempts <= 10) {
      try({  ## occasionally laGP returns a Cholesky decomposition error
        print(paste0("Making attempt #", attempts))
        attempts <- attempts + 1
        tic <- proc.time()[3]
        d <- darg(NULL, Xtrain)
        lagppreds <- aGPsep(X=Xtrain, Z=Ytrain, XX=Xtest, omp.threads=8, verb=0,
          end=50, method="nn", d=d)
        toc <- proc.time()[3]
      }, silent=TRUE)
    }
    pred_times[i,5] <- toc-tic
    rmses[i,5] <- sqrt(mean((lagppreds$mean - Ytest)^2))
    crps[i,5] <- crps(y=Ytest, mu=lagppreds$mean, s2=lagppreds$var)
    print("Finished laGP fit and predictions")
    print(paste0("Finished holdout iteration ", i))
    res <- list(fit_times=fit_times, pred_times=pred_times, rmse=rmses, crps=crps)
    saveRDS(res, paste0("surrogate_test_", format(Sys.time(), "%Y%m%d"), ".rds"))
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
    fit_times[i,6] <- toc-tic
    print("Finished deep gp fit")

    tic <- proc.time()[3]
    dgp1preds <- predict(trim(dgp1fit, burn=100, thin=5), x_new=Xtest)
    toc <- proc.time()[3]
    pred_times[i,6] <- toc-tic
    rmses[i,6] <- sqrt(mean((dgp1preds$mean - Ytest)^2))
    crps[i,6] <- crps(y=Ytest, mu=dgp1preds$mean, s2=dgp1preds$s2)
    print("Finished deep gp predictions")
    print(paste0("Finished holdout iteration ", i))
    res <- list(fit_times=fit_times, pred_times=pred_times, rmse=rmses, crps=crps)
    saveRDS(res, paste0("surrogate_test_", format(Sys.time(), "%Y%m%d"), ".rds"))
  }
}

res <- list(fit_times=fit_times, pred_times=pred_times, rmse=rmses, crps=crps)
saveRDS(res, paste0("surrogate_test_", format(Sys.time(), "%Y%m%d"), ".rds"))
