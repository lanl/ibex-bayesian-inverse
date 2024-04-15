library(doParallel)
library(parallel)
library(tidyverse)

source("../helper.R")
source("../mcmc.R")

args <- R.utils::commandArgs(asValues=TRUE)

esa_lev <- ifelse(!is.null(args$esa), as.integer(args$esa), 3)
num_steps <- ifelse(!is.null(args$steps), as.integer(args$steps), 13)
nmcs <- ifelse(!is.null(args$nmcs), as.integer(args$nmcs), 50)
nmcmcs <- ifelse(!is.null(args$nmcmcs), as.integer(args$nmcmcs), 100)
target_rate <- ifelse(!is.null(args$rate), as.numeric(args$rate), 0.239)

model_data <- readRDS(file="../data/Simulations.rds")
field_data <- readRDS(file="../data/simulated_binned_direct_events_data.rds")

steps <- 0.01*2^seq(-7, 5, length=num_steps)
cpars <- field_data %>% dplyr::distinct(parallel_mean_free_path, ratio)
rates <- data.frame(matrix(data=NA, nrow=nmcs*length(steps), ncol=2))
colnames(rates) <- c("step_size", "acceptance_rate")
row <- 1

tic <- proc.time()[3]
for (s in steps) {
  cl <- parallel::makeCluster(10)
  doParallel::registerDoParallel(cl)
  step_rates <- foreach(i = 1:nmcs, .combine='c', .packages=c("GPvecchia",
      "GpGp", "laGP", "mvtnorm", "tmvtnorm", "tidyverse")) %dopar% {
    ind <- sample(1:nrow(cpars), 1)
    fpmfp <- cpars[ind,]$parallel_mean_free_path
    fratio <- cpars[ind,]$ratio
    pd <- preprocess_data(md=model_data, fd=field_data, esa_lev=esa_lev,
                          fparams=c(fpmfp, fratio), scales=c(1, 1), quant=0)
    mcmc_res <- mcmc(Xm=pd$Xmod, Um=pd$Umod, Zm=pd$Zmod, Xf=pd$Xfield,
                     Zf=pd$Zfield, Of=pd$Ofield, nmcmcs=nmcmcs, step=s,
                     gpmeth="svecchia", end=25, thrds=1, vb=TRUE, debug=FALSE)
    length(unique(mcmc_res$u[,1]))/nmcmcs
  }
  parallel::stopCluster(cl)
  rates[row:(row+nmcs-1),1] <- s
  rates[row:(row+nmcs-1),2] <- step_rates
  row <- row + nmcs
  print(paste0("Finished step size: ", s))
  print(paste0("Elapsed time: ", proc.time()[3] - tic))
  saveRDS(rates, file=paste0("../results/tune_step_size_temp.rds"))
}
toc <- proc.time()[3]

rates$rate_logit <- log((rates$acceptance_rate)/(1-rates$acceptance_rate))
rates$step_log <- log(rates$step_size)
fit <- lm(rate_logit ~ step_log, data=rates)
s_hat <- exp((1/fit$coefficients[2])*(log(target_rate/(1-target_rate))-fit$coefficients[1]))
res <- list(rates=rates, steps=steps, step_size=s_hat, time=toc-tic)
saveRDS(res, file=paste0("../results/tuning_steps_test_",
                         format(Sys.time(), "%Y%m%d%H%M%S"), ".rds"))
