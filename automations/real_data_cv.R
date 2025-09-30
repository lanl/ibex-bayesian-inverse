###############################################################################
#### Cross-validation code for Bayesian Computer Model Calibration for data
#### collected by the Interstellar Boundary Explorer satellite.
#### Uses field data from 2009-2011. Results will be used to calculate CRPS on
#### held out data.
###############################################################################

source("../helper.R")
source("../mcmc.R")

seed <- 711930
fold_seed <- 12937120
ncvs <- 10
fid <- 1
nmcmcs <- 10000

## read in the command line arguments
## run with: R CMD BATCH '--args seed=711930 fold_seed=12937120 ncvs=10 fid=1 nmcmcs=10000' real_data_cv.R
args <- commandArgs(TRUE)
if (length(args) > 0) {
  for(i in 1:length(args)) {
    eval(parse(text=args[[i]]))
  }
}
print(args)

model_data <- read.csv(file="../data/sims.csv")
field_data <- read.csv(file="../data/ibex_real.csv")
field_data <- field_data[field_data$map %in% paste0(2009:2011, "A"),]
nf <- nrow(field_data)
colnames(field_data)[colnames(field_data)== "counts"] <- "sim_counts"
colnames(field_data)[colnames(field_data)== "ecliptic_lat"] <- "lat"
colnames(field_data)[colnames(field_data)== "ecliptic_lon"] <- "lon"
field_data$xcod <- (field_data$x - min(field_data$x)) / diff(range(field_data$x))
field_data$ycod <- (field_data$y - min(field_data$y)) / diff(range(field_data$y))
field_data$zcod <- (field_data$z - min(field_data$z)) / diff(range(field_data$z))

set.seed(seed)
field_cv_inds <- sample(1:nf, nf)
field_cv_inds <- split(field_cv_inds, cut(seq_along(field_cv_inds), ncvs,
  labels=FALSE))
field_data_fold <- field_data[-field_cv_inds[[fid]],]

set.seed(fold_seed)
fparams <- paste0(2009:2011, "A")
pd <- preprocess_data(md=model_data, fd=field_data_fold, esa_lev=4,
  fparams=fparams, scales=c(1, 1), tol=NA, quant=0.0, real=TRUE)
svecfit <- fit_scaled(y=pd$Zmod, inputs=as.matrix(cbind(pd$Xmod, pd$Umod)),
  nug=1e-4, ms=25)
mcmc_res <- mcmc(Xm=pd$Xmod, Um=pd$Umod, Zm=pd$Zmod, Xf=pd$Xfield,
  Zf=pd$Zfield, Of=pd$Ofield, m=25, nmcmcs=nmcmcs, step=0.05,
  gpmeth="svecchia", vb=TRUE, true_u=NA, true_logscl=NA,
  betashape=2, adapt=FALSE)
save.image("real_data_cv.RData")

sample_pmfps <- seq(500, 3000, length=10)
sample_pmfps_cod <- (sample_pmfps - 500)/2500
sample_ratios <- seq(0.001, 0.1, length=10)
sample_ratios_cod <- (sample_ratios - 0.001)/(0.1-0.001)
sample_grid <- as.matrix(expand.grid(seq(500, 3000, length=40),
  seq(0.001, 0.1, length=40)))
sample_grid_cod <- sample_grid
sample_grid_cod[,1] <- (sample_grid_cod[,1] - 500)/2500
sample_grid_cod[,2] <- (sample_grid_cod[,2] - 0.001)/(0.1-0.001)

crps_pmfp <- rep(NA, length(sample_pmfps)
crps_ratio <- rep(NA, length(sample_ratios)
crps_grid <- rep(NA, nrow(sample_grid)

field_data <- read.csv(file="../data/ibex_real.csv")

### Read in results
iter_res <- readRDS(paste0("../results/", files[i]))
field_fold <- field_data[field_cv_inds,]
Xfold <- field_fold[,c("xcod", "ycod", "zcod")]

### Calculate posterior mean
post_mean_cod <- apply(mcmc_res$u[seq(1001, nmcmcs, by=10),], 2, mean)
post_mean <- post_mean_cod
post_mean[1] <- post_mean_cod[1]*2500 + 500
post_mean[2] <- post_mean_cod[2]*(0.1-0.001) + 0.001

for (j in 1:length(sample_pmfps_cod)) {
  ### Pull out of sample counts
  Xfield <- cbind(Xfold, matrix(c(sample_pmfps_cod[j], post_mean_cod[2]), nrow=1))
  colnames(Xfield) <- c("x", "y", "z", "pmfp", "ratio")

  ### Predict lambda at out of sample using Scaled Vecchia
  preds <- predictions_scaled(svecfit, as.matrix(Xfield), m=25,
    joint=FALSE, predvar=FALSE)
  ### Calculate CRPS (include time + background)
  crps_pmfp[j] <- mean(crps_pois(y=field_fold$sim_counts,
    lambda=preds*field_fold$time + field_fold$background))
  if (j %% 20 == 0) {
    print(paste0("Fold ", fid, ": Finished pmfp ", j, "/", length(sample_pmfps)))
  }
}
save.image("real_data_cv.RData")

for (j in 1:length(sample_ratios_cod)) {
  ### Pull out of sample counts
  Xfield <- cbind(Xfold, matrix(c(post_mean_cod[1], sample_ratios_cod[j]), nrow=1))
  colnames(Xfield) <- c("x", "y", "z", "pmfp", "ratio")

  ### Predict lambda at out of sample using Scaled Vecchia
  preds <- predictions_scaled(svecfit, as.matrix(Xfield), m=25,
    joint=FALSE, predvar=FALSE)

  ### Calculate CRPS (include time + background)
  crps_ratio[j] <- mean(crps_pois(y=field_fold$sim_counts,
    lambda=preds*field_fold$time + field_fold$background))
  if (j %% 20 == 0) {
    print(paste0("Fold ", fid, ": Finished ratio ", j, "/", length(sample_ratios)))
  }
}
save.image("real_data_cv.RData")

for (j in 1:nrow(sample_grid_cod)) {
  ### Pull out of sample counts
  Xfield <- cbind(Xfold, matrix(c(sample_grid_cod[j,1], sample_grid_cod[j,2]), nrow=1))
  colnames(Xfield) <- c("x", "y", "z", "pmfp", "ratio")

  ### Predict lambda at out of sample using Scaled Vecchia
  preds <- predictions_scaled(svecfit, as.matrix(Xfield), m=25,
    joint=FALSE, predvar=FALSE)

  ### Calculate CRPS (include time + background)
  crps_grid[j] <- mean(crps_pois(y=field_fold$sim_counts,
    lambda=preds*field_fold$time + field_fold$background))
  if (j %% 20 == 0) {
    print(paste0("Fold ", fid, ": Finished combo ", j, "/", nrow(sample_grid)))
  }
}
save.image("real_data_cv.RData")
print(paste0("Fold ", fid, ": COMPLETE"))

res <- list(mcmc_res=mcmc_res, crps_pmfp=crps_pmfp, crps_ratio=crps_ratio,
  crps_grid=crps_grid, pmfp_grid=sample_pmfps, ratio_grid=sample_ratios,
  grid=sample_grid, post_mean=post_mean,
  settings=list(seed=seed, ncvs=ncvs, fid=fid))

saveRDS(res, file=paste0("../results/real_data_cv_fold", fid, "_seed", seed,
  format(Sys.time(), "_%Y%m%d%H%M%S"), ".rds"))
