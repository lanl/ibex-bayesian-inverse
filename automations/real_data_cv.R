###############################################################################
#### Cross-validation codee for solving a Bayesian inverse problem for counts
#### collected by the Interstellar Boundary Explorer satellite. Uses field data
#### from 2009-2011. Results will be used to calculate CRPS on held out data.
###############################################################################

library(scoringRules)

setwd("..")
source("pois_bayes_inv.R")
source("helper.R")
setwd("automations")

seed <- 711930
fold_seed <- fid <- 1
ncvs <- 10
nmcmcs <- 10000

## read in the command line arguments
## run with: R CMD BATCH '--args seed=711930 fold_seed=1 ncvs=10 fid=1 nmcmcs=10000' real_data_cv.R
args <- commandArgs(TRUE)
if (length(args) > 0) {
  for(i in 1:length(args)) {
    eval(parse(text=args[[i]]))
  }
}

# load in simulator and satellite data
model_data <- read.csv(file="../data/sims.csv")
field_data <- read.csv(file="../data/ibex_real.csv")

set.seed(fold_seed)
maps <- paste0(2009:2011, "A")
pd <- preprocess_data(md=model_data, fd=field_data, map=maps, esa_lev=4)
svecfit <- fit_scaled(y=pd$ym, inputs=cbind(pd$xm, pd$um), nug=1e-4, ms=25)
nf <- nrow(pd$xf)

set.seed(seed)
field_cv_inds <- sample(1:nf, nf)
field_cv_inds <- split(field_cv_inds, cut(seq_along(field_cv_inds), ncvs,
  labels=FALSE))

fold_xf <- pd$xf[-field_cv_inds[[fid]],]
fold_yf <- pd$yf[-field_cv_inds[[fid]]]
fold_e <- pd$e[-field_cv_inds[[fid]]]
fold_bg <- pd$bg[-field_cv_inds[[fid]]]

res <- pois_bayes_inv(xm=pd$xm, um=pd$um, ym=pd$ym, xf=fold_xf, yf=fold_yf,
  e=fold_e, fold_bg, T=nmcmcs)
save.image(paste0("real_data_cv_", fold_seed, "_", seed, ".RData"))

pmfps <- seq(500, 3000, length=200)
pmfps_unit <- (pmfps - min(pmfps))/diff(range(pmfps))
ratios <- seq(0.001, 0.1, length=200)
ratios_unit <- (ratios - min(ratios))/diff(range(ratios))
grid <- as.matrix(expand.grid(seq(500, 3000, length=30),
  seq(0.001, 0.1, length=30)))
grid_unit <- grid
grid_unit[,1] <- (grid_unit[,1] - min(pmfps))/diff(range(pmfps))
grid_unit[,2] <- (grid_unit[,2] - min(ratios))/diff(range(ratios))

crps_pmfp <- rep(NA, length(pmfps))
crps_ratio <- rep(NA, length(ratios))
crps_grid <- rep(NA, nrow(grid))

### Read in test data
fold_test_xf <- pd$xf[field_cv_inds[[fid]],]
fold_test_yf <- pd$yf[field_cv_inds[[fid]]]
fold_test_e <- pd$e[field_cv_inds[[fid]]]
fold_test_bg <- pd$bg[field_cv_inds[[fid]]]

### Calculate posterior mean
post_mean_unit <- apply(res$u[seq(nmcmcs*0.1, nmcmcs, by=10),], 2, mean)
post_mean <- post_mean_unit
post_mean[1] <- post_mean_unit[1]*diff(range(pmfps)) + min(pmfps)
post_mean[2] <- post_mean_unit[2]*diff(range(ratios)) + min(ratios)

for (j in 1:length(pmfps_unit)) {
  ### Pull out of sample counts
  xf <- cbind(fold_test_xf, matrix(c(pmfps_unit[j], post_mean_unit[2]),
    nrow=nrow(fold_test_xf), ncol=length(post_mean), byrow=TRUE))

  ### Predict lambda at out of sample using Scaled Vecchia
  preds <- predictions_scaled(svecfit, xf, m=25, joint=FALSE)
  ### Calculate CRPS (include time + background)
  crps_pmfp[j] <- mean(crps_pois(y=fold_test_yf,
    lambda=(preds+fold_test_bg)*fold_test_e))
  if (j %% 10 == 0) {
    print(paste0("Fold ", fid, ": Finished pmfp ", j, "/", length(pmfps)))
  }
}
save.image(paste0("real_data_cv_", fold_seed, "_", seed, ".RData"))

for (j in 1:length(ratios_unit)) {
  ### Pull out of sample counts
  xf <- cbind(fold_test_xf, matrix(c(post_mean_unit[1], ratios_unit[j]),
    nrow=nrow(fold_test_xf), ncol=length(post_mean), byrow=TRUE))

  ### Predict lambda at out of sample using Scaled Vecchia
  preds <- predictions_scaled(svecfit, xf, m=25, joint=FALSE)

  ### Calculate CRPS (include time + background)
  crps_ratio[j] <- mean(crps_pois(y=fold_test_yf,
    lambda=(preds+fold_test_bg)*fold_test_e))
  if (j %% 10 == 0) {
    print(paste0("Fold ", fid, ": Finished ratio ", j, "/", length(ratios)))
  }
}
save.image(paste0("real_data_cv_", fold_seed, "_", seed, ".RData"))

for (j in 1:nrow(grid_unit)) {
  ### Pull out of sample counts
  xf <- cbind(fold_test_xf, matrix(grid_unit[j,],
    nrow=nrow(fold_test_xf), ncol=ncol(grid_unit), byrow=TRUE))

  ### Predict lambda at out of sample using Scaled Vecchia
  preds <- predictions_scaled(svecfit, xf, m=25, joint=FALSE)

  ### Calculate CRPS (include time + background)
  crps_grid[j] <- mean(crps_pois(y=fold_test_yf,
    lambda=(preds+fold_test_bg)*fold_test_e))
  if (j %% 10 == 0) {
    print(paste0("Fold ", fid, ": Finished combo ", j, "/", nrow(grid)))
  }
}
save.image(paste0("real_data_cv_", fold_seed, "_", seed, ".RData"))
print(paste0("Fold ", fid, ": COMPLETE"))

res <- list(res=res, crps_pmfp=crps_pmfp, crps_ratio=crps_ratio,
  crps_grid=crps_grid, pmfp_grid=pmfps, ratio_grid=ratios, grid=grid,
  post_mean=post_mean, settings=list(seed=seed, ncvs=ncvs, fid=fid))

saveRDS(res, file=paste0("real_data_cv_fold", fid, "_seed", seed,
  format(Sys.time(), "_%Y%m%d%H%M%S"), ".rds"))
## remove saved RData files
unlink(paste0("real_data_cv_", fold_seed, "_", seed, ".RData"))
