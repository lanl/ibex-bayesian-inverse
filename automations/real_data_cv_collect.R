library(scoringRules)

source("../helper.R")
source("../vecchia_scaled.R")

seed <- NA

args <- commandArgs(TRUE)
if (length(args) > 0) {
  for (i in 1:length(args)) {
    eval(parse(text=args[[i]]))
  }
}
if (is.na(seed)) {
  stop("Must provide a seed to collect results!")
}

files <- list.files(path="results",
  pattern=paste0("real_data_cv_fold[1-9][0-9]{0,1}_seed", seed, "_[0-9]*.rds"))
res <- readRDS(paste0("results/", files[1]))

us <- array(NA, dim=c(nrow(res$res$u), ncol(res$res$u), length(files)))
logscls <- matrix(NA, nrow=length(res$res$logscl), ncol=length(files))

crps_pmfp <- matrix(NA, ncol=length(files), nrow=length(res$crps_pmfp))
crps_ratio <- matrix(NA, ncol=length(files), nrow=length(res$crps_ratio))
crps_grid <- matrix(NA, ncol=length(files), nrow=length(res$crps_grid))

pmfp_grid <- res$pmfp_grid
ratio_grid <- res$ratio_grid
grid <- res$grid

post_means <- matrix(NA, ncol=2, nrow=length(files))
colnames(post_means) <- c("pmfp", "ratio")

settings <- list(seed=res$settings$seed, fids=rep(NA, length(files)))

for (i in 1:length(files)) {
  res <- readRDS(paste0("results/", files[i]))
  crps_pmfp[,i] <- res$crps_pmfp
  crps_ratio[,i] <- res$crps_ratio
  crps_grid[,i] <- res$crps_grid
  post_means[i,] <- res$post_mean
  settings$fids[i] <- res$settings$fid

  us[,,i] <- res$res$u
  logscls[,i] <- res$res$logscl
}

metrics <- list(crps_pmfp=crps_pmfp, crps_ratio=crps_ratio,
  crps_grid=crps_grid, pmfp_grid=pmfp_grid, ratio_grid=ratio_grid, grid=grid,
  post_means=post_means, us=us, logscls=logscls)
saveRDS(metrics, file=paste0("../results/real_data_cv_metrics",
  format(Sys.time(), "_%Y%m%d%H%M%S"), ".rds"))
