library(laGP)
library(tidyverse)

source("../helper.R")

args <- R.utils::commandArgs(asValues=TRUE)
numruns <- ifelse(!is.null(args$runs), as.numeric(args$runs), 1)

md <- readRDS("../data/Simulations.rds") %>%
  dplyr::rename(pmfp=parallel_mean_free_path)
md[,c("x", "y", "z")] <- geo_to_spher_coords(md$lat, md$lon)
md <- md %>%
  dplyr::select(ESA, pmfp, ratio, x, y, z, blurred_ena_rate)
fd <- readRDS("../data/simulated_binned_direct_events_data.rds") %>%
  dplyr::rename(pmfp=parallel_mean_free_path)
fd[,c("x", "y", "z")] <- geo_to_spher_coords(fd$lat, fd$lon)
fd <- fd %>%
  dplyr::select(esa, pmfp, ratio, x, y, z)
cpars <- md %>% dplyr::distinct(pmfp, ratio)

times <- matrix(data=NA, nrow=numruns, ncol=8)

for (i in 1:numruns) {
  rcol <- 1
  esa_lev <- sample(2:6, 1)
  fpars <- cpars[sample(1:nrow(cpars), 1),]
  fpmfp <- fpars$pmfp
  fratio <- fpars$ratio

  X <- md %>% dplyr::filter(ESA==esa_lev) %>% dplyr::select(!c(ESA, blurred_ena_rate))
  Z <- md %>% dplyr::filter(ESA==esa_lev) %>% dplyr::pull(blurred_ena_rate)
  XX <- fd %>% dplyr::filter(esa==esa_lev, pmfp==fpmfp, ratio==fratio) %>%
    dplyr::select(!esa)
  Xall <- rbind(X, XX)

  for (j in 1:ncol(Xall)) {
    X[,j] <- (X[,j] - min(Xall[,j]))/diff(range(Xall[,j]))
    XX[,j] <- (XX[,j] - min(Xall[,j]))/diff(range(Xall[,j]))
  }

  ## Isotropic  vs. Separable Covariance Function
  tic <- proc.time()[3]
  isofit <- aGP(X=X, Z=Z, XX=XX, omp.threads=10, verb=0)
  toc <- proc.time()[3]
  times[i,rcol] <- toc-tic
  rcol <- rcol + 1
  print("Finished fitting isotropic covariance")
  
  tic <- proc.time()[3]
  sepfit <- aGPsep(X=X, Z=Z, XX=XX, omp.threads=10, verb=0)
  toc <- proc.time()[3]
  times[i,rcol] <- toc-tic
  rcol <- rcol + 1
  print("Finished fitting separable covariance")

  ## Neighborhood Size
  ns <- seq(10, 50, by=10)
  for (n in ns) {
    tic <- proc.time()[3]
    aGP(X=X, Z=Z, XX=XX, end=n, omp.threads=10, verb=0)
    toc <- proc.time()[3]
    times[i,rcol] <- toc-tic
    rcol <- rcol + 1
    print(paste0("Finished fitting neighborhood size ", n))
  }

  ## Splitting data into XX% chunks and making 100/XX aGPsep calls
  xsplits <- seq(0, 1, length=11)
  split_time <- 0
  for (j in 1:(length(xsplits)-1)) {
    subi <- which(X$x > (xsplits[j]-0.01) & X$x < (xsplits[j+1]+0.01))
    Xsub <- X[subi,]
    Zsub <- Z[subi]
    XXsub <- XX %>% dplyr::filter(x > xsplits[j], x < xsplits[j+1])
    tic <- proc.time()[3]
    subfit <- aGP(X=Xsub, Z=Zsub, XX=XXsub, omp.threads=10, verb=0)
    toc <- proc.time()[3]
    split_time <- split_time + toc-tic
    print(paste0("Finished fitting split number ", j))
  }
  times[i,rcol] <- split_time
  print(paste0("Finished iteration ", i, "/", numruns))
}
colnames(times) <- c("isotropic", "separable", "iso_n10", "iso_n20", "iso_n30",
                     "iso_n40", "iso_n50", "iso_split")

res <- list(times=times)
saveRDS(res, file=paste0("../results/laGP_speed_test_runs", numruns,
                         format(Sys.time(), "_%Y%m%d%H%M%S"), ".rds"))
