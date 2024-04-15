library(laGP)
library(R.utils)
library(tidyverse)

source("../helper.R")

args <- R.utils::commandArgs(asValues=TRUE)
esa_lev <- ifelse(!is.null(args$esa), as.numeric(args$esa), 3)

md <- readRDS(file="../data/Simulations.rds") %>%
  dplyr::rename(pmfp = parallel_mean_free_path) %>%
  dplyr::mutate(lratio = log(ratio)) %>%
  dplyr::mutate(sqratio = sqrt(ratio))
md[,c("x", "y", "z")] <- geo_to_spher_coords(md$lat, md$lon)
mdesa <- md %>% dplyr::filter(ESA==esa_lev)
cpars <- mdesa %>%
  dplyr::distinct(pmfp, sqratio)

hout <- sample(1:nrow(cpars), 1)
train <- mdesa %>% dplyr::filter(pmfp != cpars[hout,]$pmfp | sqratio != cpars[hout,]$sqratio)
test <- mdesa %>% dplyr::filter(pmfp == cpars[hout,]$pmfp, sqratio == cpars[hout,]$sqratio)

Xorig <- mdesa[,c("x", "y", "z", "pmfp", "sqratio")]
ns <- seq(10, 50, by=5)
scales <- 1e-4*2^(1:13)
Xscales <- matrix(NA, nrow=length(scales), ncol=ncol(Xorig))
Xscales[,1:3] <- 1
Xscales[,4] <- Xscales[,5] <- scales

metrics <- matrix(data=NA, nrow=length(ns)*nrow(Xscales), ncol=22)
colnames(metrics) <- c("nsize", "scale", "avgx", "minx", "maxx", "xlssd",
                       "avgy", "miny", "maxy", "ylssd", "avgz", "minz", "maxz",
                       "zlssd", "avgpmfps", "minpmfp", "maxpmfp", "pmfplssd",
                       "avgratio", "minratio", "maxratio", "ratiolssd")
tic <- proc.time()[3]
row <- 1
for (j in 1:nrow(Xscales)) {
  for (n in ns) {
    X <- train %>% dplyr::select(x, y, z, pmfp, sqratio)
    Z <- train %>% dplyr::pull(blurred_ena_rate)
    XX <- test %>% dplyr::select(x, y, z, pmfp, sqratio)
    ZZ <- test %>% dplyr::pull(blurred_ena_rate)

    ## Code variables
    for (i in 1:ncol(X)) {
      X[,i] <- (X[,i] - min(Xorig[,i]))/(diff(range(Xorig[,i]))/Xscales[j,i])
      XX[,i] <- (XX[,i] - min(Xorig[,i]))/(diff(range(Xorig[,i]))/Xscales[j,i])
    }
    d <- darg(d=NULL, X=X)
    print(paste("Lengthscale prior information:", d))
    fit <- aGPsep(X, Z, XX, end=n, omp.threads=14, d=d, verb=1,
                  method="nn")
    xs <- ys <- zs <- pmfps <- ratios <- rep(NA, nrow(XX))
    for (i in 1:nrow(XX)) {
      xs[i] <- length(unique(X[fit$Xi[i,],]$x))
      ys[i] <- length(unique(X[fit$Xi[i,],]$y))
      zs[i] <- length(unique(X[fit$Xi[i,],]$z))
      pmfps[i] <- length(unique(X[fit$Xi[i,],]$pmfp))
      ratios[i] <- length(unique(X[fit$Xi[i,],]$sqratio))
    }
    ls <- apply(fit$mle, 2, sd)
    metrics[row,] <- c(n, Xscales[j,4], mean(xs), range(xs), ls[1], mean(ys),
                       range(ys), ls[2], mean(zs), range(zs), ls[3], mean(pmfps),
                       range(pmfps), ls[4], mean(ratios), range(ratios), ls[5])
    print(paste0("Finished combination ", row, "/", nrow(metrics)))
    print(paste0("Neighborhood size: ", n, "; Scale: ", Xscales[j,4]))
    row <- row + 1
    saveRDS(metrics, file="../results/temp_scale_inputs_test.rds")
  }
}
toc <- proc.time()[3]
res <- list(metrics=metrics, time=toc-tic)
saveRDS(res, file=paste0("../results/input_scales_test_",
                         format(Sys.time(), "%Y%m%d%H%M%S"), ".rds"))
