library(laGP)
library(R.utils)
library(tidyverse)

source("../helper.R")

args <- R.utils::commandArgs(asValues=TRUE)
nmcs <- ifelse(!is.null(args$nmcs), as.numeric(args$nmcs), 100)

md <- readRDS(file="../data/Simulations.rds") %>%
  dplyr::rename(pmfp = parallel_mean_free_path) %>%
  dplyr::mutate(lratio = log(ratio)) %>%
  dplyr::mutate(sqratio = sqrt(ratio))
md[,c("x", "y", "z")] <- geo_to_spher_coords(md$lat, md$lon)
cpars <- md %>% dplyr::distinct(pmfp, sqratio)

for (i in 1:nmcs) {
  tic <- proc.time()[3]
  esa_lev <- sample(2:6, 1)
  mdesa <- md %>% dplyr::filter(ESA==esa_lev)
  n <- sample(10:50, 1)
  scales <- seq(1e-4, 1, by=1e-3)
  scale <- scales[sample(1:length(scales), 1)]
  Xscales <- c(1, 1, 1, rep(scale, 2))
  hout <- sample(1:nrow(cpars), 1)
  train <- mdesa %>% dplyr::filter(pmfp != cpars[hout,]$pmfp | sqratio != cpars[hout,]$sqratio)
  test <- mdesa %>% dplyr::filter(pmfp == cpars[hout,]$pmfp, sqratio == cpars[hout,]$sqratio)

  Xorig <- mdesa[,c("x", "y", "z", "pmfp", "sqratio")]
  X <- train %>% dplyr::select(x, y, z, pmfp, sqratio)
  Z <- train %>% dplyr::pull(blurred_ena_rate)
  XX <- test %>% dplyr::select(x, y, z, pmfp, sqratio)
  ZZ <- test %>% dplyr::pull(blurred_ena_rate)
  ## Code variables
  for (j in 1:ncol(X)) {
    X[,j] <- (X[,j] - min(Xorig[,j]))/(diff(range(Xorig[,j]))/Xscales[j])
    XX[,j] <- (XX[,j] - min(Xorig[,j]))/(diff(range(Xorig[,j]))/Xscales[j])
  }
  myseeds <- sample(1000000000:.Machine$integer.max, 2)
  set.seed(myseeds[1])
  d <- darg(d=NULL, X=X)
  set.seed(myseeds[2])
  g <- garg(g=NULL, Z)
  print(paste0("ESA: ", esa_lev, "; Holdout index: ", hout, "; PMFP: ",
               cpars[hout,]$pmfp, "; SqRatio: ", cpars[hout,]$sqratio,
               "; Neighborhood size: ", n, "; Scale: ", scale, "; dseed: ",
               myseeds[1], "; dstart: ", d$start, "; dmax: ", d$max, "; dmin: ",
               d$min, "; da: ", d$ab[1], "; db: ", d$ab[2], "; gseed: ",
               myseeds[2], "; gstart: ", g$start, "; gmax: ", g$max, "; gmin: ",
               g$min, "; ga: ", g$ab[1], "; db: ", g$ab[2]))
  fit <- aGPsep(X, Z, XX, end=n, omp.threads=1, d=d, g=g, verb=1, method="nn")
  print(paste("Finished iteration", i))
  print(paste("Time elapsed:", proc.time()[3]-tic))
}
