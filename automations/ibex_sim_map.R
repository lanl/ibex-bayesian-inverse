library(R.utils)
library(doParallel)
library(foreach)
library(parallel)

source("../helper.R")

args <- R.utils::commandArgs(asValues=TRUE)

map <- ifelse(!is.null(args[["map"]]), as.character(args[["map"]]), "2020A")
row <- ifelse(!is.null(args[["row"]]), as.integer(args[["row"]]), 0)
vb <- ifelse(!is.null(args[["v"]]), TRUE, FALSE)
settings <- list(map=map, row=row)
if (vb) print(settings)

if (row==0) {
  stop("Must specify a row to find the nearest simulation value for")
}

sims <- read.csv("../data/sims.csv")
sims_sub <- sims[sims$parallel_mean_free_path==500 & sims$ratio==0.001,]
ibex_real <- read.csv("../data/ibex_real.csv")
ibex_real <- ibex_real[ibex_real$map==map,]
sims_sub[,c("x","y","z")] <- geo_to_spher_coords(lat=sims_sub$lat, lon=sims_sub$lon)
ibex_real[,c("x","y","z")] <- geo_to_spher_coords(lat=ibex_real$ecliptic_lat,
  lon=ibex_real$ecliptic_lon)

min_dist <- .Machine$integer.max
min_ind <- NA
ibex_x <- ibex_real[row,]$x
ibex_y <- ibex_real[row,]$y
ibex_z <- ibex_real[row,]$z
for (i in 1:nrow(sims_sub)) {
  sim_x <- sims_sub[i,]$x
  sim_y <- sims_sub[i,]$y
  sim_z <- sims_sub[i,]$z
  iter_dist <- sqrt((ibex_x-sim_x)^2 + (ibex_y-sim_y)^2 + (ibex_z-sim_z)^2)
  if (iter_dist < min_dist) {
    min_dist <- iter_dist
    min_ind <- i
  }
}

saveRDS(matrix(c(row, min_ind), nrow=1), file=paste0("near_sim_", map, "_", row, ".rds"))
