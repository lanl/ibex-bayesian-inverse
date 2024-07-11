library(R.utils)
library(doParallel)
library(foreach)
library(parallel)

source("../helper.R")

args <- R.utils::commandArgs(asValues=TRUE)

map <- ifelse(!is.null(args[["map"]]), as.character(args[["map"]]), "2020A")
ncores <- ifelse(!is.null(args[["ncores"]]), as.integer(args[["ncores"]]),
 ifelse(startsWith(Sys.info()[["nodename"]], "tc"), Sys.getenv("SLURM_NPROCS"),
 parallel::detectCores()-1))
vb <- ifelse(!is.null(args[["v"]]), TRUE, FALSE)
settings <- list(map=map, ncores=ncores)
if (vb) print(settings)

sims <- read.csv("../data/sims.csv")
sims_sub <- sims[sims$parallel_mean_free_path==500 & sims$ratio==0.001,]
ibex_real <- read.csv("../data/ibex_real.csv")
ibex_real <- ibex_real[ibex_real$map==map,]
sims_sub[,c("x","y","z")] <- geo_to_spher_coords(lat=sims_sub$lat, lon=sims_sub$lon)
ibex_real[,c("x","y","z")] <- geo_to_spher_coords(lat=ibex_real$ecliptic_lat,
  lon=ibex_real$ecliptic_lon)

cl <- makeCluster(ncores)
registerDoParallel(cl)

near_sims <- foreach(i = 1:nrow(ibex_real), .combine='c') %dopar% {
  min_dist <- .Machine$integer.max
  min_ind <- NA
  ibex_x <- ibex_real[i,]$x
  ibex_y <- ibex_real[i,]$y
  ibex_z <- ibex_real[i,]$z
  for (j in 1:nrow(sims_sub)) {
    sim_x <- sims_sub[j,]$x
    sim_y <- sims_sub[j,]$y
    sim_z <- sims_sub[j,]$z
    iter_dist <- sqrt((ibex_x-sim_x)^2 + (ibex_y-sim_y)^2 + (ibex_z-sim_z)^2)
    if (iter_dist < min_dist) {
      min_dist <- iter_dist
      min_ind <- j
    }
  }
  min_ind
}
ibex_real_mod <- cbind(ibex_real, near_sims)
colnames(ibex_real_mod) <- c(colnames(ibex_real), "near_sim")
write.csv(ibex_real_mod, file=paste0("../data/ibex_real_", map, ".csv"),
 row.names=FALSE)
