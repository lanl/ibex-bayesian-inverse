library(R.utils)
library(doParallel)
library(foreach)
library(parallel)

args <- R.utils::commandArgs(asValues=TRUE)

map <- ifelse(!is.null(args[["map"]]), as.character(args[["map"]]), "2020A")
ncores <- ifelse(!is.null(args[["ncores"]]), as.integer(args[["ncores"]]),
 parallel::detectCores()-1)
vb <- ifelse(!is.null(args[["v"]]), TRUE, FALSE)
settings <- list(map=map, ncores=ncores)
if (vb) print(settings)

sims <- read.csv("../data/sims.csv")
sims_sub <- sims[sims$parallel_mean_free_path==500 & sims$ratio==0.001,]
ibex_real <- read.csv("../data/ibex_real.csv")
ibex_real <- ibex_real[ibex_real$map==map,]

cl <- makeCluster(ncores)
registerDoParallel(cl)

near_sims <- foreach(i = 1:nrow(ibex_real), .combine='c') %dopar% {
  min_dist <- .Machine$integer.max
  min_ind <- NA
  ibex_lat <- ibex_real[i,]$ecliptic_lat
  ibex_lon <- ibex_real[i,]$ecliptic_lon
  for (j in 1:nrow(sims_sub)) {
    sim_lat <- sims_sub[j,]$lat
    sim_lon <- sims_sub[j,]$lon
    iter_dist <- acos(sin(ibex_lat)*sin(sim_lat)+cos(ibex_lat)*cos(sim_lat)*cos(sim_lon-ibex_lon))
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
