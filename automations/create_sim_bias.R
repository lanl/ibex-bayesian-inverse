library(doParallel)
library(foreach)
library(parallel)

pmfps <- seq(500,3000, by=250)
ratios <- c(0.001, 0.005, 0.01, 0.02, 0.05, 0.1)
grid <- as.matrix(expand.grid(pmfps, ratios))
colnames(grid) <- c("pmfp", "ratio")

ncores <- parallel::detectCores()-1
cl <- parallel::makeCluster(ncores, outfile="../temp/log.txt")
doParallel::registerDoParallel(cl)

res <- foreach(i = 1:nrow(grid)) %dopar% {
  scls <- c(seq(0.1, 1, by=0.1), seq(2, 20, by=2))
  ibex_2020A <- read.csv("../data/ibex_real_2020A.csv")
  sims <- read.csv("../data/sims.csv")
  sims_sub <- sims[sims$parallel_mean_free_path==grid[i,1] & sims$ratio==grid[i,2],]
  sims_bias <- data.frame(matrix(NA, nrow=nrow(ibex_2020A)*length(scls),
   ncol=10))
  colnames(sims_bias) <- c("parallel_mean_free_path", "ratio", "lon", "lat",
    "esa", "time", "blurred_ena_rate", "background", "sim_counts", "scl")
  row <- 1
  for (s in scls) {
    for (j in 1:nrow(ibex_2020A)) {
      sim_val <- sims_sub[ibex_2020A[j,]$near_sim,]
      sim_count <- rpois(1,(sim_val$blurred_ena_rate*s+ibex_2020A[j,]$background)*ibex_2020A[j,]$time)
      sims_bias[row,] <- c(grid[i,1],grid[i,2],ibex_2020A[j,]$ecliptic_lon,
        ibex_2020A[j,]$ecliptic_lat,4,ibex_2020A[j,]$time,
        sim_val$blurred_ena_rate,ibex_2020A[j,]$background, sim_count,s)
      row <- row+1
    }
    print(s)
  }

  write.csv(sims_bias, file=paste0("../data/sims_bias_2020A_p", grid[i,1],
   "_r", grid[i,2], ".csv"), row.names=FALSE)
}
stopCluster(cl)
