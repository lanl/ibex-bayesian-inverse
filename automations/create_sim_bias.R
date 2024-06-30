scls <- c(seq(0.1, 1, by=0.1), seq(2, 20, by=2))
ibex_2020A <- read.csv("ibex_real_2020A.csv")
sims <- read.csv("sims.csv")
sims_sub <- sims[sims$parallel_mean_free_path==1750 & sims$ratio==0.02,]
sims_bias <- data.frame(matrix(NA, nrow=nrow(ibex_2020A)*length(scls),
 ncol=10))
colnames(sims_bias) <- c("parallel_mean_free_path", "ratio", "lon", "lat",
  "esa", "time", "blurred_ena_rate", "background", "sim_counts", "scl")
row <- 1
for (s in scls) {
  for (i in 1:nrow(ibex_2020A)) {
    sim_val <- sims_sub[ibex_2020A[i,]$near_sim,]
    sim_count <- rpois(1,(sim_val$blurred_ena_rate*s+ibex_2020A[i,]$background)*ibex_2020A[i,]$time)
    sims_bias[row,] <- c(1750,0.02,ibex_2020A[i,]$ecliptic_lon,
      ibex_2020A[i,]$ecliptic_lat,4,ibex_2020A[i,]$time,
      sim_val$blurred_ena_rate,ibex_2020A[i,]$background, sim_count,s)
    row <- row+1
  }
  print(s)
}

write.csv(sims_bias, file="sims_bias.csv", row.names=FALSE)
