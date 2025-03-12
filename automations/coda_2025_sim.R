library(R.utils)

source("../helper.R")

args <- R.utils::commandArgs(asValues=TRUE)

map <- "2020A"
row1 <- ifelse(!is.null(args[["row1"]]), as.integer(args[["row1"]]), 0)
row2 <- ifelse(!is.null(args[["row2"]]), as.integer(args[["row2"]]), 0)
vb <- TRUE
settings <- list(map=map, row=row)
if (vb) print(settings)

sims <- read.csv("../data/sims.csv")
ibex_real <- read.csv("../data/ibex_real.csv")
sims_sub <- sims[sims$parallel_mean_free_path==500 & sims$ratio==0.001,]
ibex_real <- ibex_real[ibex_real$map=="2020A",]
sims_sub[,c("x","y","z")] <- geo_to_spher_coords(lat=sims_sub$lat, lon=sims_sub$lon)
ibex_real[,c("x","y","z")] <- geo_to_spher_coords(lat=ibex_real$ecliptic_lat,
  lon=ibex_real$ecliptic_lon)

res <- matrix(NA, nrow=row2-row1+1, ncol=2)
storage_row <- 1
for (row in row1:row2) {
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
  res[storage_row,] <- c(row, min_ind)
  if (row %% 10 == 0) print(paste0("Finished row ", row))
  storage_row <- storage_row + 1
}
saveRDS(res, file=paste0("near_sim_", map, "_", row1, "_", row2, ".rds"))
