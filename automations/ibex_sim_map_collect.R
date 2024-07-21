res_files <- list.files(pattern="near_sim_[0-9]*.rds")

if (length(res_files) != 14154) {
  stop("Not all rows of the IBEX real data were executed")
}
res <- matrix(NA, nrow=14154, ncol=2)

for (row in 1:length(res_files)) {
  near_sim <- readRDS(res_files[row])
  res[row,] <- c(near_sim[1], near_sim[2])
#  file.remove(res_files[row])
}
colnames(res) <- c("row", "near_sim")

write.csv(res, paste0('ibex_sim_map.csv'), row.names=FALSE)
