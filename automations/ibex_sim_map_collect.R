map <- "2020A"

args <- commandArgs(TRUE)
if (length(args) > 0) {
  for (i in 1:length(args)) {
    eval(parse(text=args[[i]]))
  }
}

ibex_real <- read.csv("../data/ibex_real.csv")
ibex_real <- ibex_real[ibex_real$map==map,]
n <- nrow(ibex_real)

res_files <- list.files(pattern=paste0("near_sim_", map, "_[0-9]*.rds"))

if (length(res_files) != n) {
  stop("Not all rows of the IBEX real data were executed")
}
res <- matrix(NA, nrow=length(res_files), ncol=2)

for (row in 1:length(res_files)) {
  near_sim <- readRDS(res_files[row])
  res[row,] <- c(near_sim[1], near_sim[2])
#  file.remove(res_files[row])
}
colnames(res) <- c("row", "near_sim")

write.csv(res, paste0('ibex_sim_map_', map, '.csv'), row.names=FALSE)
