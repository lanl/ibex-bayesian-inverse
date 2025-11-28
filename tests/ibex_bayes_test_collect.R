## list files for each simulated map
files <- list.files(pattern="pois_bayes_inv_res_pmfp[0-9]{3,4}_rat[0-9.]{4,5}_[0-9]{14}.rds")
res <- list()
res1 <- readRDS(files[1])
synth <- is.null(res1$year)
## read in results from each file
for (i in 1:length(files)) {
  iter_res <- readRDS(files[i])
  if (synth) {
    res[[i]] <- list(u=iter_res$u, truth=iter_res$truth)
  } else {
    res[[i]] <- list(u=iter_res$u, year=iter_res$year)
  }
}
saveRDS(res, file=paste0("../papers/sim_calib_results",
  format(Sys.time(), "_%Y%m%d%H%M%S"), ".rds"))
file.remove(files)
