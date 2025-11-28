## list files for each simulated map
files <- list.files(pattern="pois_bayes_inv_res_pmfp[0-9]{3,4}_rat[0-9.]{4,5}_[0-9]{14}.rds")
if (length(files)==0) {
  files <- list.files(pattern="pois_bayes_inv_res_(2009|201[0-9]|202[0-1])+_[0-9]{14}.rds")
}
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

if (synth) {
  saveRDS(res, file=paste0("../papers/sim_calib_results",
    format(Sys.time(), "_%Y%m%d%H%M%S"), ".rds"))
} else if (length(files)==1) {
  saveRDS(res, file=paste0("../papers/real_calib_results_091011",
    format(Sys.time(), "_%Y%m%d%H%M%S"), ".rds"))
} else {
  saveRDS(res, file=paste0("../papers/real_calib_results_all_years",
    format(Sys.time(), "_%Y%m%d%H%M%S"), ".rds"))
}
file.remove(files)
