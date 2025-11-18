files <- list.files(path="results", pattern=paste0("scale_disc_test_scale_.*.rds"))

res <- list()
for (i in 1:length(files)) {
  iter_res <- readRDS(paste0("results/", files[i]))
  res[[i]] <- list(logscls=iter_res$mcmc_res$logscl, truth=iter_res$settings$scale)
}
saveRDS(metrics, file=paste0("results/scale_disc_test_results",
  format(Sys.time(), "_%Y%m%d%H%M%S"), ".rds"))
