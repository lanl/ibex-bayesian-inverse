files <- list.files(pattern="scale_disc_test_scale_.*.rds")

res <- list()
for (i in 1:length(files)) {
  iter_res <- readRDS(files[i])
  res[[i]] <- list(logscls=iter_res$res$logscl, truth=iter_res$settings$scale)
}
saveRDS(res, file=paste0("scale_disc_test_results",
  format(Sys.time(), "_%Y%m%d%H%M%S"), ".rds"))
