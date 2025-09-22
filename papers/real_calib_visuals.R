## Visuals for bivariate posteriors of model parameters
library(MASS)
library(coda)
library(ks)

ind_files <- FALSE
if (ind_files) {
  calib_files <- list.files(
    pattern="mcmc_res_esa4_nmcmc25000_end25_scNA_quantNA_svecchia_real20*")
  res <- list()
  for (i in 1:length(calib_files)) {
    iter_res <- readRDS(calib_files[[i]])
    res[[i]] <- list(
      u=data.frame(pmfp=iter_res$mcmc_res$u[,1]*2500+500,
        ratio=iter_res$mcmc_res$u[,2]*(0.1-0.001)+0.001),
      year=as.integer(substr(iter_res$params, start=1, stop=4)))
  }
} else {
  res <- readRDS("real_calib_results_20250916.rds")
}

years <- 2009:2022
pmfp_labs <- seq(500, 2500, by=500)
ratio_labs <- seq(0.02, 0.1, length=5)

pdf("real_bayes_inv_res.pdf", width=7, height=4)
par(mfrow=c(2, 7), mar=c(0.25,0.25,1.25,0.15),
  oma=c(7,5,0.5,0.5))
for (i in 1:length(years)) {
  yticks <- i %% (length(years)/2)-1 == 0
  xticks <- i > (length(years)/2)

  year <- years[i]
  index <- which(years==year)
  iter_pmfp <- res[[index]]$u[seq(5001, 25000, by=10),1]
  iter_ratio <- res[[index]]$u[seq(5001, 25000, by=10),2]

  xy <- cbind(iter_pmfp, iter_ratio)
  H <- Hpi(xy)/5
  fhat <- kde(x=xy, H=H, xmin=c(500, 0.001), xmax=c(3000, 0.1),
    compute.cont=TRUE, gridsize=rep(401, ncol(xy)))
  fhat$estimate <- pmax(fhat$estimate, 0)
  dx <- diff(fhat$eval.points[[1]][1:2])
  dy <- diff(fhat$eval.points[[2]][1:2])

  # Flatten density values
  dens_vals <- sort(as.vector(fhat$estimate), decreasing=TRUE)
  cum_prob <- cumsum(dens_vals)*dx*dy

  # Threshold for 95% HPD
  thresh <- dens_vals[which(cum_prob >= 0.95)[1]]
  if (is.na(thresh)) {
    thresh <- dens_vals[length(cum_prob)]
    if (thresh==0) {
      thresh <- dens_vals[max(which(dens_vals > 0))]
    }
  }

  cls <- contourLines(fhat$eval.points[[1]],
    fhat$eval.points[[2]], fhat$estimate, levels=thresh)

  # If multiple, keep the largest (by number of vertices)
  largest <- cls[[which.max(sapply(cls, function(cl) length(cl$x)))]]

    # Plot contour at HPD threshold
  image(fhat$eval.points[[1]], fhat$eval.points[[2]], fhat$estimate,
    col=rev(heat.colors(128)), xaxt="n", yaxt="n", xlab="", ylab="", main=year,
    xlim=c(500, 3000), ylim=c(0, 0.1))
  abline(v=seq(500, 3000, by=500), col="lightgrey", lty=3)
  abline(h=seq(0, 0.1, length=6), col="lightgrey", lty=3)
  if (xticks) {
    axis(1, labels=FALSE, tck=-0.05)
    text(pmfp_labs, par("usr")[3]-0.0075, labels=pmfp_labs, srt=90, adj=1, xpd=NA,
      cex=1.05)
  }
  if (yticks) {
    axis(2, labels=FALSE, tck=-0.05)
    text(x=par("usr")[1]-225, y=ratio_labs, labels=ratio_labs, adj=1,
      cex=1.05, xpd=NA)
  }
  lines(largest$x, largest$y, lty=2)
}
mtext("Parallel Mean Free Path", side=1, outer=TRUE, line=4.25, cex=1.2)
mtext("Ratio", side=2, outer=TRUE, line=3.0, cex=1.2)
dev.off()


## Visuals for bivariate posteriors of model parameters (years 2009-2011)
library(MASS)
library(coda)
library(ks)

res <- readRDS("final_results/mcmc_res_esa4_nmcmc30000_real2009A2010A2011A_20250917145959.rds")

pdf("suncyc_real_bayes_inv_res.pdf", width=7, height=4)
par(mfrow=c(2, 7), mar=c(0.25,0.25,1.25,0.15),
  oma=c(7,5,0.5,0.5))

pmfps <- res$mcmc_res$u[seq(20001, 30000, by=10),1]*2500+500
ratios <- res$mcmc_res$u[seq(20001, 30000, by=10),2]*(0.1-0.001)+0.001

xy <- cbind(pmfps, ratios)
H <- Hpi(xy)*2
fhat <- kde(x=xy, H=H, xmin=c(500, 0.001), xmax=c(3000, 0.1),
  compute.cont=TRUE, gridsize=rep(401, ncol(xy)))
fhat$estimate <- pmax(fhat$estimate, 0)
dx <- diff(fhat$eval.points[[1]][1:2])
dy <- diff(fhat$eval.points[[2]][1:2])

# Flatten density values
dens_vals <- sort(as.vector(fhat$estimate), decreasing=TRUE)
cum_prob <- cumsum(dens_vals)*dx*dy

# Threshold for 95% HPD
thresh <- dens_vals[which(cum_prob >= 0.95)[1]]
cls <- contourLines(fhat$eval.points[[1]],
  fhat$eval.points[[2]], fhat$estimate, levels=thresh)[[1]]

# Plot contour at HPD threshold
image(fhat$eval.points[[1]], fhat$eval.points[[2]], fhat$estimate,
  col=rev(heat.colors(128)), xlab="Parallel Mean Free Path", ylab="Ratio",
  xlim=c(500, 3000), ylim=c(0, 0.1))
abline(v=seq(500, 3000, by=500), col="lightgrey", lty=3)
abline(h=seq(0, 0.1, length=6), col="lightgrey", lty=3)
lines(cls$x, cls$y, lty=2)
dev.off()
