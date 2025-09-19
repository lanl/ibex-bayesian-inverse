source("../helper.R")
source('../vecchia_scaled.R')

library(ggplot2)
library(MASS)
library(coda)
library(ks)

## Visuals for comparing simulated counts, "true" simulator output
## estimated simulator output via surrogate predictions
# res <- readRDS("final_results/sim_calib_results_20250902.rds")
res <- readRDS("final_results/sim_calib_results_20250915.rds")
single_index <- NA
single_pmfp <- 1750
single_ratio <- 0.02

for (i in 1:length(res)) {
  if (res[[i]]$truth[1]==single_pmfp && res[[i]]$truth[2]==single_ratio) {
    single_index <- i
    break
  }
}

pred_params <- apply(res[[single_index]]$u[seq(10001, 20000, by=10),], 2, mean)
pred_params[1] <- (pred_params[1] - 500)/2500
pred_params[2] <- (pred_params[2] - 0.001)/(0.1-0.001)

model_data <- read.csv(file="../data/sims.csv")
field_data <- read.csv(file="../data/sims_real.csv")
pd <- preprocess_data(md=model_data, fd=field_data, esa_lev=4,
  fparams=c(single_pmfp, single_ratio), scales=c(1, 1), tol=NA, quant=0.0,
  real=FALSE, disc=FALSE)

field_data <- field_data[field_data$parallel_mean_free_path==single_pmfp &
  field_data$ratio==single_ratio,]
field_data$est_rate <- field_data$sim_counts/field_data$time - field_data$background
field_data$nlon <- nose_center_lons(field_data$lon)
field_data <- field_data[which(!is.nan(field_data$est_rate)),]
model_data <- model_data[model_data$parallel_mean_free_path==single_pmfp &
  model_data$ratio==single_ratio,]
model_data$nlon <- nose_center_lons(model_data$lon)

fit <- fit_scaled(y=pd$Zmod, inputs=as.matrix(cbind(pd$Xmod, pd$Umod)),
 nug=1e-4, ms=25)
XX_ll <- cbind(unique(model_data[,c("lon", "lat")]), matrix(pred_params, nrow=1))
XX_ll[,c("x", "y", "z")] <- geo_to_spher_coords(lat=XX_ll$lat, lon=XX_ll$lon)
XX_ll$x <- (XX_ll$x - min(XX_ll$x)) / diff(range(XX_ll$x))
XX_ll$y <- (XX_ll$y - min(XX_ll$y)) / diff(range(XX_ll$y))
XX_ll$z <- (XX_ll$z - min(XX_ll$z)) / diff(range(XX_ll$z))
XX <- XX_ll[,c("x", "y", "z", "1", "2")]
colnames(XX) <- c("x", "y", "z", "pmfp", "ratio")
lhat_curr <- predictions_scaled(fit, as.matrix(XX), m=25, joint=FALSE,
  predvar=FALSE)
pred_data <- data.frame(XX_ll, lhat_curr)
pred_data$nlon <- nose_center_lons(pred_data$lon)

predrange <- range(model_data$blurred_ena_rate, na.rm=TRUE)
cols <- colorRampPalette(c("blue", "cyan", "green", "yellow", "red", "magenta"))(500)
bks <- seq(predrange[1], predrange[2], length=length(cols)+1)
ylims <- range(model_data$lat)
xlims <- rev(range(model_data$nlon))

model_lons <- sort(unique(model_data$nlon))
model_lats <- sort(unique(model_data$lat))
model_zmat <- xtabs(blurred_ena_rate ~ nlon + lat, data=model_data)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 0.2, 0.2))
pdf("ibex_sim_mod.pdf", width=7, height=5)
image(x=model_lons, y=model_lats, z=model_zmat, col=cols, xlab="Longitude",
  xaxt="n", ylab="Latitude", breaks=bks, cex.lab=1.1, ylim=ylims, xlim=xlims)
axis(1, at=seq(325, 25, by=-60),
  labels=c(60, 0, 300, 240, 180, 120))
dev.off()

field_lons <- sort(unique(field_data$nlon))
field_lats <- sort(unique(field_data$lat))
field_rates <- cut(field_data$est_rate, breaks=bks,
  labels=FALSE)
field_rates[which(field_data$est_rate <= predrange[1])] <- 1
field_rates[which(field_data$est_rate >= predrange[2])] <- length(cols)
field_cols <- cols[field_rates]
par(mfrow=c(1,1), mar=c(5.1, 4.1, 0.2, 0.2))
pdf("ibex_sim_field.pdf", width=7, height=5)
plot(x=field_data$nlon, y=field_data$lat, col=field_cols, pch=16, cex=0.7,
  xlab="Longitude", xaxt="n", ylab="Latitude", xlim=xlims, ylim=ylims,
  cex.lab=1.1)
axis(1, at=seq(325, 25, by=-60),
  labels=c(60, 0, 300, 240, 180, 120))
dev.off()

pred_lons <- sort(unique(pred_data$nlon))
pred_lats <- sort(unique(pred_data$lat))
pred_zmat <- xtabs(lhat_curr ~ nlon + lat, data=pred_data)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 0.2, 0.2))
pdf("ibex_sim_est.pdf", width=7, height=5)
image(x=pred_lons, y=pred_lats, z=pred_zmat, col=cols, breaks=bks,
  xlab="Longitude", xaxt="n", ylab="Latitude", useRaster=TRUE,
  xlim=xlims, ylim=ylims, cex.lab=1.1)
axis(1, at=seq(325, 25, by=-60),
  labels=c(60, 0, 300, 240, 180, 120))


dev.off()

iter_pmfp <- res[[single_index]]$u[seq(10001, 20000, by=10),1]
iter_ratio <- res[[single_index]]$u[seq(10001, 20000, by=10),2]

xy <- cbind(iter_pmfp, iter_ratio)
H <- Hpi(xy)
fhat <- kde(x=xy, H=H, xmin=c(500, 0.001), xmax=c(3000, 0.1),
  compute.cont=TRUE)
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
par(mfrow=c(1,1), mar=c(5.1, 4.1, 0.2, 0.2))
pdf("ibex_post_est.pdf", width=7, height=5)
image(fhat$eval.points[[1]], fhat$eval.points[[2]], fhat$estimate,
  col=rev(heat.colors(128)), xlab=expression("Parallel Mean Free Path ("~u[1]~")"),
  ylab=expression("Ratio ("~u[2]~")"), xlim=c(500, 3000), ylim=c(0, 0.1),
  cex.lab=1.1)
abline(v=seq(500, 3000, by=500), col="lightgrey", lty=3)
abline(h=seq(0, 0.1, length=6), col="lightgrey", lty=3)
lines(cls$x, cls$y, lty=2)
points(x=1750, 0.02, col=4, pch=8, cex=1.5, lwd=1.25)
dev.off()

#########################################################################
#########################################################################
#########################################################################

## Visuals for bivariate posteriors of model parameters
library(MASS)
library(coda)
library(ks)

ind_files <- FALSE
if (ind_files) {
  calib_files <- list.files(pattern="mcmc_res_esa4_nmcmc*")
  res <- list()
  for (i in 1:length(calib_files)) {
    iter_res <- readRDS(calib_files[[i]])
    res[[i]] <- list(
      u=data.frame(pmfp=iter_res$mcmc_res$u[,1]*2500+500,
        ratio=iter_res$mcmc_res$u[,2]*(0.1-0.001)+0.001),
      truth=iter_res$params)
  }
} else {
  res <- readRDS("sim_calib_results_20250915.rds")
}

pmfps <- seq(750, 2750, by=250)
ratios <- c(0.005, 0.01, 0.02, 0.05)
pmfp_rat_grid <- data.frame(matrix(NA, ncol=2, nrow=length(res)))
colnames(pmfp_rat_grid) <- c("pmfp", "ratio")
for (i in 1:length(res)) {
  pmfp_rat_grid[i,1] <- res[[i]]$truth[1]
  pmfp_rat_grid[i,2] <- res[[i]]$truth[2]
}
pmfp_labs <- seq(500, 2500, by=500)
ratio_labs <- seq(0.02, 0.1, length=5)

pdf("sim_bayes_inv_res.pdf", width=7, height=5)
par(mfrow=c(length(ratios), length(pmfps)),
  mar=c(0.25,0.25,0.25,0.15), oma=c(7,5,0.5,0.5))
for (i in 1:length(ratios)) {
  for (j in 1:length(pmfps)) {
    yticks <- j == 1
    xticks <- i == length(ratios)

    r <- ratios[i]
    p <- pmfps[j]

    index <- which(pmfp_rat_grid$pmfp==p & pmfp_rat_grid$ratio==r)

    iter_pmfp <- res[[index]]$u[seq(10001, 20000, by=10),1]
    iter_ratio <- res[[index]]$u[seq(10001, 20000, by=10),2]

    xy <- cbind(iter_pmfp, iter_ratio)
    H <- Hpi(xy)
    fhat <- kde(x=xy, H=H, xmin=c(500, 0.001), xmax=c(3000, 0.1),
      compute.cont=TRUE)
    dx <- diff(fhat$eval.points[[1]][1:2])
    dy <- diff(fhat$eval.points[[2]][1:2])

    # Flatten density values
    dens_vals <- sort(as.vector(fhat$estimate), decreasing=TRUE)
    cum_prob <- cumsum(dens_vals)*dx*dy

    # Threshold for 95% HPD
    thresh <- dens_vals[which(cum_prob >= 0.95)[1]]
    if (is.na(thresh)) {
      thresh <- dens_vals[length(cum_prob)]
    }

    cls <- contourLines(fhat$eval.points[[1]],
      fhat$eval.points[[2]], fhat$estimate, levels=thresh)

    # If multiple, keep the largest (by number of vertices)
    largest <- cls[[which.max(sapply(cls, function(cl) length(cl$x)))]]

    # Plot contour at HPD threshold
    image(fhat$eval.points[[1]], fhat$eval.points[[2]], fhat$estimate,
      col=rev(heat.colors(128)), xaxt="n", yaxt="n", xlab="", ylab="", main="",
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
    points(x=p, r, col=4, pch=8, cex=1.25, lwd=1.25)
  }
}
mtext("Parallel Mean Free Path", side=1, outer=TRUE, line=4.25, cex=1.2)
mtext("Ratio", side=2, outer=TRUE, line=3.0, cex=1.2)
dev.off()
