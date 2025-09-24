
###############################################################################
###############################################################################
## Bivariate visual of model parameter posterior given real counts from the
## satellite. Bivariate posterior is generated for each year.
## Arranged in a 2x7 grid.
###############################################################################
###############################################################################

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
mtext(expression("Parallel Mean Free Path ("~u[1]~")"), side=1, outer=TRUE, line=4.25, cex=1.2)
mtext(expression("Ratio ("~u[2]~")"), side=2, outer=TRUE, line=3.0, cex=1.2)
dev.off()


###############################################################################
###############################################################################
## One bivariate visual of model parameter posterior given real counts from the
## satellite. Field data is satellite data from 2009-2011. This aligns with
## the GDF from the two-parameter simulation
###############################################################################
###############################################################################

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
  col=rev(heat.colors(128)), xlab=expression("Parallel Mean Free Path ("~u[1]~")"),
  ylab=expression("Ratio ("~u[2]~")"), xlim=c(500, 3000), ylim=c(0, 0.1))
abline(v=seq(500, 3000, by=500), col="lightgrey", lty=3)
abline(h=seq(0, 0.1, length=6), col="lightgrey", lty=3)
lines(cls$x, cls$y, lty=2)
dev.off()

###############################################################################
###############################################################################
## Real data visual containing six plots:
## - One plot of each year's satellite data from 2009-2011
## - One plot of predicted surrogate output at posterior mean of model
##   parameters given 2009-2011 data
## - Bivariate posterior of model parameters
###############################################################################
###############################################################################

library(MASS)
library(coda)
library(ks)

source("../helper.R")
source("../vecchia_scaled.R")

## Read in all model and field data
model_data <- read.csv(file="../data/sims.csv")
field_data <- read.csv(file="../data/ibex_real.csv")

## Rename field name columns to align with preprocessing step
colnames(field_data)[colnames(field_data)== "counts"] <- "sim_counts"
colnames(field_data)[colnames(field_data)== "ecliptic_lat"] <- "lat"
colnames(field_data)[colnames(field_data)== "ecliptic_lon"] <- "lon"
field_data$est_rate <- field_data$sim_counts/field_data$time - field_data$background
field_data <- field_data[which(!is.nan(field_data$est_rate)),]
field_data$nlon <- nose_center_lons(field_data$lon)
pd <- preprocess_data(md=model_data, fd=field_data, esa_lev=4,
  fparams=c("2020A"), scales=c(1, 1), tol=NA, quant=0.0,
  real=TRUE, disc=FALSE)
model_data$nlon <- nose_center_lons(model_data$lon)

## Load results of run on real data
res <- readRDS("final_results/sun_cycle_index_24.rds")
post_mean <- apply(res$mcmc_res$u[seq(1001, 10000, by=10),], 2, mean)

fit <- fit_scaled(y=pd$Zmod, inputs=as.matrix(cbind(pd$Xmod, pd$Umod)),
 nug=1e-4, ms=25)

XX_ll <- cbind(unique(model_data[,c("lon", "lat")]), matrix(post_mean, nrow=1))
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
cols <- colorRampPalette(c("blue", "cyan", "green", "yellow", "red", "magenta"))(128)
bks <- seq(predrange[1], predrange[2], length=length(cols)+1)
ylims <- range(model_data$lat)
xlims <- rev(range(model_data$nlon))

## Determine colors for plots based on breaks
field_rates <- cut(field_data$est_rate, breaks=bks, labels=FALSE)
field_rates[which(field_data$est_rate <= predrange[1])] <- 1
field_rates[which(field_data$est_rate >= predrange[2])] <- length(cols)
field_cols <- cols[field_rates]
field_data$col <- field_cols

field_data_09 <- field_data[field_data$map=="2009A",]
field_data_10 <- field_data[field_data$map=="2010A",]
field_data_11 <- field_data[field_data$map=="2011A",]

par(mfrow=c(1,1), mar=c(5.1, 4.1, 0.2, 0.2))
pdf("ibex_field_09.pdf", width=6.0, height=5)
plot(x=field_data_09$nlon, y=field_data_09$lat, col=field_data_09$col, pch=16,
  cex=0.7, xlab="Longitude", xaxt="n", ylab="Latitude", xlim=xlims,
  ylim=ylims, cex.lab=1.1)
axis(1, at=seq(325, 25, by=-60),
  labels=c(60, 0, 300, 240, 180, 120))
dev.off()

par(mfrow=c(1,1), mar=c(5.1, 4.1, 0.2, 0.2))
pdf("ibex_field_10.pdf", width=6.0, height=5)
plot(x=field_data_10$nlon, y=field_data_10$lat, col=field_data_10$col,
  pch=16, cex=0.7, xlab="Longitude", xaxt="n", ylab="Latitude", xlim=xlims,
  ylim=ylims, cex.lab=1.1)
axis(1, at=seq(325, 25, by=-60),
  labels=c(60, 0, 300, 240, 180, 120))
dev.off()

par(mfrow=c(1,1), mar=c(5.1, 4.1, 0.2, 0.2))
pdf("ibex_field_11.pdf", width=6.0, height=5)
plot(x=field_data_11$nlon, y=field_data_11$lat, col=field_data_11$col,
  pch=16, cex=0.7, xlab="Longitude", xaxt="n", ylab="Latitude", xlim=xlims,
  ylim=ylims, cex.lab=1.1)
axis(1, at=seq(325, 25, by=-60),
  labels=c(60, 0, 300, 240, 180, 120))
dev.off()

pred_lons <- sort(unique(pred_data$nlon))
pred_lats <- sort(unique(pred_data$lat))
pred_zmat <- xtabs(lhat_curr ~ nlon + lat, data=pred_data)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 0.2, 0.2))
pdf("ibex_surr_pred_real.pdf", width=6.0, height=5)
image(x=pred_lons, y=pred_lats, z=pred_zmat, col=cols, xlab="Longitude",
  xaxt="n", ylab="Latitude", breaks=bks, cex.lab=1.1, ylim=ylims, xlim=xlims)
axis(1, at=seq(325, 25, by=-60),
  labels=c(60, 0, 300, 240, 180, 120))
dev.off()

pmfps <- res$mcmc_res$u[seq(1001, 10000, by=10),1]*2500+500
ratios <- res$mcmc_res$u[seq(1001, 10000, by=10),2]*(0.1-0.001)+0.001

xy <- cbind(pmfps, ratios)
H <- Hpi(xy)*2
fhat <- kde(x=xy, H=H, xmin=c(500, 0), xmax=c(3000, 0.1),
  compute.cont=TRUE, gridsize=rep(1501, ncol(xy)))
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
par(mfrow=c(1,1), mar=c(5.1, 4.1, 0.2, 0.2))
pdf("ibex_real_post_est.pdf", width=5, height=5)
image(fhat$eval.points[[1]], fhat$eval.points[[2]], fhat$estimate,
  col=rev(heat.colors(128)), xlab=expression("Parallel Mean Free Path ("~u[1]~")"),
  ylab=expression("Ratio ("~u[2]~")"), xlim=c(500, 3000), ylim=c(0, 0.1))
abline(v=seq(500, 3000, by=500), col="lightgrey", lty=3)
abline(h=seq(0, 0.1, length=6), col="lightgrey", lty=3)
lines(cls$x, cls$y, lty=2)
dev.off()

# Plot contour at HPD threshold (zoomed in)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 0.2, 0.2))
pdf("ibex_real_post_est_zoom.pdf", width=5, height=5)
image(fhat$eval.points[[1]], fhat$eval.points[[2]], fhat$estimate,
  col=rev(heat.colors(128)), xlab=expression("Parallel Mean Free Path ("~u[1]~")"),
  ylab=expression("Ratio ("~u[2]~")"), xlim=c(2700, 3000), ylim=c(0, 0.004))
lines(cls$x, cls$y, lty=2)
dev.off()

###############################################################################
###############################################################################
## Plots showing discrepancy between real data and surrogate output at
## estimated model parameters
###############################################################################
###############################################################################

source("../helper.R")
source("../vecchia_scaled.R")

## Read in all model and field data
model_data <- read.csv(file="../data/sims.csv")
field_data <- read.csv(file="../data/ibex_real.csv")

## Rename field name columns to align with preprocessing step
colnames(field_data)[colnames(field_data)== "counts"] <- "sim_counts"
colnames(field_data)[colnames(field_data)== "ecliptic_lat"] <- "lat"
colnames(field_data)[colnames(field_data)== "ecliptic_lon"] <- "lon"
field_data$est_rate <- field_data$sim_counts/field_data$time - field_data$background
field_data <- field_data[which(!is.nan(field_data$est_rate)),]
field_data$nlon <- nose_center_lons(field_data$lon)
pd <- preprocess_data(md=model_data, fd=field_data, esa_lev=4,
  fparams=c("2020A"), scales=c(1, 1), tol=NA, quant=0.0,
  real=TRUE, disc=FALSE)
model_data$nlon <- nose_center_lons(model_data$lon)

## Load results of run on real data
res <- readRDS("final_results/sun_cycle_index_24.rds")
post_mean <- apply(res$mcmc_res$u[seq(1001, 10000, by=10),], 2, mean)

fit <- fit_scaled(y=pd$Zmod, inputs=as.matrix(cbind(pd$Xmod, pd$Umod)),
 nug=1e-4, ms=25)

XX_ll <- cbind(unique(model_data[,c("lon", "lat")]), matrix(post_mean, nrow=1))
XX_ll[,c("x", "y", "z")] <- geo_to_spher_coords(lat=XX_ll$lat, lon=XX_ll$lon)
XX_ll$x <- (XX_ll$x - min(XX_ll$x)) / diff(range(XX_ll$x))
XX_ll$y <- (XX_ll$y - min(XX_ll$y)) / diff(range(XX_ll$y))
XX_ll$z <- (XX_ll$z - min(XX_ll$z)) / diff(range(XX_ll$z))
XX <- XX_ll[,c("x", "y", "z", "1", "2")]
colnames(XX) <- c("x", "y", "z", "pmfp", "ratio")

lhat_curr <- predictions_scaled(fit, as.matrix(XX), m=25, joint=FALSE,
  predvar=FALSE)
pred_data <- data.frame(XX_ll, lhat_curr)
pred_data[,c("x", "y", "z")] <- geo_to_spher_coords(lat=pred_data$lat, lon=pred_data$lon)
pred_data$nlon <- nose_center_lons(pred_data$lon)

field_data <- field_data[field_data$map %in% c("2009A", "2010A", "2011A"),]
field_data$sim_rate <- NA
for (i in 1:nrow(field_data)) {
  fx <- field_data$x[i]
  fy <- field_data$y[i]
  fz <- field_data$z[i]
  min_dist <- 100000
  closest_ind <- NA

  for (j in 1:nrow(pred_data)) {
    dist <- sqrt((pred_data$x[j]-fx)^2 + (pred_data$y[j]-fy)^2 + (pred_data$z[j]-fz)^2)
    if (dist < min_dist) {
      min_dist <- dist
      closest_ind <- j
    }
  }
  field_data$sim_rate[i] <- pred_data$lhat_curr[closest_ind]

  if (i %% 100 == 0) {
    print(paste0("Finished point ", i, " out of ", nrow(field_data)))
  }
}
field_data$disc <- field_data$est_rate - field_data$sim_rate
max_disc <- max(abs(quantile(field_data$disc, probs=c(0.05, 0.95))))
disc_cols <- colorRampPalette(c("blue", "white", "red"))(128)
bks <- seq(-max_disc, max_disc, length=length(disc_cols)+1)
field_discs <- cut(field_data$disc, breaks=bks, labels=FALSE)
field_disc_cols <- disc_cols[field_discs]
field_data$disc_col <- field_disc_cols

field_data_09 <- field_data[field_data$map=="2009A",]
field_data_10 <- field_data[field_data$map=="2010A",]
field_data_11 <- field_data[field_data$map=="2011A",]

ylims <- range(model_data$lat)
xlims <- rev(range(model_data$nlon))

pdf("ibex_field_disc_09.pdf", width=7, height=5)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 0.2, 6))
plot(x=field_data_09$nlon, y=field_data_09$lat, col=field_data_09$disc_col,
  pch=16, cex=0.7, xlab="Longitude", xaxt="n", ylab="Latitude", xlim=xlims,
  ylim=ylims, cex.lab=1.1)
axis(1, at=seq(325, 25, by=-60),
  labels=c(60, 0, 300, 240, 180, 120))
fields::image.plot(zlim=c(-max_disc, max_disc), col=disc_cols,
  legend.only=TRUE, side=4, line=2, smallplot=c(0.86, 0.9, 0.3, 0.9))
dev.off()

pdf("ibex_field_disc_10.pdf", width=7, height=5)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 0.2, 6))
plot(x=field_data_10$nlon, y=field_data_10$lat, col=field_data_10$disc_col,
  pch=16, cex=0.7, xlab="Longitude", xaxt="n", ylab="Latitude", xlim=xlims,
  ylim=ylims, cex.lab=1.1)
axis(1, at=seq(325, 25, by=-60),
  labels=c(60, 0, 300, 240, 180, 120))
fields::image.plot(zlim=c(-max_disc, max_disc), col=disc_cols,
  legend.only=TRUE, side=4, line=2, smallplot=c(0.86, 0.9, 0.3, 0.9))
dev.off()

pdf("ibex_field_disc_11.pdf", width=7, height=5)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 0.2, 6))
plot(x=field_data_11$nlon, y=field_data_11$lat, col=field_data_11$disc_col,
  pch=16, cex=0.7, xlab="Longitude", xaxt="n", ylab="Latitude", xlim=xlims,
  ylim=ylims, cex.lab=1.1)
axis(1, at=seq(325, 25, by=-60),
  labels=c(60, 0, 300, 240, 180, 120))
fields::image.plot(zlim=c(-max_disc, max_disc), col=disc_cols,
  legend.only=TRUE, side=4, line=2, smallplot=c(0.86, 0.9, 0.3, 0.9))
dev.off()

pdf("ibex_mult_scale_disc.pdf", width=7, height=5)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 0.2, 6))
hist(exp(res$mcmc_res$logscl[seq(1001, 10000, 10)]),
  xlab="multiplicative scale", main="")
dev.off()

## Another option:
pdf("ibex_field_disc_09.pdf", width=7, height=5)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 0.2, 5))
plot(x=field_data_09$nlon, y=field_data_09$lat, col=field_data_09$disc_col,
  pch=16, cex=0.7, xlab="Longitude", xaxt="n", ylab="Latitude", xlim=xlims,
  ylim=ylims, cex.lab=1.1)
axis(1, at=seq(325, 25, by=-60),
  labels=c(60, 0, 300, 240, 180, 120))
fields::image.plot(zlim=c(-max_disc, max_disc), col=disc_cols,
  legend.only=TRUE, side=4, line=2, smallplot=c(0.86, 0.9, 0.3, 0.9))
dev.off()

pdf("ibex_mult_scale_disc.pdf", width=7, height=5)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 0.2, 5))
hist(exp(res$mcmc_res$logscl[seq(1001, 10000, 10)]),
  xlab="multiplicative scale", main="", breaks=12)
dev.off()
