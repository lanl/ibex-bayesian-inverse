###############################################################################
###############################################################################
## Figures for running our Poisson Bayesian inverse problem framework on
## real data collected from the IBEX satellite
###############################################################################
###############################################################################

###############################################################################
## FIGURE 13: Bivariate visual of model parameter posterior given real counts
## from the satellite. Bivariate posterior is generated for each year, arranged
## in a 2x7 grid
## DATA NEEDED: real_calib_results_all_years.rds
###############################################################################

library(MASS)
library(coda)
library(ks)

res <- readRDS("final_results/real_calib_results_all_years.rds")

years <- 2012:2021
pmfp_labs <- seq(500, 2500, by=500)
ratio_labs <- seq(0.02, 0.1, length=5)

## Figure 13
pdf("real_bayes_inv_res.pdf", width=7, height=4)
par(mfrow=c(2, 5), mar=c(0.25,0.25,1.25,0.15),
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
    compute.cont=TRUE, gridsize=rep(201, ncol(xy)))
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
mtext(expression("Parallel Mean Free Path ("~u[1]~")"), side=1, outer=TRUE,
  line=4.25, cex=1.2)
mtext(expression("Ratio ("~u[2]~")"), side=2, outer=TRUE, line=3.0, cex=1.2)
dev.off()

###############################################################################
## FIGURE 11: Real data visual containing six plots:
## - One plot of each year's satellite data from 2009-2011
## - One plot of predicted surrogate output at posterior mean of model
##   parameters given 2009-2011 data
## - Bivariate posterior of model parameters
## DATA NEEDED: sims.csv, ibex_real.csv, real_calib_results_091011.rds
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
res <- readRDS("final_results/real_calib_results_091011.rds")
post_mean <- apply(res$mcmc_res$u[seq(1001, 10000, by=10),], 2, mean)

## Fit surrogate on simulator output
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

predrange <- quantile(model_data$blurred_ena_rate, probs=c(0.00015, 0.9985))
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

## Figure 11 (top left panel)
pdf("ibex_field_09.pdf", width=6.0, height=5)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 7.1), mgp=c(2.4, 0.6, 0))
plot(x=field_data_09$nlon, y=field_data_09$lat, col=field_data_09$col, pch=16,
  cex=0.7, xlab="Longitude", xaxt="n", ylab="Latitude", xlim=xlims,
  ylim=ylims, cex.lab=1.1)
axis(1, at=seq(325, 25, by=-60),
  labels=c(60, 0, 300, 240, 180, 120))
fields::image.plot(zlim=predrange, col=cols, legend.lab="ENAs/sec", legend.line=3,
  legend.only=TRUE, side=4, line=2, smallplot=c(0.8, 0.84, 0.3, 0.75))
dev.off()

## Figure 11 (top middle panel)
pdf("ibex_field_10.pdf", width=6.0, height=5)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 7.1), mgp=c(2.4, 0.6, 0))
plot(x=field_data_10$nlon, y=field_data_10$lat, col=field_data_10$col,
  pch=16, cex=0.7, xlab="Longitude", xaxt="n", ylab="Latitude", xlim=xlims,
  ylim=ylims, cex.lab=1.1)
axis(1, at=seq(325, 25, by=-60),
  labels=c(60, 0, 300, 240, 180, 120))
fields::image.plot(zlim=predrange, col=cols, legend.lab="ENAs/sec", legend.line=3,
  legend.only=TRUE, side=4, line=2, smallplot=c(0.8, 0.84, 0.3, 0.75))
dev.off()

## Figure 11 (bottom left panel)
pdf("ibex_field_11.pdf", width=6.0, height=5)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 7.1), mgp=c(2.4, 0.6, 0))
plot(x=field_data_11$nlon, y=field_data_11$lat, col=field_data_11$col,
  pch=16, cex=0.7, xlab="Longitude", xaxt="n", ylab="Latitude", xlim=xlims,
  ylim=ylims, cex.lab=1.1)
axis(1, at=seq(325, 25, by=-60),
  labels=c(60, 0, 300, 240, 180, 120))
fields::image.plot(zlim=predrange, col=cols, legend.lab="ENAs/sec", legend.line=3,
  legend.only=TRUE, side=4, line=2, smallplot=c(0.8, 0.84, 0.3, 0.75))
dev.off()

pred_lons <- sort(unique(pred_data$nlon))
pred_lats <- sort(unique(pred_data$lat))
pred_zmat <- xtabs(lhat_curr ~ nlon + lat, data=pred_data)

## Figure 11 (bottom middle panel)
pdf("ibex_surr_pred_real.pdf", width=6.0, height=5)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 7.1), mgp=c(2.4, 0.6, 0))
## If NOT using pdf(), image will be flipped because of useRaster=TRUE
image(x=pred_lons, y=pred_lats, z=pred_zmat, col=cols, xlab="Longitude",
  xaxt="n", ylab="Latitude", breaks=bks, cex.lab=1.1, ylim=ylims,
  xlim=xlims, useRaster=TRUE)
axis(1, at=seq(325, 25, by=-60),
  labels=c(60, 0, 300, 240, 180, 120))
dev.off()

## create kernel density of posterior samples
pmfps <- res$mcmc_res$u[seq(1001, 10000, by=10),1]*2500+500
ratios <- res$mcmc_res$u[seq(1001, 10000, by=10),2]*(0.1-0.001)+0.001
xy <- cbind(pmfps, ratios)
H <- Hpi(xy)*2
fhat <- kde(x=xy, H=H, xmin=c(500, 0), xmax=c(3000, 0.1),
  compute.cont=TRUE, gridsize=rep(301, ncol(xy)))
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

## Figure 11 (top right panel)
# Plot contour at HPD threshold
pdf("ibex_real_post_est.pdf", width=5, height=5)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(2.4, 0.6, 0))
image(fhat$eval.points[[1]], fhat$eval.points[[2]], fhat$estimate,
  col=rev(heat.colors(128)), xlab=expression("Parallel Mean Free Path ("~u[1]~")"),
  ylab=expression("Ratio ("~u[2]~")"), xlim=c(500, 3000), ylim=c(0, 0.1))
abline(v=seq(500, 3000, by=500), col="lightgrey", lty=3)
abline(h=seq(0, 0.1, length=6), col="lightgrey", lty=3)
lines(cls$x, cls$y, lty=2)
dev.off()

## create kernel density of posterior samples, but zoomed in
fhat <- kde(x=xy, H=H, xmin=c(500, 0), xmax=c(3000, 0.1),
  compute.cont=TRUE, gridsize=rep(1001, ncol(xy)))
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

## Figure 11 (bottom right panel)
# Plot contour at HPD threshold (zoomed in)
pdf("ibex_real_post_est_zoom.pdf", width=5, height=5)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(2.4, 0.6, 0))
## If NOT using pdf(), image will be flipped because of useRaster=TRUE
image(fhat$eval.points[[1]], fhat$eval.points[[2]], fhat$estimate,
  col=rev(heat.colors(128)), useRaster=TRUE,
  xlab=expression("Parallel Mean Free Path ("~u[1]~")"),
  ylab=expression("Ratio ("~u[2]~")"), xlim=c(2700, 3000), ylim=c(0, 0.004))
lines(cls$x, cls$y, lty=2)
dev.off()

###############################################################################
## FIGURE 12: Plots of CRPS for 10-fold cross validation on IBEX real data
## DATA NEEDED: real_calib_results_091011.rds, real_data_cv_metrics.rds
###############################################################################

real_dat_res <- readRDS("final_results/real_calib_results_091011.rds")
post_mean <- apply(real_dat_res$mcmc_res$u[seq(1001, 10000, by=10),], 2, mean)
post_mean[1] <- post_mean[1]*2500+500
post_mean[2] <- post_mean[2]*(0.1-0.001)+0.001

cv_res <- readRDS("final_results/real_data_cv_metrics.rds")

crps_range <- range(c(apply(cv_res$crps_pmfp, 1, mean),
  apply(cv_res$crps_ratio, 1, mean)))

## Figure 12 (left panel)
pdf("crps_pmfp.pdf", width=5, height=5)
par(mgp=c(2.25, 0.8, 0))
plot(x=cv_res$pmfp_grid, y=apply(cv_res$crps_pmfp, 1, mean), type="l",
  xlab=expression("Parallel Mean Free Path ("~u[1]~")"), ylab="crps",
  ylim=crps_range)
abline(v=post_mean[1], col=2, lty=2, lwd=2.5)
legend("topleft", c("mean across folds", expression(hat(u[i])~"(2009-2011)")),
  col=1:2, lty=1:2, lwd=c(1, 2.5), bg="white", cex=1.15)
dev.off()

## Figure 12 (middle panel)
pdf("crps_ratio.pdf", width=5, height=5)
par(mgp=c(2.25, 0.8, 0))
plot(x=cv_res$ratio_grid, y=apply(cv_res$crps_ratio, 1, mean), type="l",
  xlab=expression("Ratio ("~u[2]~")"), ylab="crps", ylim=crps_range)
abline(v=post_mean[2], col=2, lty=2, lwd=2.5)
dev.off()

## Figure 12 (right panel)
pdf("crps_grid.pdf", width=5, height=5)
par(mgp=c(2.25, 0.8, 0))
unique_pmfps <- sort(unique(cv_res$grid[,1]))
unique_ratios <- sort(unique(cv_res$grid[,2]))
## If NOT using pdf(), image will be flipped because of useRaster=TRUE
image(x=unique_pmfps, y=unique_ratios, useRaster=TRUE,
  z=matrix(apply(cv_res$crps_grid, 1, mean), ncol=length(unique_pmfps)),
  col=heat.colors(128), xlab=expression("Parallel Mean Free Path ("~u[1]~")"),
  ylab=expression("Ratio ("~u[2]~")"))
dev.off()

###############################################################################
## FIGURE 14: Plots showing discrepancy between real data and surrogate output
## at estimated model parameters
## DATA NEEDED: sims.csv, ibex_real.csv, real_calib_results_091011.rds
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
res <- readRDS("final_results/real_calib_results_091011.rds")
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
XXF <- field_data[,c("x","y","z")]
XXF$x <- (field_data$x - min(field_data$x))/diff(range(field_data$x))
XXF$y <- (field_data$y - min(field_data$y))/diff(range(field_data$y))
XXF$z <- (field_data$z - min(field_data$z))/diff(range(field_data$z))
XXF <- cbind(XXF, matrix(post_mean, nrow=1))
colnames(XXF) <- c("x", "y", "z", "pmfp", "ratio")
XXF <- as.matrix(XXF)
field_data$sim_rate <- predictions_scaled(fit, XXF, m=25, joint=FALSE, predvar=FALSE)

field_data$disc <- field_data$est_rate - field_data$sim_rate
max_disc <- max(abs(quantile(field_data$disc, probs=c(0.05, 0.95))))
disc_cols <- colorRampPalette(c("blue", "white", "red"))(128)
bks <- seq(-max_disc, max_disc, length=length(disc_cols)+1)
field_discs <- cut(field_data$disc, breaks=bks, labels=FALSE)
field_discs[which(field_data$disc <= -max_disc)] <- 1
field_discs[which(field_data$disc >= max_disc)] <- length(disc_cols)
field_disc_cols <- disc_cols[field_discs]
field_data$disc_col <- field_disc_cols

field_data_09 <- field_data[field_data$map=="2009A",]
field_data_10 <- field_data[field_data$map=="2010A",]
field_data_11 <- field_data[field_data$map=="2011A",]

ylims <- range(model_data$lat)
xlims <- rev(range(model_data$nlon))

## Figure 14 (top left panel)
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

## Figure 14 (top right panel)
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

## Figure 14 (bottom left panel)
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

## Figure 14 (bottom right panel)
pdf("ibex_mult_scale_disc.pdf", width=7, height=5)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 0.2, 6))
hist(exp(res$mcmc_res$logscl[seq(1001, 10000, 10)]),
  xlab="multiplicative scale", main="")
dev.off()

###############################################################################
## FIGURE 16: PIT histograms for estimated sky maps generated from estimates
## of u for actual IBEX satellite data
## DATA NEEDED: real_calib_results_091011.rds, sims.csv, ibex_real.csv,
## real_calib_results_all_years.rds
###############################################################################

source("../helper.R")
source('../vecchia_scaled.R')

library(MASS)
library(coda)
library(ks)

# Create plot for years 2009-2011
### Read in calibration results
res <- readRDS("final_results/real_calib_results_091011.rds")

### Fit Scaled Vecchia GP surrogate
model_data <- read.csv(file="../data/sims.csv")
field_data <- read.csv(file="../data/ibex_real.csv")
field_data <- field_data[field_data$map %in% paste0(2009:2011, "A"),]
field_data$lat <- field_data$ecliptic_lat
field_data$lon <- field_data$ecliptic_lon
field_data$sim_counts <- field_data$counts
Xmod <- model_data[,c("lat", "lon")]
Xmod[,c("x", "y", "z")] <- geo_to_spher_coords(Xmod$lat, Xmod$lon)
Xmod <- Xmod[,c("x", "y", "z")]
Xfield <- unique(field_data[c("lon", "lat")])
Xfield[,c("x", "y", "z")] <- geo_to_spher_coords(Xfield$lat, Xfield$lon)
Xfield <- Xfield[,c("x", "y", "z")]
Xall <- rbind(Xmod, Xfield)
for (i in 1:ncol(Xall)) {
  Xmod[,i] <- (Xmod[,i] - min(Xall[,i]))/diff(range(Xall[,i]))
  Xfield[,i] <- (Xfield[,i] - min(Xall[,i]))/diff(range(Xall[,i]))
}

pd <- preprocess_data(md=model_data, fd=field_data, esa_lev=4,
  fparams=paste0(2009:2011, "A"), scales=c(1, 1), tol=NA, quant=0.0,
  real=TRUE, disc=FALSE)
fit <- fit_scaled(y=pd$Zmod, inputs=as.matrix(cbind(pd$Xmod, pd$Umod)),
 nug=1e-4, ms=25)

##### Predict using the GP surrogate at the estimated calibration parameters and XF
u_est <- apply(res$mcmc_res$u[seq(1001, 10000, by=10),], 2, mean)
XX <- as.matrix(cbind(Xfield, matrix(u_est, nrow=1)))
colnames(XX) <- NULL
preds <- predictions_scaled(fit, as.matrix(XX), m=25, joint=FALSE,
  predvar=FALSE)
##### Calculate PIT values
Fy  <- ppois(field_data$sim_counts, (preds+field_data$background)*field_data$time)
Fy1 <- ppois(field_data$sim_counts - 1, (preds+field_data$background)*field_data$time)
# Randomized PIT for discrete distributions
pit <- Fy1 + runif(length(field_data$sim_counts)) * (Fy - Fy1)

# Create data for plots for all years
### Read in calibration results
res <- readRDS("final_results/real_calib_results_all_years.rds")
### Read in data
model_data <- read.csv(file="../data/sims.csv")
Xmod <- model_data[,c("lat", "lon")]
Xmod[,c("x", "y", "z")] <- geo_to_spher_coords(Xmod$lat, Xmod$lon)
Xmod <- Xmod[,c("x", "y", "z")]
field_data <- read.csv(file="../data/ibex_real.csv")
field_data$lat <- field_data$ecliptic_lat
field_data$lon <- field_data$ecliptic_lon
field_data$sim_counts <- field_data$counts

pits <- list()

ymax <- 1.05
for (i in 4:(length(res)-1)) {
  iter_res <- res[[i]]
  ### Fit Scaled Vecchia GP surrogate
  Xmod_iter <- Xmod
  iter_field <- field_data[field_data$map==paste0(iter_res$year, "A"),]
  Xfield <- iter_field[c("lon", "lat")]
  Xfield[,c("x", "y", "z")] <- geo_to_spher_coords(Xfield$lat, Xfield$lon)
  Xfield <- Xfield[,c("x", "y", "z")]
  Xall <- rbind(Xmod_iter, Xfield)
  for (j in 1:ncol(Xall)) {
    Xmod_iter[,j] <- (Xmod_iter[,j] - min(Xall[,j]))/diff(range(Xall[,j]))
    Xfield[,j] <- (Xfield[,j] - min(Xall[,j]))/diff(range(Xall[,j]))
  }
  pd <- preprocess_data(md=model_data, fd=iter_field, esa_lev=4,
    fparams=c(paste0(iter_res$year, "A")), scales=c(1, 1), tol=NA, quant=0.0,
    real=TRUE, disc=FALSE)
  fit <- fit_scaled(y=pd$Zmod, inputs=as.matrix(cbind(pd$Xmod, pd$Umod)),
   nug=1e-4, ms=25)

  u_est <- apply(iter_res$u[seq(15001, 25000, by=10),], 2, mean)
  u_unit <- c((u_est[1] - 500)/2500, (u_est[2] - 0.001)/(0.1-0.001))
  XX <- as.matrix(cbind(Xfield, matrix(u_unit, nrow=1)))
  colnames(XX) <- NULL
  preds <- predictions_scaled(fit, as.matrix(XX), m=25, joint=FALSE,
    predvar=FALSE)

  ##### Calculate PIT values
  Fy  <- ppois(iter_field$sim_counts, (preds+iter_field$background)*iter_field$time)
  Fy1 <- ppois(iter_field$sim_counts - 1, (preds+iter_field$background)*iter_field$time)
  # Randomized PIT for discrete distributions
  pits[[i-3]] <- Fy1 + runif(length(iter_field$sim_counts)) * (Fy - Fy1)
  if (max(hist(pits[[i-3]], plot=FALSE)$density) > ymax) {
    ymax <- max(hist(pits[[i-3]], plot=FALSE)$density)
  }
}

pdf("real_pit_hists.pdf", width=13.5, height=5.625)
layout(matrix(c(1,1,2:6,1,1,7:11), nrow=2, byrow=TRUE))
par(mar=c(0.9, 1.2, 1.0, 1.8), mgp=c(1.0, 0.7, 0), oma=c(4.0, 2.75, 0.5, 0.5),
  cex.main=2.0, cex.axis=1.1)
# par(mar=c(5.1, 4.1, 2.1, 2.1), oma=c(0, 0, 0, 0),
  # mgp=c(3.0, 0.7, 0))
##### Display in histogram
hist(pit, breaks=20, main="", xlab="", ylab="",
  col="lightgray", border="white", freq=FALSE)
title(main="2009-2011", line=-1.5)
abline(h=1, col="red", lwd=2, lty=2)
### For each calibration result
par(mar=c(0.9, 1.2, 1.0, 0.3))
for (i in 1:length(pits)) {
  ##### Display in histogram
  hist(pits[[i]], breaks=20, main="", ylab="",
    xlab="", col="lightgray", border="white", freq=FALSE,
    axes=FALSE, ylim=c(0, ymax))
  title(main=res[[i+3]]$year, line=-1.5)
  abline(h=1, col="red", lwd=1, lty=2)
  if ((i-1) %% 5==0) {
    axis(2, at=seq(0, ymax, by=0.5))
  } else {
    axis(2, at=seq(0, ymax, by=0.5), labels=FALSE)
  }
  if (i > 5*(2-1)) {
    axis(1, at=seq(0, 1, by=0.25))
  } else {
    axis(1, at=seq(0, 1, by=0.25), labels=FALSE)
  }
}
mtext("Probability Integral Transform", side=1, outer=TRUE, line=2.75, cex=1.1)
mtext("Density", side=2, outer=TRUE, line=1.5, cex=1.1)
dev.off()
