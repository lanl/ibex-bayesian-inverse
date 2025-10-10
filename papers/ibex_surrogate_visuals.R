library(deepgp)
library(ggplot2)
library(laGP)

source("../helper.R")
source('../vecchia_scaled.R')

model_data <- read.csv(file="../data/sims.csv")
model_data <- model_data[order(model_data$parallel_mean_free_path, model_data$ratio,
  model_data$lat, model_data$lon),]
model_data[,c("x", "y", "z")] <- geo_to_spher_coords(lat=model_data$lat,
  lon=model_data$lon)
model_data_ll <- model_data[,c("lat", "lon")]
model_data_ll$nlon <- nose_center_lons(model_data_ll$lon)
model_data <- model_data[,c("x", "y", "z", "parallel_mean_free_path", "ratio",
 "blurred_ena_rate")]
# predrange <- c(0.04174779, 0.18489323)
predrange <- range(model_data$blurred_ena_rate)

md_ranges <- data.frame(matrix(NA, nrow=2, ncol=ncol(model_data)-1))
for (i in 1:(ncol(model_data)-1)) {
  md_ranges[,i] <- range(model_data[,i])
  model_data[,i] <- (model_data[,i] - md_ranges[1,i])/diff(md_ranges[,i])
}
colnames(md_ranges) <- colnames(model_data)[1:ncol(md_ranges)]
model_data <- cbind(model_data, model_data_ll)

## Hold one out - 66 model runs
pmfps <- unique(model_data$parallel_mean_free_path)
ratios <- unique(model_data$ratio)
unique_runs <- as.matrix(expand.grid(pmfps, ratios))
colnames(unique_runs) <- NULL

pmfp_unit <- pmfps[ceiling(length(pmfps)/2)]
ratio_unit <- ratios[floor(length(ratios)/2)]

Xtrain <- model_data[,c("parallel_mean_free_path", "ratio", "x", "y", "z")]
Ytrain <- model_data[,c("blurred_ena_rate")]

## takes 2-3 minutes
set.seed(2349837)
svecfit <- fit_scaled(y=Ytrain, inputs=as.matrix(Xtrain), nug=1e-4, ms=100)
all_inputs <- data.frame(svecfit$inputs.ord)
all_inputs$x <- all_inputs$x * diff(md_ranges$x) + md_ranges$x[1]
all_inputs$y <- all_inputs$y * diff(md_ranges$y) + md_ranges$y[1]
all_inputs$z <- all_inputs$z * diff(md_ranges$z) + md_ranges$z[1]
all_inputs[,c("lat", "lon")] <-
  spher_to_geo_coords(x=all_inputs$x, y=all_inputs$y, z=all_inputs$z)
all_inputs$nlon <- nose_center_lons(all_inputs$lon)
all_inputs$parallel_mean_free_path <- all_inputs$parallel_mean_free_path * diff(md_ranges$parallel_mean_free_path) +
   md_ranges$parallel_mean_free_path[1]
all_inputs$ratio <- all_inputs$ratio * diff(md_ranges$ratio) +
   md_ranges$ratio[1]

pmfp <- pmfp_unit * diff(md_ranges$parallel_mean_free_path) +
   md_ranges$parallel_mean_free_path[1]
ratio <- ratio_unit * diff(md_ranges$ratio) + md_ranges$ratio[1]

ribbon_points <- which(all_inputs$lat > -30 & all_inputs$lat < 30 &
  all_inputs$lon < 300 & all_inputs$lon > 270 &
  all_inputs$parallel_mean_free_path==pmfp & all_inputs$ratio==ratio)
## Should be index 4473
ref_point <- ribbon_points[length(ribbon_points)]

plot_data <- model_data[model_data$parallel_mean_free_path==pmfp_unit &
  model_data$ratio==ratio_unit,]
ref_neighbors <- all_inputs[svecfit$NNarray[ref_point,-1],]

cols <- colorRampPalette(c("blue", "cyan", "green", "yellow", "red", "magenta"))(500)
bks <- seq(predrange[1], predrange[2], length=length(cols)+1)
ylims <- range(model_data$lat)
xlims <- rev(range(model_data$nlon))

lons <- sort(unique(plot_data$nlon))
lats <- sort(unique(plot_data$lat))
zmat <- xtabs(blurred_ena_rate ~ nlon + lat, data=plot_data)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 0.2, 0.2))
pdf("ibex_nbr_latlon.pdf", width=7, height=5)
image(x=lons, y=lats, z=zmat, col=cols, xlab="Longitude", xaxt="n",
  ylab="Latitude", breaks=bks, xlim=rev(range(lons)), useRaster=TRUE)
axis(1, at=seq(325, 25, by=-60),
  labels=c(60, 0, 300, 240, 180, 120))
usr <- par("usr")  # plotting region: c(xmin, xmax, ymin, ymax)
rect(usr[1], usr[3], usr[2], usr[4],
     col = rgb(0.5, 0.5, 0.5, 0.8), border = NA)
points(lat ~ nlon, data=ref_neighbors[1:25,], pch=21, col=1, bg=3)
points(lat ~ nlon, data=ref_neighbors[26:50,], pch=22, col=1, bg=4)
points(lat ~ nlon, data=ref_neighbors[51:75,], pch=23, col=1, bg=5)
points(lat ~ nlon, data=ref_neighbors[76:100,], pch=24, col=1, bg=6)
points(lat ~ nlon, data=all_inputs[ref_point,], pch=8, col=2, cex=2, lwd=3)
dev.off()

par(mfrow=c(1,1), mar=c(5.1, 4.1, 0.2, 0.2))
pdf("ibex_nbr_params.pdf", width=7, height=5)
plot(x=jitter(ref_neighbors[1:25,c("parallel_mean_free_path")]),
  y=jitter(ref_neighbors[1:25,c("ratio")]), type="n",
  xlim=range(all_inputs$parallel_mean_free_path),
  ylim=range(all_inputs$ratio), pch=21, col=1, bg=3,
  xlab="parallel mean free path", ylab="ratio")
abline(v=seq(500, 3000, by=500), col=1, lty=2)
abline(h=c(0.001, 0.005, 0.01, 0.02, 0.05, 0.1), col=1, lty=2)
points(x=jitter(ref_neighbors[1:25,c("parallel_mean_free_path")]),
  y=jitter(ref_neighbors[1:25,c("ratio")]), pch=21, col=1, bg=3)
points(x=jitter(ref_neighbors[26:50,c("parallel_mean_free_path")]),
  y=jitter(ref_neighbors[26:50,c("ratio")]), pch=22, col=1, bg=4)
points(x=jitter(ref_neighbors[51:75,c("parallel_mean_free_path")]),
  y=jitter(ref_neighbors[51:75,c("ratio")]), pch=23, col=1, bg=5)
points(x=jitter(ref_neighbors[76:100,c("parallel_mean_free_path")]),
  y=jitter(ref_neighbors[76:100,c("ratio")]), pch=24, col=1, bg=6)
points(x=all_inputs[ref_point,c("parallel_mean_free_path")],
  y=all_inputs[ref_point,c("ratio")], pch=8, col=2, cex=2, lwd=3)
legend("topleft", c("point of interest", paste0("m=", c(25,50,75,100))),
  pch=c(8, 21:24), col=c(2, rep(1, 4)), pt.bg=c(NA, 3:6),
  lwd=c(2, NA, NA, NA, NA), lty=rep(NA, 5), bg="white", cex=1.05)
dev.off()


###############################################################################
###############################################################################
## Illustration showing simulated output next to surrogate output to
## demonstrate effectiveness of surrogate
###############################################################################
###############################################################################

source("../helper.R")
source('../vecchia_scaled.R')

model_data <- read.csv(file="../data/sims.csv")
model_data$nlon <- nose_center_lons(model_data$lon)
field_data <- read.csv(file="../data/sims_real.csv")
pd <- preprocess_data(md=model_data, fd=field_data, esa_lev=4,
  fparams=c(1750, 0.02), scales=c(1, 1), tol=NA, quant=0.0,
  real=FALSE, disc=FALSE)
# predrange <- c(0.04174779, 0.18489323)
predrange <- range(model_data$blurred_ena_rate)
Xtrain <- as.matrix(cbind(pd$Xmod, pd$Umod))
fit <- fit_scaled(y=pd$Zmod, inputs=Xtrain, nug=1e-4, ms=25)

pmfps <- c(1500, 1625, 1750)
ratios <- c(0.005, 0.0075, 0.01)
grid <- as.matrix(expand.grid(pmfps, ratios))
colnames(grid) <- c("pmfp", "ratio")
model_pmfps <- unique(model_data$parallel_mean_free_path)
model_ratios <- unique(model_data$ratio)
cols <- colorRampPalette(c("blue", "cyan", "green", "yellow", "red", "magenta"))(500)
bks <- seq(predrange[1], predrange[2], length=length(cols)+1)
ylims <- range(model_data$lat)
xlims <- rev(range(model_data$nlon))

pdf("ibex_surr_vis_check.pdf", width=14, height=11)
par(mfrow=c(length(ratios), length(pmfps)),
  mar=c(0, 0, 1.25, 1.25), oma=c(5, 5, 1, 1))
for (i in 1:nrow(grid)) {
  if (grid[i,1] %in% model_pmfps && grid[i,2] %in% model_ratios) {

    iter_data <- model_data[model_data$parallel_mean_free_path==grid[i,1] &
      model_data$ratio==grid[i,2],]
    lons <- sort(unique(model_data$nlon))
    lats <- sort(unique(model_data$lat))
    zmat <- xtabs(blurred_ena_rate ~ nlon + lat, data=iter_data)
    zmat[zmat > predrange[2]] <- predrange[2]
    zmat[zmat < predrange[1]] <- predrange[1]
    image(x=lons, y=lats, z=zmat, col=cols, xlab="", xaxt="n", yaxt="n",
      ylab="", breaks=bks, ylim=ylims, xlim=xlims, useRaster=TRUE)
  } else {
    pred_params <- grid[i,]
    pred_params[1] <- (pred_params[1] - 500)/2500
    pred_params[2] <- (pred_params[2] - 0.001)/(0.1-0.001)

    XX_ll <- cbind(unique(model_data[,c("lon", "lat")]), matrix(pred_params, nrow=1))
    XX_ll[,c("x", "y", "z")] <- geo_to_spher_coords(lat=XX_ll$lat, lon=XX_ll$lon)
    XX_ll$x <- (XX_ll$x - min(XX_ll$x)) / diff(range(XX_ll$x))
    XX_ll$y <- (XX_ll$y - min(XX_ll$y)) / diff(range(XX_ll$y))
    XX_ll$z <- (XX_ll$z - min(XX_ll$z)) / diff(range(XX_ll$z))
    XX <- XX_ll[,c("x", "y", "z", "1", "2")]
    colnames(XX) <- c("x", "y", "z", "pmfp", "ratio")
    XX <- as.matrix(XX)
    lhat_curr <- predictions_scaled(fit, XX, m=25, joint=FALSE, predvar=FALSE)
    XX_ll$preds <- lhat_curr
    XX_ll$nlon <- nose_center_lons(XX_ll$lon)
    lons <- sort(unique(XX_ll$nlon))
    lats <- sort(unique(XX_ll$lat))
    zmat <- xtabs(preds ~ nlon + lat, data=XX_ll)
    zmat[zmat > predrange[2]] <- predrange[2]
    zmat[zmat < predrange[1]] <- predrange[1]
    image(x=lons, y=lats, z=zmat, col=cols, xaxt="n", yaxt="n", xlab="",
      ylab="", breaks=bks, ylim=ylims, xlim=xlims, bty="n", useRaster=TRUE)
    box(lwd=4, lty=2)
  }
  if ((i-1) %% length(pmfps)==0) {
    axis(2, at=seq(-50, 50, by=50), cex.axis=1.25)
  }
  if (i > length(pmfps)*(length(ratios)-1)) {
    axis(1, at=seq(325, 25, by=-60), cex.axis=1.25,
      labels=c(60, 0, 300, 240, 180, 120))
  }
}
mtext("Longitude", side=1, line=3.5, outer=TRUE, cex=1.25)
mtext("Latitude", side=2, line=3.5, outer=TRUE, cex=1.25)
dev.off()
