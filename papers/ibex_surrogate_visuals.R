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
ref_point <- ribbon_points[length(ribbon_points)]

plot_data <- model_data[model_data$parallel_mean_free_path==pmfp_unit &
  model_data$ratio==ratio_unit,]
ref_neighbors <- all_inputs[svecfit$NNarray[ref_point,-1],]

lons <- sort(unique(plot_data$nlon))
lats <- sort(unique(plot_data$lat))
zmat <- xtabs(blurred_ena_rate ~ nlon + lat, data=plot_data)
cols <- colorRampPalette(c("blue", "cyan", "green", "yellow", "red", "magenta"))
par(mfrow=c(1,1), mar=c(5.1, 4.1, 0.2, 0.2))
pdf("ibex_nbr_latlon.pdf", width=7, height=5)
image(x=lons, y=lats, z=zmat, col=cols(500), xlab="Longitude", xaxt="n",
  ylab="Latitude", xlim=rev(range(lons)))
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
  y=jitter(ref_neighbors[1:25,c("ratio")]), xlim=range(all_inputs$parallel_mean_free_path),
  ylim=range(all_inputs$ratio), pch=21, col=1, bg=3,
  xlab="parallel mean free path", ylab="ratio")
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
