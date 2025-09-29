library(ggplot2)

source("../helper.R")

map <- "2020A"
ibex_real <- read.csv("../data/ibex_real.csv")
ibex_real <- ibex_real[ibex_real$map==map,]
ibex_real$est_rate <- ibex_real$counts / ibex_real$time - ibex_real$background

ibex_sim <- read.csv("../data/sims.csv")
ibex_sim1 <- ibex_sim[ibex_sim$parallel_mean_free_path==500 &
  ibex_sim$ratio==0.001,]
ibex_sim2 <- ibex_sim[ibex_sim$parallel_mean_free_path==3000 &
  ibex_sim$ratio==0.001,]

# ena_range <- c(0.04174779, 0.18489323)
ena_range <- range(ibex_sim$blurred_ena_rate)

ibex_sim$nlon <- nose_center_lons(ibex_sim$lon)
ibex_real$nlon <- nose_center_lons(ibex_real$ecliptic_lon)
ibex_sim1$nlon <- nose_center_lons(ibex_sim1$lon)
ibex_sim2$nlon <- nose_center_lons(ibex_sim2$lon)

cols <- colorRampPalette(c("blue", "cyan", "green", "yellow", "red", "magenta"))(500)
bks <- seq(ena_range[1], ena_range[2], length=length(cols)+1)
ylims <- range(ibex_sim$lat)
xlims <- rev(range(ibex_sim$nlon))

ibex_real_lons <- sort(unique(ibex_real$nlon))
ibex_real_lats <- sort(unique(ibex_real$lat))
ibex_real_rates <- cut(ibex_real$est_rate, breaks=bks, labels=FALSE)
ibex_real_rates[which(ibex_real$est_rate <= ena_range[1])] <- 1
ibex_real_rates[which(ibex_real$est_rate >= ena_range[2])] <- length(cols)
ibex_real_cols <- cols[ibex_real_rates]
pdf("ibex_real1.pdf", width=5, height=3.25)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 0.2, 6))
plot(x=ibex_real$nlon, y=ibex_real$ecliptic_lat, col=ibex_real_cols, pch=16, cex=0.7,
  xlab="Longitude", xaxt="n", ylab="Latitude", xlim=xlims, ylim=ylims,
  cex.lab=1.1)
axis(1, at=seq(325, 25, by=-60),
  labels=c(60, 0, 300, 240, 180, 120))
fields::image.plot(zlim=ena_range, col=cols, legend.lab="ENAs/sec", legend.line=3,
  legend.only=TRUE, side=4, line=2, smallplot=c(0.8, 0.84, 0.4, 0.9))
dev.off()

ibex_sim1_lons <- sort(unique(ibex_sim1$nlon))
ibex_sim1_lats <- sort(unique(ibex_sim1$lat))
ibex_sim1_zmat <- xtabs(blurred_ena_rate ~ nlon + lat, data=ibex_sim1)
ibex_sim1_zmat[ibex_sim1_zmat < ena_range[1]] <- ena_range[1]
ibex_sim1_zmat[ibex_sim1_zmat > ena_range[2]] <- ena_range[2]
pdf("sim1.pdf", width=5, height=3.25)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 0.2, 6))
image(x=ibex_sim1_lons, y=ibex_sim1_lats, z=ibex_sim1_zmat, col=cols, xlab="Longitude",
  xaxt="n", ylab="Latitude", breaks=bks, cex.lab=1.1, ylim=ylims, xlim=xlims)
axis(1, at=seq(325, 25, by=-60),
  labels=c(60, 0, 300, 240, 180, 120))
dev.off()

ibex_sim2_lons <- sort(unique(ibex_sim2$nlon))
ibex_sim2_lats <- sort(unique(ibex_sim2$lat))
ibex_sim2_zmat <- xtabs(blurred_ena_rate ~ nlon + lat, data=ibex_sim2)
ibex_sim2_zmat[ibex_sim2_zmat < ena_range[1]] <- ena_range[1]
ibex_sim2_zmat[ibex_sim2_zmat > ena_range[2]] <- ena_range[2]
pdf("sim2.pdf", width=5, height=3.25)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 0.2, 6))
image(x=ibex_sim2_lons, y=ibex_sim2_lats, z=ibex_sim2_zmat, col=cols, xlab="Longitude",
  xaxt="n", ylab="Latitude", breaks=bks, cex.lab=1.1, ylim=ylims, xlim=xlims)
axis(1, at=seq(325, 25, by=-60),
  labels=c(60, 0, 300, 240, 180, 120))
fields::image.plot(zlim=ena_range, col=cols, legend.lab="ENAs/sec", legend.line=3,
  legend.only=TRUE, side=4, line=2, smallplot=c(0.8, 0.84, 0.4, 0.9))
dev.off()
