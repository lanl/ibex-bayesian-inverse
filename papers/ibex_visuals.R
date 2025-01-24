library(ggplot2)

source("../helper.R")

map <- "2020A"
ibex_real <- read.csv("../data/ibex_real.csv")
ibex_real <- ibex_real[ibex_real$map==map,]
ibex_real$ena_rate <- ibex_real$counts / ibex_real$time - ibex_real$background

ibex_sim <- read.csv("../data/sims.csv")

ibex_sim1 <- ibex_sim[ibex_sim$parallel_mean_free_path==500 &
  ibex_sim$ratio==0.001,]

ibex_sim2 <- ibex_sim[ibex_sim$parallel_mean_free_path==3000 &
  ibex_sim$ratio==0.1,]

#ena_range <- range(c(ibex_real$ena_rate, ibex_sim1$blurred_ena_rate,
#  ibex_sim2$blurred_ena_rate))
ena_range <- c(0, 0.25)

ibex_real$nlon <- nose_center_lons(ibex_real$ecliptic_lon)
ibex_sim1$nlon <- nose_center_lons(ibex_sim1$lon)
ibex_sim2$nlon <- nose_center_lons(ibex_sim2$lon)

ggplot(ibex_real, aes(x=nlon, y=ecliptic_lat)) +
  geom_point(aes(color=ena_rate), size=0.4) +
  scale_color_gradientn(colours=rev(rainbow(6)[c(6,1:5)]),
    limits=ena_range, oob=scales::oob_squish) +
  scale_x_continuous(trans="reverse", breaks=seq(325, 25, by=-60),
    label=c(60, 0, 300, 240, 180, 120)) +
  theme_bw() +
  theme(axis.title.x=element_text(margin=margin(t=10, r=0, b=0, l=0))) +
  labs(x="  ", y="Latitude")
ggsave("ibex_real1.pdf", dpi=320, width=4, height=3.25)

ggplot(ibex_sim1, aes(x=nlon, y=lat)) +
  geom_raster(aes(fill=ena_rate)) +
  scale_fill_gradientn(colours=rev(rainbow(6)[c(6,1:5)]),
    limits=ena_range, oob=scales::oob_squish) +
  scale_x_continuous(trans="reverse", breaks=seq(325, 25, by=-60),
    label=c(60, 0, 300, 240, 180, 120)) +
  theme_bw() +
  theme(axis.ticks.y=element_blank(),
        axis.title.x=element_text(margin=margin(t=10, r=0, b=0, l=0))) +
  labs(x="Longitude", y="Latitude")
ggsave("sim1.pdf", dpi=320, width=4, height=3.25)

ggplot(ibex_sim2, aes(x=nlon, y=lat)) +
  geom_raster(aes(fill=ena_rate)) +
  scale_fill_gradientn(colours=rev(rainbow(6)[c(6,1:5)]),
    limits=ena_range, oob=scales::oob_squish) +
  scale_x_continuous(trans="reverse", breaks=seq(325, 25, by=-60),
    label=c(60, 0, 300, 240, 180, 120)) +
  theme_bw() +
  theme(axis.ticks.y=element_blank(),
        axis.title.x=element_text(margin=margin(t=10, r=0, b=0, l=0)),
        legend.spacing.y=unit(5.5, 'mm')) +
  labs(x="  ", y="Latitude", fill="ENA rate\n")
ggsave("sim2.pdf", dpi=320, width=4, height=3.25)
