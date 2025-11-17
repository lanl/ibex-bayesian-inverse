###############################################################################
## FIGURE 6: Performance metrics for surrogate modeling of the IBEX simulator
## DATA NEEDED: surrogate_test.rds, sepia_metrics_X.csv
###############################################################################

surr1_res <- readRDS("surrogate_test.rds")
sepia_files <- list.files(pattern="sepia_metrics_[3-6].csv")

rmses <- surr1_res$rmse
crps <- surr1_res$crps
for (i in 1:length(sepia_files)) {
  sep_res <- read.csv(sepia_files[i], header=FALSE)
  rmses <- cbind(rmses, sep_res[,1])
  crps <- cbind(crps, sep_res[,2])
  colnames(crps)[ncol(crps)] <- colnames(rmses)[ncol(rmses)] <- paste0("sepia", i+2)
}

## Figure 6 (left panel)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 0.2, 0.2))
pdf("ibex_surr_rmse.pdf", width=7, height=5)
rmses_ord <- rmses[,c(1,5:7,2:4,8:10)]
labels <- c("SVEC (m=25)", "laGP", "deepgp", expression(SEPIA~"("*n[k]*"="*3*")"),
  paste0("SVEC (m=", c(50, 75, 100), ")"), expression(SEPIA~"("*n[k]*"="*4*")"),
  expression(SEPIA~"("*n[k]*"="*5*")"), expression(SEPIA~"("*n[k]*"="*6*")"))
boxplot(rmses_ord, ylab="rmse", xaxt="n", cex.lab=1.05)
axis(1, at=1:ncol(rmses_ord), labels=FALSE, tck=-0.02)
text(x=1:ncol(rmses_ord), y=par("usr")[3]-0.0002, labels=labels, srt=45,
  xpd=TRUE, adj=1, cex=1.05)
abline(v=4.5, lty=2)
dev.off()

## Figure 6 (right panel)
par(mfrow=c(1,1), mar=c(5.1, 4.3, 0.2, 0.2))
pdf("ibex_surr_crps.pdf", width=7, height=5)
crps_ord <- crps[,c(1,5:7,2:4,8:10)]
boxplot(crps_ord, ylab="crps", xaxt="n", cex.lab=1.05)
axis(1, at=1:ncol(crps_ord), labels = FALSE, tck=-0.02)
text(x=1:ncol(crps_ord), y=par("usr")[3]-0.0002, labels=labels, srt=45,
  xpd=TRUE, adj=1, cex=1.05)
abline(v=4.5, lty=2)
dev.off()
