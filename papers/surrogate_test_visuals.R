res <- readRDS("surrogate_test_202050801.rds")

par(mfrow=c(1,1), mar=c(5.1, 4.1, 0.2, 0.2))
pdf("ibex_surr_rmse1.pdf", width=4, height=5)
rmse1 <- res$rmse[,c(1,5,6,7)]
labels <- c("SVEC (m=25)", "laGP", "deepgp", "SEPIA (pc=3)")
boxplot(rmse1, ylab="rmse", xaxt="n")
axis(1, at=1:ncol(rmse1), labels = FALSE)
text(x=1:ncol(rmse1), y=par("usr")[3]-0.0003, labels=labels, srt=45,
  xpd=TRUE, adj=1, cex=0.9)
dev.off()

par(mfrow=c(1,1), mar=c(5.1, 4.1, 0.2, 0.2))
pdf("ibex_surr_rmse2.pdf", width=4, height=5)
rmse2 <- res$rmse[,c(2:4, 8:10)]
labels <- c("SVEC (m=50)", "SVEC (m=75)", "SVEC (m=100)", "SEPIA (pc=4)", "SEPIA (pc=5)", "SEPIA (pc=6)")
boxplot(rmse2, ylab="rmse", xaxt="n")
axis(1, at=1:ncol(rmse2), labels = FALSE)
text(x=1:ncol(rmse2), y=par("usr")[3]-0.00003, labels=labels, srt=45,
  xpd=TRUE, adj=1, cex=0.9)
dev.off()

par(mfrow=c(1,1), mar=c(5.1, 4.1, 0.2, 0.2))
pdf("ibex_surr_rmse1.pdf", width=4, height=5)
crps1 <- res$crps[,c(1,5,6,7)]
labels <- c("SVEC (m=25)", "laGP", "deepgp", "SEPIA (pc=3)")
boxplot(crps1, ylab="crps", xaxt="n")
axis(1, at=1:ncol(crps1), labels = FALSE)
text(x=1:ncol(crps1), y=par("usr")[3]-0.00015, labels=labels, srt=45,
  xpd=TRUE, adj=1, cex=0.9)
dev.off()

par(mfrow=c(1,1), mar=c(5.1, 4.1, 0.2, 0.2))
pdf("ibex_surr_crps2.pdf", width=4, height=5)
crps2 <- res$crps[,c(2:4, 8:10)]
labels <- c("SVEC (m=50)", "SVEC (m=75)", "SVEC (m=100)", "SEPIA (pc=4)", "SEPIA (pc=5)", "SEPIA (pc=6)")
boxplot(crps2, ylab="crps", xaxt="n")
axis(1, at=1:ncol(crps2), labels = FALSE)
text(x=1:ncol(crps2), y=par("usr")[3]-0.000015, labels=labels, srt=45,
  xpd=TRUE, adj=1, cex=0.9)
dev.off()
