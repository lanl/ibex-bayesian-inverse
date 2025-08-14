## Visuals for varying the dimension of the response
surr_times <- readRDS("surrogate_time_test_20250811.rds")
surr_times$fit_times <- surr_times$fit_times[,c(1:38, 40:dim(surr_times$fit_times)[2]),]
surr_times$pred_times <- surr_times$pred_times[,c(1:38, 40:dim(surr_times$pred_times)[2]),]
surr_fit_times <- apply(surr_times$fit_times, c(2,3), mean)
surr_pred_times <- apply(surr_times$pred_times, c(2,3), mean)
sepia_fit_files <- list.files(pattern="sepia_fit_times_[3-6]_dim.csv")
sepia_pred_files <- list.files(pattern="sepia_pred_times_[3-6]_dim.csv")

for (i in 1:length(sepia_fit_files)) {
  sep_fit <- read.csv(sepia_fit_files[i], header=FALSE)
  sep_pred <- read.csv(sepia_pred_files[i], header=FALSE)
  surr_fit_times <- cbind(surr_fit_times, apply(sep_fit, 2, mean))
  surr_pred_times <- cbind(surr_pred_times, apply(sep_pred, 2, mean))
  colnames(surr_fit_times)[ncol(surr_fit_times)] <-
    colnames(surr_pred_times)[ncol(surr_pred_times)] <- paste0("sepia", i+2)
}

surr_total_times <- surr_fit_times + surr_pred_times

par(mfrow=c(1,1), mar=c(5.1, 4.1, 0.2, 0.2))
pdf("ibex_surr_fit_times.pdf", width=4, height=5)
fit_times_ord <- surr_fit_times[,c(1,6,7)]
exp_pows <- 7:44
matplot(x=exp_pows, y=fit_times_ord[1:length(exp_pows),]/60, type="l", ylim=c(0, 30),
  xlab="dim of response = 10 + 1.25^x", ylab="fitting time (minutes)",
  lwd=3)
legend("topleft", c("SVEC (m=25)", "deepgp", "SEPIA (pc=3)", "laGP"),
  col=1:4, lty=1:4, lwd=2, cex=0.75)
dev.off()

par(mfrow=c(1,1), mar=c(5.1, 4.1, 0.2, 0.2))
pdf("ibex_surr_pred_times.pdf", width=4, height=5)
pred_times_ord <- surr_pred_times[,c(1,6,7,5)]
matplot(x=exp_pows, y=pred_times_ord[1:length(exp_pows),]/60, type="l", ylim=c(0, 10),
  xlab="dim of response = 10 + 1.25^x", ylab="prediction time (minutes)",
  lwd=3)
dev.off()

large_ns <- seq(20000, 75000, by=5000)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 0.2, 0.2))
pdf("ibex_surr_total_times.pdf", width=4, height=5)
long_times_ord <- surr_total_times[,c(1:4,7:10)]
matplot(x=large_ns, y=long_times_ord[(length(exp_pows)+1):nrow(long_times_ord),]/60, type="l",
  xlab="dim of response", ylab="fitting + prediction time (minutes)",
  lwd=3, col=rep(c(1,3), each=4), lty=rep(1:4, 2))
legend("topleft", c("SVEC (m=25,50,75,100)", "SEPIA (pc=3,4,5,6)"), col=c(1,3), lty=1, lwd=2, cex=0.75)
dev.off()

## Visuals for varying the number of computer experiment runs
surr_times <- readRDS("surrogate_time_test_n_20250806.rds")
surr_fit_times <- apply(surr_times$fit_times, c(2,3), mean)
surr_pred_times <- apply(surr_times$pred_times, c(2,3), mean)
sepia_fit_files <- list.files(pattern="sepia_fit_times_3_ns.csv")
sepia_pred_files <- list.files(pattern="sepia_pred_times_3_ns.csv")

for (i in 1:length(sepia_fit_files)) {
  sep_fit <- read.csv(sepia_fit_files[i], header=FALSE)
  sep_pred <- read.csv(sepia_pred_files[i], header=FALSE)
  surr_fit_times <- cbind(surr_fit_times, apply(sep_fit, 2, mean))
  surr_pred_times <- cbind(surr_pred_times, apply(sep_pred, 2, mean))
  colnames(surr_fit_times)[ncol(surr_fit_times)] <-
    colnames(surr_pred_times)[ncol(surr_pred_times)] <- paste0("sepia", i+2)
}

surr_total_times <- surr_fit_times + surr_pred_times

par(mfrow=c(1,1), mar=c(5.1, 4.1, 0.2, 0.2))
pdf("ibex_surr_fit_times_ns.pdf", width=4, height=5)
fit_times_ord <- surr_fit_times[,c(1,6,7)]
ns <- seq(10, 100, by=10)
matplot(x=ns, y=fit_times_ord[1:length(ns),]/60, type="l", ylim=c(0, 30),
  xlab="# simulator runs", ylab="fitting time (minutes)",
  lwd=3)
legend("topright", c("SVEC (m=25)", "deepgp", "SEPIA (pc=3)", "laGP"),
  col=1:4, lty=1:4, lwd=2, cex=0.75)
dev.off()

par(mfrow=c(1,1), mar=c(5.1, 4.1, 0.2, 0.2))
pdf("ibex_surr_pred_times_ns.pdf", width=4, height=5)
pred_times_ord <- surr_pred_times[,c(1,6,7,5)]
matplot(x=ns, y=pred_times_ord[1:length(ns),]/60, type="l", ylim=c(0, 10),
  xlab="# simulator runs", ylab="prediction time (minutes)",
  lwd=3)
dev.off()

match_cols <- c(1,3:4)
large_ns <- seq(500, 2500, by=500)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 0.2, 0.2))
pdf("ibex_surr_total_times_ns.pdf", width=4, height=5)
long_times_ord <- surr_total_times[,c(1,5,7)]
matplot(x=large_ns, y=long_times_ord[(length(ns)+1):nrow(long_times_ord),]/60, type="l",
  xlab="# simulator runs", ylab="fitting + prediction time (minutes)",
  lwd=3, col=match_cols, lty=match_cols, ylim=c(0, 30))
dev.off()
