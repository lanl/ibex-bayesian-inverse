

###############################################################################
#### Conducts Bayesian Computer Model Calibration for the data collected by
#### the Interstellar Boundary Explorer satellite. Uses field data from
#### 2009-2011.
###############################################################################

source("../helper.R")
source("../mcmc.R")

index <- NA
ncvs <- NA
scale <- NA
seed <- 7198701
map <- "2020A"

scales <- c(seq(0.2, 1.0, by=0.2), 1.5, 2, 4, 8, 10)

## read in the command line arguments
## run with: R CMD BATCH '--args index=1 scale=0.2 ncvs=10 seed=719 map="2020A"' scale_disc_test.R
args <- commandArgs(TRUE)
if (length(args) > 0) {
  for(i in 1:length(args)) {
    eval(parse(text=args[[i]]))
  }
}
print(args)

if (!is.na(index) && !is.na(ncvs)) {
  design <- rep(scales, ncvs)
  scale <- design[index] 
}
settings <- data.frame(index=index, ncvs=ncvs, scale=scale, seed=seed)
print(settings)

## Load model and field data
model_data <- read.csv(file="../data/sims.csv")
sub_model_data <- model_data[model_data$parallel_mean_free_path==1750 & model_data$ratio==0.02,]
field_data <- read.csv(file="../data/ibex_real.csv")
colnames(field_data)[colnames(field_data)== "counts"] <- "sim_counts"
colnames(field_data)[colnames(field_data)== "ecliptic_lat"] <- "lat"
colnames(field_data)[colnames(field_data)== "ecliptic_lon"] <- "lon"
field_data <- field_data[field_data$map==map,]

field_sim_map <- read.csv(paste0("../data/ibex_sim_map_", map, ".csv"))
field_sim_map <- field_sim_map[order(field_sim_map$row),]
fd_true_rate <- scale*sub_model_data[field_sim_map$near_sim,c("blurred_ena_rate")]
rates_wbg_exp <- (fd_true_rate + field_data$background)*field_data$time
field_data$sim_counts <- rpois(rates_wbg_exp, rates_wbg_exp)

## Create simulated field data with multiplicative discrepancy
pd <- preprocess_data(md=model_data, fd=field_data, esa_lev=4,
  fparams=map, scales=c(1, 1), tol=NA, quant=0.0, real=TRUE)
mcmc_res <- mcmc(Xm=pd$Xmod, Um=pd$Umod, Zm=pd$Zmod, Xf=pd$Xfield,
  Zf=pd$Zfield, Of=pd$Ofield, m=25, nmcmcs=50, step=0.05,
  gpmeth="svecchia", vb=TRUE, true_u=NA, true_logscl=NA,
  betashape=2, adapt=FALSE)
res <- list(mcmc_res=mcmc_res, settings=settings)

out_fn_sub <- ifelse(!is.na(index), paste0("index_", index), paste0("scale_", scale))
saveRDS(res, file=paste0("../results/scale_disc_test_", out_fn_sub,
  format(Sys.time(), "_%Y%m%d%H%M%S"), ".rds"))
