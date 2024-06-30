library(doParallel)
library(parallel)

source("../helper.R")
source("../mcmc.R")

###############################################################################
#### Conducts Bayesian Computer Model Calibration for the data collected by
#### the Interstellar Boundary Explorer satellite. Uses scaled Vecchia
#### approximation or locally approximated Gaussian processes to fit surrogate
#### model from the computer model output.
###############################################################################

args <- R.utils::commandArgs(asValues=TRUE)

## energy level of data to work with
esa_lev <- ifelse(!is.null(args[["esa"]]), as.integer(args[["esa"]]), 4)
## number of MCMC iterations to run
nmcmcs <- ifelse(!is.null(args[["nmcmcs"]]), as.integer(args[["nmcmcs"]]), 10000)
## specifies type of GP to fit (e.g. svecchia, laGP w/ nearest neighbors)
gp <- ifelse(!is.null(args[["gp"]]), as.character(args[["gp"]]), "svecchia")
## if using laGP, end is the neighborhood size for each local GP
end <- ifelse(!is.null(args[["end"]]), as.integer(args[["end"]]), 25)
## if using laGP, number of threads to run during each fit
thrds <- ifelse(!is.null(args[["threads"]]), as.integer(args[["threads"]]), 10)
## maximum number of calibrations to perform at once (for the dopar loop)
maxprocs <- ifelse(!is.null(args[["procs"]]), as.integer(args[["procs"]]), 6)
## flag indicating if simulated field data should have discrepancy
disc <- ifelse(!is.null(args[["d"]]), as.logical(args[["d"]]), FALSE)
## flag indicating if field data in calibration should be real or simulated
real <- ifelse(!is.null(args[["r"]]), as.logical(args[["r"]]), FALSE)
## if using simulated field data, what the true PMFP and ratio should be
fpmfp <- ifelse(!is.null(args[["fpmfp"]]), as.numeric(args[["fpmfp"]]), 2500)
fratio <- ifelse(!is.null(args[["fratio"]]), as.numeric(args[["fratio"]]), 0.05)
## if using real field data, what year it should come from
fyear <- ifelse(!is.null(args[["fyear"]]), as.numeric(args[["fyear"]]), 2009)
## file that contains multiple parameter combinations for field data
infile <- ifelse(!is.null(args[["infile"]]), as.character(args[["infile"]]), NA)
## amount the PMFP and ratio should be scaled before fitting (used for laGP)
psc <- ifelse(!is.null(args[["psc"]]), as.numeric(args[["psc"]]), 1)
rsc <- ifelse(!is.null(args[["rsc"]]), as.numeric(args[["rsc"]]), 1)
## step size for random walk in McMC
step_size <- ifelse(!is.null(args[["step"]]), as.numeric(args[["step"]]), 0.06517432)
## quantile of changes in y below which to remove points
quant <- ifelse(!is.null(args[["quant"]]), as.numeric(args[["quant"]]), 0.0)
## tolerance of changes in y below which to remove points
tol <- ifelse(!is.null(args[["tol"]]), as.numeric(args[["tol"]]), NA)
## flag to turn on "debug" mode, which will save more output
debug <- ifelse(!is.null(args[["debug"]]), as.logical(args[["debug"]]), FALSE)
## flag to print more output to screen
vb <- ifelse(!is.null(args[["v"]]), as.logical(args[["v"]]), FALSE)
settings <- list(esa_lev=esa_lev, nmcmcs=nmcmcs, gp=gp, end=end, thrds=thrds,
  maxprocs=maxprocs, real=real, fpmfp=fpmfp, fratio=fratio, fyear=fyear,
  infile=infile, psc=psc, rsc=rsc, step_size=step_size, quant=quant, tol=tol,
  debug=debug, vb=vb)
if (vb) print(settings)

model_data <- read.csv(file="../data/sims.csv")
if (!real) {
  if (disc) {
    field_data <- read.csv(file="../data/sims_bias.csv")
  } else {
    field_data <- read.csv(file="../data/sims_real.csv")
  }
} else {
  field_data <- read.csv(file="../data/ibex_real.csv")
  field_data <- field_data %>% dplyr::rename(esa=ESA, sim_counts=counts)
}

if (!is.na(infile)) {
  cpars <- read.csv(file=infile, head=TRUE)
} else if (disc) {
  cpars <- data.frame(matrix(data=unique(field_data$scl)))
  names(cpars) <- c("scl")
} else {
  if (!real) {
    cpars <- data.frame(matrix(data=c(fpmfp, fratio), nrow=1))
    names(cpars) <- c("pmfp", "ratio")
  } else {
    cpars <- data.frame(matrix(data=c(fyear)))
    names(cpars) <- c("year")
  }
}

cl <- parallel::makeCluster(ifelse(nrow(cpars) < maxprocs, nrow(cpars),
                                   maxprocs), outfile="../temp/log.txt")
doParallel::registerDoParallel(cl)
foreach(i = 1:nrow(cpars), .packages=c("GpGp", "GPvecchia", "laGP", "tidyverse",
                                       "tmvtnorm")) %dopar% {
  if (!real && !disc) {
    fparams <- c(cpars[i,]$pmfp, cpars[i,]$ratio)
  } else {
    fparams <- c(cpars[i,])
  }
  pd <- preprocess_data(md=model_data, fd=field_data, esa_lev=esa_lev,
    fparams=fparams, scales=c(psc, rsc), tol=tol, quant=quant, real=real,
    disc=disc)
  mcmc_res <- mcmc(Xm=pd$Xmod, Um=pd$Umod, Zm=pd$Zmod, Xf=pd$Xfield,
    Zf=pd$Zfield, Of=pd$Ofield, nmcmcs=nmcmcs, step=step_size,
    gpmeth=gp, end=end, thrds=thrds, vb=vb, debug=debug)
  res <- list(mcmc_res=mcmc_res,
    settings=list(esa=esa_lev, nmcmcs=nmcmcs, threads=thrds,
                  params=fparams, real=real, step=step_size,
                  pmfp_scale=psc, ratio_scale=rsc))
  saveRDS(res, file=paste0("../results/mcmc_res_esa", esa_lev, "_nmcmc", nmcmcs,
    "_end", end, "_sc", strsplit(as.character(psc), split='\\.')[[1]][2],
    "_quant", strsplit(as.character(quant), split='\\.')[[1]][2],
    ifelse(gp=="svecchia", "_svecchia", ""),
    ifelse(real, paste0("_real", fparams[1]),
           paste0("_pmfp", fparams[1], "_rat",
                  strsplit(as.character(fparams[2]), split='\\.')[[1]][2])),
    format(Sys.time(), "_%Y%m%d%H%M%S"), ".rds"))
}
parallel::stopCluster(cl)
