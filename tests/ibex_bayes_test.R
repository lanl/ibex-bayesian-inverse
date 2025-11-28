library(doParallel)
library(parallel)

###############################################################################
#### Solves a Bayesian inverse problem for particle counts collected by
#### the Interstellar Boundary Explorer satellite. Uses the Scaled Vecchia GP
#### approximation to fit surrogate model from the computer model output.
#### DATA NEEDED: sims.csv, ibex_real.csv, synth_sat_data.csv
###############################################################################

setwd("..")
source("pois_bayes_inv.R")
source("helper.R")
setwd("tests")

###############################################################################

args <- R.utils::commandArgs(asValues=TRUE)

## flag indicating if field data in inverse problem should be real or synthetic
real <- ifelse(!is.null(args[["r"]]), as.logical(args[["r"]]), FALSE)
## if using real field data, what year it should come from
year <- args[["y"]]
## file that contains multiple parameter combinations for synthetic field data
infile <- args[["if"]]
## flag to print more output to screen
vb <- ifelse(!is.null(args[["v"]]), as.logical(args[["v"]]), FALSE)
settings <- list(real=real, year=year, infile=infile, vb=vb)
if (vb) print(settings)

model_data <- read.csv(file="../data/sims.csv")
if (real) {
  field_data <- read.csv(file="../data/ibex_real.csv")
  if (is.null(year)) {
    stop("must specify a year")
  }
  if (year=="all") {
    cpars <- matrix(2009:2022, ncol=1)
  } else if (year=="mod_align") {
    cpars <- matrix(2009)
  } else {
    cpars <- matrix(year)
  }
} else {
  field_data <- read.csv(file="../data/synth_sat_data.csv")
  if (is.null(infile)) {
    stop("must specify a file with unique model parameter combinations")
  }
  cpars <- read.csv(file=infile, head=TRUE)
}

ncores <- detectCores()-1
cl <- parallel::makeCluster(ifelse(nrow(cpars) < ncores, nrow(cpars), ncores),
  outfile="log.txt")
doParallel::registerDoParallel(cl)
foreach(i = 1:nrow(cpars), .packages=c("GpGp", "GPvecchia", "laGP", "tmvtnorm")) %dopar% {
  if (real) {
    if (year=="mod_align") {
      map <- paste0(2009:2011, "A")
    } else {
      map <- paste0(cpars[i], "A")
    }
  } else {
    thread_fd <- field_data[field_data$parallel_mean_free_path==cpars[i,c("pmfp")] &
      field_data$ratio==cpars[i,c("ratio"),],]
    thread_fd$counts <- thread_fd$sim_counts
    thread_fd$ecliptic_lon <- thread_fd$lon
    thread_fd$ecliptic_lat <- thread_fd$lat
    map <- NULL
  }
  pd <- preprocess_data(md=model_data, fd=thread_fd, map=map)
  res <- pois_bayes_inv(xm=pd$xm, um=pd$um, ym=pd$ym, xf=pd$xf, yf=pd$yf, e=pd$e,
    lam0=pd$bg, T=10000)
  if (real) {
    res$year <- map
  } else {
    res$truth <- cpars[i,]
  }
  saveRDS(res, file=paste0("pois_bayes_inv_res_",
    ifelse(real, map, paste0("pmfp", cpars[i,c("pmfp")], "_rat", cpars[i,c("ratio")])),
    format(Sys.time(), "_%Y%m%d%H%M%S"), ".rds"))
}
parallel::stopCluster(cl)
