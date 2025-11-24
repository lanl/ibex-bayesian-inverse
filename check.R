check_settings <- function(settings) {

  if (is.null(settings$step)) {
    settings$step <- 0.05
  }

  if (is.null(settings$beta)) {
    settings$beta <- 1.1
  }

  if (is.null(settings$m)) {
    settings$m <- 25
  }

  if (is.null(settings$sample_scl)) {
    settings$sample_scl <- FALSE
  }

  return(settings)
}

check_init <- function(init, um) {

  if (is.null(init$u)) {
    init$u <- apply(um, 2, mean)
  }
  return(init)
}
