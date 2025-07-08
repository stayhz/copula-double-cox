#----------------------------------------------------------------------
# Simulation under dependent censoring data.
# Then fit data using the copula-single-Cox model.
#
# Arguments:
#   vv          : iteration index (for seed shifting)
#   n           : sample size
#   p, q        : number of covariates in X and W
#   truepara    : true parameter vector
#   copula_name : copula name vector (copula_name[2] used)
#   dist        : distribution ("Weibull", "GenGompertz", or "GenExp")
#   numseed0    : base seed number
#   pA          : parameter for administrative censoring mechanism
#
# Returns:
#   A list with parameter estimates, SE, data, time, and retry_count
#----------------------------------------------------------------------

simu_single <- function(vv, n, p, q, truepara, copula_name, dist, numseed0, pA) {
  start_time <- Sys.time()
  
  ## === Load libraries === ##
  library(copula)
  library(ucminf)
  library(parallel)
  library(MASS)
  
  ## === Source all required helper functions === ##
  source("func_marg.R")
  source("func_partialcop.R")
  source("likelihood.R")
  source("copulacdf2.R")
  source("func_revkendall.R")
  source("gendata_weibull.R")
  source("gendata_genexp.R")
  source("easydoublecox.R")
  source("method_doublecox.R")
  source("method_indep.R")
  source("method_singlecox.R")
  
  ## === Retry mechanism in case of estimation failure === ##
  retry_count <- 0
  success <- FALSE
  
  while (!success && retry_count <= 5) {
    try({
      numseed <- numseed0 + vv + retry_count
      
      # === Generate data based on specified distribution === #
      gendata <- switch(dist,
                        "Weibull"     = gendata_weibull(numseed, num1, n, p, q, truepara, copula_name, dist, pA),
                        "GenExp"      = gendata_genexp(numseed, num1, n, p, q, truepara, copula_name, dist, pA),
                        stop("Unsupported distribution type")
      )
      
      # === Extract variables === #
      covaX <- gendata$covaX
      covaW <- gendata$covaW
      trunc <- gendata$trunc
      time  <- gendata$time
      Delta <- gendata$Delta
      Xi    <- gendata$Xi
      
      # === Construct data frame for estimation === #
      Ddata <- data.frame(
        trunc = trunc,
        time = time,
        covaX1 = covaX[, 1],
        covaX2 = covaX[, 2],
        covaW1 = covaW[, 1],
        covaW2 = covaW[, 2],
        Delta = Delta,
        Xi = Xi
      )
      
      # === Fit copula-single-Cox model === #
      estresult <- method_singlecox(n, p, q, copula_name, dist, numseed0, Ddata)
      estpara   <- estresult$params
      SE        <- estresult$SE
      likvalue  <- -estresult$value  # Likelihood (negated for maximization)
      
      success <- TRUE
    }, silent = TRUE)
    retry_count <- retry_count + 1
  }
  
  etime <- Sys.time() - start_time
  
  ## === Return result list === ##
  if (success) {
    print(vv)  # Print progress
    return(list(
      params = estpara,
      SE = SE,
      Z = time,
      X1 = covaX[, 1],
      X2 = covaX[, 2],
      Delta = Delta,
      Xi = Xi,
      etime = etime,
      retry_count = retry_count,
      likvalue = likvalue
    ))
  } else {
    return(list(
      params = NA,
      SE = NA,
      Z = time,
      X1 = covaX[, 1],
      X2 = covaX[, 2],
      Delta = Delta,
      Xi = Xi,
      etime = etime,
      retry_count = retry_count,
      likvalue = NA
    ))
  }
}
