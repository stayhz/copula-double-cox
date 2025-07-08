#----------------------------------------------------------------------
# Simulate data under the dependent model and fit using double-Cox method
# Generate data using log-normal baseline
#----------------------------------------------------------------------

simu_lognorm <- function(vv, n, p, q, truepara, copula_name, dist, numseed0, pA) {
  start_time <- Sys.time()
  
  # Load required libraries
  library(copula)
  library(ucminf)
  library(parallel)
  library(MASS)
  
  # Load source code for functions
  source("func_marg.R")
  source("func_partialcop.R")
  source("likelihood.R")
  source("copulacdf2.R")
  source("func_revkendall.R")
  source("gendata_weibull.R")
  source("gendata_genexp.R")
  source("easydoublecox.R")
  source("method_doublecox.R")
  source("gendata_lognorm.R")
  
  # Initialization
  retry_count <- 0
  success <- FALSE
  
   while (!success && retry_count <= 5) {
    try({
      numseed <- numseed0 + vv + retry_count
      
      # Generate data using log-normal baseline
      gendata <- gendata_lognorm(numseed, num1, n, p, q, truepara, copula_name, dist, pA)
      
      # Extract generated data
      covaX  <- gendata$covaX
      covaW  <- gendata$covaW
      trunc  <- gendata$trunc
      time   <- gendata$time
      Delta  <- gendata$Delta
      Xi     <- gendata$Xi
      
      # Extract individual covariates
      covaX1 <- covaX[, 1]
      covaX2 <- covaX[, 2]
      covaW1 <- covaW[, 1]
      covaW2 <- covaW[, 2]
      
      # Combine into a data frame
      Ddata <- data.frame(trunc, time, covaX1, covaX2, covaW1, covaW2, Delta, Xi)
      
      # Inspect censoring/truncation rates
      prop_Delta     <- mean(Delta == 1)
      prop_Xi        <- mean(Xi == 1)
      prop_both_zero <- mean((Delta + Xi) == 0)
      
      # Fit the double-Cox model
      estresult <- method_doublecox(n, p, q, copula_name, dist, numseed0, Ddata)
      
      estpara   <- estresult$params
      SE        <- estresult$SE
      likvalue  <- -estresult$value  # Flip sign of log-likelihood value
      
      success <- TRUE
    }, silent = TRUE)
    
    retry_count <- retry_count + 1
  }
  
  etime <- Sys.time() - start_time
  
  # Return simulation result
  if (success) {
    cat("Simulation ID:", vv, "\n")
    return(list(
      params      = estpara,
      Z           = time,
      X1          = covaX1,
      X2          = covaX2,
      W1          = covaW1,
      W2          = covaW2,
      Delta       = Delta,
      Xi          = Xi,
      etime       = etime,
      retry_count = retry_count,
      SE          = SE,
      likvalue    = likvalue
    ))
  } else {
    return(list(
      params      = NA,
      Z           = time,
      X1          = covaX[, 1],
      X2          = covaX[, 2],
      W1          = covaW[, 1],
      W2          = covaW[, 2],
      Delta       = Delta,
      Xi          = Xi,
      etime       = etime,
      retry_count = retry_count,
      SE          = NA,
      likvalue    = NA
    ))
  }
}
