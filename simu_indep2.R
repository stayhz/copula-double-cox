#---------------------------------------------------------------
# Simulation under independent censoring data.
# Then fit data using the independent-double-Cox model.
#---------------------------------------------------------------

simu_indep2 <- function(vv, n, p, q, truepara, copula_name, dist, numseed0, pA) {
  start_time <- Sys.time()
  
  ## === Load required libraries === ##
  library(copula)
  library(ucminf)
  library(parallel)
  library(MASS)  # For ginv(), mvrnorm()
  
  ## === Load required external function files === ##
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
  
  ## === Retry mechanism in case of optimization failure === ##
  retry_count <- 0
  success <- FALSE
  
  while (!success && retry_count <= 5) {
    try({
      # --- Step 1: Set seed ---
      numseed <- numseed0 + vv + retry_count
      
      # --- Step 2: Simulate data under independence assumption ---
      gendata <- switch(
        dist,
        "Weibull"      = gendata_weibull_indep(numseed, num1, n, p, q, truepara, copula_name, dist, pA),
        "GenExp"       = gendata_genexp_indep(numseed, num1, n, p, q, truepara, copula_name, dist, pA),
        stop("Unsupported distribution type.")
      )
      
      covaX  <- gendata$covaX
      covaW  <- gendata$covaW
      trunc  <- gendata$trunc
      time   <- gendata$time
      Delta  <- gendata$Delta
      Xi     <- gendata$Xi
      
      # --- Step 3: Assemble data frame ---
      Ddata <- data.frame(
        trunc     = trunc,
        time      = time,
        covaX1    = covaX[, 1],
        covaX2    = covaX[, 2],
        covaW1    = covaW[, 1],
        covaW2    = covaW[, 2],
        Delta     = Delta,
        Xi        = Xi
      )
      
      
      ## --- (Optional) Check censoring proportions ---
      prop_Delta      <- sum(gendata$Delta == 1) / n
      prop_Xi         <- sum(gendata$Xi == 1) / n
      prop_both_zero  <- sum((gendata$Delta + gendata$Xi) == 0) / n
      
      # --- Step 4: Fit independent double-Cox model ---
      estresult <- method_indep(n, p, q, copula_name, dist, numseed0, Ddata)
      
      # --- Step 5: Extract estimation results ---
      estpara <- estresult$params
      SE      <- estresult$SE
      likvalue  <- -estresult$value  # log-Likelihood value 
      
      success <- TRUE
      
    }, silent = TRUE)
    
    retry_count <- retry_count + 1
  }
  
  etime <- Sys.time() - start_time
  
  # --- Step 6: Return results (successful or not) ---
  if (success) {
    return(list(
      params       = estpara,
      Z            = time,
      X1           = covaX[, 1],
      X2           = covaX[, 2],
      Delta        = Delta,
      Xi           = Xi,
      etime        = etime,
      retry_count  = retry_count,
      SE           = SE,
      likvalue     = likvalue
    ))
  } else {
    warning(sprintf("Simulation %d failed after %d retries.", vv, retry_count))
    return(list(
      params       = NA,
      Z            = time,
      X1           = covaX[, 1],
      X2           = covaX[, 2],
      Delta        = Delta,
      Xi           = Xi,
      etime        = etime,
      retry_count  = retry_count,
      SE           = NA,
      likvalue     = NA
    ))
  }
}
