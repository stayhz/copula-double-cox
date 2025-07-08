#---------------------------------------------------------------
# Simulation under independent censoring data,
# then fit data using the Cox model for T and C.
#---------------------------------------------------------------

simu_indepsingle <- function(vv, n, p, q, truepara, copula_name, dist, numseed0, pA) {
  start_time <- Sys.time()
  
  # --- Load required libraries ---
  library(copula)
  library(ucminf)
  library(parallel)
  library(MASS)       # For multivariate normal distributions
  library(survival)   # For Cox model
  
  # --- Load external function definitions ---
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
  
  retry_count <- 0
  success <- FALSE
  
  while (!success && retry_count <= 5) {
    try({
      numseed <- numseed0 + vv + retry_count
      
      # --- Generate data based on selected distribution ---
      gendata <- switch(
        dist,
        Weibull     = gendata_weibull(numseed, num1, n, p, q, truepara, copula_name, dist, pA),
        GenExp      = gendata_genexp(numseed, num1, n, p, q, truepara, copula_name, dist, pA),
        stop("Unsupported distribution type.")
      )
      
      covaX  <- gendata$covaX
      covaW  <- gendata$covaW
      trunc  <- gendata$trunc
      time   <- gendata$time
      Delta  <- gendata$Delta
      Xi     <- gendata$Xi
      
      # --- Assemble data frame ---
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
      
      # --- Check censoring proportions ---
      prop_Delta      <- mean(Delta == 1)
      prop_Xi         <- mean(Xi == 1)
      prop_both_zero  <- mean((Delta + Xi) == 0)
      
      # =======================================================
      # Fit Cox model for T ~ X
      # =======================================================
      cox_model1 <- coxph(
        formula = Surv(time, Delta) ~ covaX1 + covaX2,
        data = Ddata,
        method = "efron"
      )
      
      estpara1 <- coef(cox_model1)
      SE1      <- sqrt(diag(cox_model1$var))
      breslow1 <- basehaz(cox_model1, centered = FALSE)
      
      # =======================================================
      # Fit Cox model for C ~ W
      # =======================================================
      cox_model2 <- coxph(
        formula = Surv(time, Xi) ~ covaW1 + covaW2,
        data = Ddata,
        method = "efron"
      )
      
      estpara2 <- coef(cox_model2)
      SE2      <- sqrt(diag(cox_model2$var))
      breslow2 <- basehaz(cox_model2, centered = FALSE)
      
      # --- Combine parameters and SEs into full-length vectors ---
      estpara <- c(0, 0, estpara1, 0, 0, 0, 0, estpara2, 0, 0, 0)
      SE      <- c(0, 0, SE1,      0, 0, 0, 0, SE2,      0, 0, 0)
      
      success <- TRUE
    }, silent = TRUE)
    
    retry_count <- retry_count + 1
  }
  
  etime <- Sys.time() - start_time
  
  # --- Return results ---
  if (success) {
    print(vv)  # progress indicator
    return(list(
      params      = estpara,
      Z           = time,
      X1          = covaX[, 1],
      X2          = covaX[, 2],
      Delta       = Delta,
      Xi          = Xi,
      etime       = etime,
      retry_count = retry_count,
      SE          = SE,
      breslow1    = breslow1,
      cox_model1  = cox_model1,
      breslow2    = breslow2,
      cox_model2  = cox_model2
    ))
  } else {
    warning(sprintf("Simulation %d failed after %d retries.", vv, retry_count))
    return(list(
      params      = NA,
      Z           = time,
      X1          = covaX[, 1],
      X2          = covaX[, 2],
      Delta       = Delta,
      Xi          = Xi,
      etime       = etime,
      retry_count = retry_count,
      SE          = NA,
      breslow1    = NA,
      cox_model1  = NA,
      breslow2    = NA,
      cox_model2  = NA
    ))
  }
}
