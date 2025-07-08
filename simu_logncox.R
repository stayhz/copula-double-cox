# Generate data from a copula-based log-normal model and fit Cox model

simu_logncox <- function(vv, n, p, q, truepara, copula_name, dist, numseed0, pA) {
  start_time <- Sys.time()
  
  # Load necessary libraries
  library(copula)
  library(ucminf)
  library(parallel)
  library(MASS)
  library(survival)
  
  # Load source functions
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
  source("gendata_lognorm.R")
  
  retry_count <- 0
  success <- FALSE
  
  while (!success && retry_count <= 5) {
    try({
      numseed <- numseed0 + vv + retry_count
      
      # Generate data from log-normal model
      gendata <- gendata_lognorm(numseed, num1, n, p, q, truepara, copula_name, dist, pA)
      
      covaX  <- gendata$covaX
      covaW  <- gendata$covaW
      trunc  <- gendata$trunc
      time   <- gendata$time
      Delta  <- gendata$Delta
      Xi     <- gendata$Xi
      
      # Extract covariates
      covaX1 <- covaX[, 1]
      covaX2 <- covaX[, 2]
      covaW1 <- covaW[, 1]
      covaW2 <- covaW[, 2]
      
      # Construct data frame
      Ddata <- data.frame(trunc, time, covaX1, covaX2, covaW1, covaW2, Delta, Xi)
      
      ###### Fit first Cox model for Delta ######
      cox_model1 <- coxph(Surv(time, Delta) ~ covaX1 + covaX2, data = Ddata, method = "efron")
      estpara1   <- coef(cox_model1)
      SE1        <- sqrt(diag(cox_model1$var))
      breslow1   <- basehaz(cox_model1, centered = FALSE)
      
      ###### Fit second Cox model for Xi ######
      cox_model2 <- coxph(Surv(time, Xi) ~ covaW1 + covaW2, data = Ddata, method = "efron")
      estpara2   <- coef(cox_model2)
      SE2        <- sqrt(diag(cox_model2$var))
      breslow2   <- basehaz(cox_model2, centered = FALSE)
      
      ###### Assemble full parameter vector and SE to match doublecox layout ######
      estpara <- c(0, 0, estpara1, 0, 0, 0, 0, estpara2, 0, 0, 0)
      SE      <- c(0, 0, SE1,      0, 0, 0, 0, SE2,      0, 0, 0)
      
      success <- TRUE
    }, silent = TRUE)
    retry_count <- retry_count + 1
  }
  
  etime <- Sys.time() - start_time
  
  if (success) {
    print(vv)
    return(list(
      params      = estpara,
      Z           = time,
      X1          = covaX1,
      X2          = covaX2,
      Delta       = Delta,
      Xi          = Xi,
      etime       = etime,
      retry_count = retry_count,
      SE          = SE,
      breslow1    = breslow1,
      breslow2    = breslow2,
      cox_model1  = cox_model1,
      cox_model2  = cox_model2
    ))
  } else {
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
      breslow2    = NA,
      cox_model1  = NA,
      cox_model2  = NA
    ))
  }
}
