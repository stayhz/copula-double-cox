#----------------------------------------------------------------------
# Fit the copula-single-Cox model using data Ddata.
# Supports Weibull and Generalized Exponential (GenExp) marginal distributions.
#----------------------------------------------------------------------

method_singlecox <- function(n, p, q, copula_name, dist, numseed0, Ddata) {
  
  ## === Load required libraries and source files === ##
  library(copula)
  library(ucminf)
  source("func_marg.R")
  source("func_partialcop.R")
  source("likelihood.R")
  source("copulacdf2.R")
  source("func_revkendall.R")
  source("easydoublecox.R")
  
  ### Step 1: Copula parameter transformation based on target Kendall's tau ###
  cop_para0 <- func_revkendall(copula_name[2], 0.5)
  
  if (copula_name[2] == "frank") {
    # no transformation needed
  } else if (copula_name[2] == "gumbel") {
    cop_para0 <- log(cop_para0 - 1)
  } else if (copula_name[2] == "gaussian") {
    cop_para0 <- atanh(cop_para0)
  } else if (copula_name[2] == "clayton") {
    cop_para0 <- log(cop_para0)
  }
  
  ### Step 2: Initial parameter values for simplified  likelihood fitting ###
  if (dist == "Weibull") {
    initial_params <- c(3, 0, 3, 0, cop_para0)
  } else if (dist == "GenExp") {
    initial_params <- c(-3, 0, -3, 0, cop_para0)
  } else {
    stop("Unsupported distribution type.")
  }
  
  # Define simplified objective functions (no covariate effect)
  objective_easy <- function(y) {
    # cat("Current raw params:", y, "\n")
    likelihoodeasy(y, Ddata, p, q, dist, copula_name)
  }
  
  ### Step 3: Fit simplified 'easy' model ###
  result_easy <- ucminf(par = initial_params, fn = objective_easy, gr = NULL, hessian = 1)
  param_easy <- result_easy$par
  
  ### Step 4: Use estimated parameters as baseline for full model ###
  initial_full <- c(
    param_easy[1:2], rep(0, p),
    param_easy[3:4], rep(0, p),
    cop_para0
  )
  
  # Full model log-likelihood
  objective_full <- function(y) {
    # cat("Current raw params:", y, "\n")
    likelihood_singlecox(y, Ddata, p, q, dist, copula_name)
  }
  
  ### Step 5: Fit full likelihood model ###
  result_full <- ucminf(par = initial_full, fn = objective_full, gr = NULL, hessian = 1)
  est_params <- result_full$par
  hessian_matrix <- result_full$hessian
  
  ### Step 6: Compute standard errors from Hessian ###
  covariance_matrix <- ginv(hessian_matrix)
  SE <- sqrt(diag(covariance_matrix))
  
  loglik_value <- objective_full(est_params)
  
  
  #### Output results ###
  return(list(
    params = est_params,
    SE = SE,
    value = loglik_value
  ))
}
