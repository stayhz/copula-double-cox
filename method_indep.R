#----------------------------------------------------------------------
# Fit the independent-double-Cox model using data Ddata.
# Supports Weibull and Generalized Exponential (GenExp) marginal distributions.
#----------------------------------------------------------------------

method_indep <- function(n, p, q, copula_name, dist, numseed0, Ddata) {
  
  ## === Load required libraries and source files === ##
  library(copula)
  library(ucminf)
  source("func_marg.R")
  source("func_partialcop.R")
  source("likelihood.R")
  source("copulacdf2.R")
  source("func_revkendall.R")
  source("easydoublecox.R")
  
  ## === Step 1: Initialize baseline parameters === ##
  if (dist == "Weibull") {
    initial_params <- c(rep(0, 4), 0)  
  } else if (dist == "GenExp") {
    initial_params <- c(rep(0, 4), 0)
  } else {
    stop("Unsupported distribution type.")
  }
  
  
  ## === Step 2: Define simplified objective functions (no covariate effect) === ##
  objective_easy <- function(y) {
    # cat("Current raw params:", y, "\n")
    likelihoodeasy_indep(y, Ddata, p, q, dist, copula_name)
  }
  
  
  ### Step 3: Fit simplified 'easy' model ###
  result_easy <- ucminf(par = initial_params, fn = objective_easy, gr = NULL, hessian = 1)
  param_easy <- result_easy$par
  
  
  ### Step 4: Use estimated parameters as baseline for full model ###
  initial_full <- c(
    param_easy[1:2], rep(0, p + q),
    param_easy[3:4], rep(0, p + q),
    0
  )
  
  
  ## === Step 5: Define full log-likelihood function === ##
  # Full model log-likelihood
  objective_full <- function(y) {
    # cat("Current raw params:", y, "\n")
    likelihood_indep(y, Ddata, p, q, dist, copula_name)
  }
  
  ### Step 5: Fit full likelihood model ###
  result_full <- ucminf(par = initial_full, fn = objective_full, gr = NULL, hessian = 1)
  est_params <- result_full$par
  hessian_matrix <- result_full$hessian
  
  ### Step 6: Compute standard errors from Hessian ###
  covariance_matrix <- ginv(hessian_matrix)
  SE <- sqrt(diag(covariance_matrix))
  
  loglik_value <- objective_full(est_params)
  
  ### Output results ###
  return(list(
    params = est_params,
    SE = SE,
    value = loglik_value
  ))
}
