#----------------------------------------------------------------------
# Main workflow: data generation, initial estimation, and full model fitting
# for the copula-double-Cox model under dependent censoring.
# Supports Weibull and Generalized Exponential (GenExp) marginal distributions.
#----------------------------------------------------------------------

method_doublecox <- function(n, p, q, copula_name, dist, numseed0, Ddata) {
  # Load required packages and source auxiliary functions
  library(copula)
  library(ucminf)
  library(MASS)  # for ginv()
  
  source("func_marg.R")
  source("func_partialcop.R")
  source("likelihood.R")
  source("copulacdf2.R")
  source("func_revkendall.R")
  source("easydoublecox.R")
  
  ### Step 1: Copula parameter transformation based on target Kendall's tau ###
  cop_para0 <- func_revkendall(copula_name[2], 0.5)
  cop_para0 <- switch(copula_name[2],
                      "frank"    = cop_para0,
                      "gumbel"   = log(cop_para0 - 1),
                      "gaussian" = atanh(cop_para0),
                      "clayton"  = log(cop_para0),
                      stop("Unsupported copula type")
  )
  
  ### Step 2: Initial parameter values for simplified  likelihood fitting ###
  initial_params <- switch(dist,
                           "Weibull" = c(3, 0, 3, 0, cop_para0),
                           "GenExp"  = c(-3, 0, -3, 0, cop_para0),
                           stop("Unsupported distribution")
  )
  
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
    param_easy[1:2], rep(0, p + q),
    param_easy[3:4], rep(0, p + q),
    cop_para0
  )
  
  # Full model log-likelihood
  objective_full <- function(y) {
    # cat("Current raw params:", y, "\n")
    likelihood(y, Ddata, p, q, dist, copula_name)
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
