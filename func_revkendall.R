
# Obtain the copula parameters from Kendall correlation
func_revkendall <- function(cop_name, tau) {
  
  # Kendall correlation for Frank copula
  fun <- function(x) {
    return(1 + 4 * integrate(function(t) t / (exp(t) - 1), 0, x)$value / (x^2) - 4 / x)
  }
  
  # Define the result variable to store the copula parameter
  if (cop_name == "frank") {
    # Solve the equation using uniroot function (equivalent to fzero in MATLAB)
    para <- uniroot(function(x) fun(x) - tau, c(0.0001, 100))$root
  } else if (cop_name == "clayton") {
    para <- (2 * tau) / (1 - tau)
  } else if (cop_name == "gumbel") {
    para <- 1 / (1 - tau)
  } else if (cop_name == "gaussian") {
    para <- 2 * sin(pi / 6) * sin(pi / 2 * tau)
  } else if (cop_name == "fgm") {
    para <- 9 * tau / 2
  } else {
    stop("Unknown copula name")
  }
  
  return(para)
}

