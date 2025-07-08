####################################################################
###################################################################
# Log-likelihood function under the copula-double-Cox model.

likelihood <- function(y, D, p, q, dist, copula_name) {
  allpara <- y
  
  # Extract observed data
  Z <- D$time
  Delta <- D$Delta
  Xi <- D$Xi
  
  # Extract covariates and reshape into matrices
  X <- matrix(unlist(D[, 3:(p + 2)]), ncol = p)
  W <- matrix(unlist(D[, (p + 3):(p + q + 2)]), ncol = q)
  n <- length(Z)
  
  # Extract model parameters
  paraT <- allpara[1:(p * 2 + 2)]
  paraC <- allpara[(p * 2 + 3):((p + q) * 2 + 4)]
  cop_para0 <- allpara[length(allpara)]
  
  # Transform copula parameter to valid domain
  copula_type <- copula_name[2]
  cop_para <- switch(copula_type,
                     "frank" = cop_para0,
                     "gumbel" = exp(cop_para0) + 1,
                     "gaussian" = tanh(cop_para0),
                     "clayton" = exp(cop_para0),
                     stop("Unsupported copula type")
  )
  
  # Define copula CDF and partial derivatives
  cop <- function(u, v) copulacdf2(copula_type, cbind(u, v), cop_para)
  partials <- func_partialcop(copula_type)
  dcop1 <- partials$dcop1
  dcop2 <- partials$dcop2
  
  # Initialize log-likelihood contributions
  Lvec <- numeric(n)
  Lvec[1] <- 0  # initial dummy value (if needed)
  
  for (i in 1:n) {
    z <- Z[i]
    
    margT <- func_marg(z, X[i, ], p, dist, paraT)
    margC <- func_marg(z, W[i, ], q, dist, paraC)
    fT <- margT$densf
    fC <- margC$densf
    u <- margT$distF
    v <- margC$distF
    
    # Compute likelihood contributions based on observation type
    if (Delta[i] == 1 && Xi[i] == 0) {
      Lvec[i] <- log(fT) + log(abs(1 - dcop1(u, v, cop_para)))
    } else if (Delta[i] == 0 && Xi[i] == 1) {
      Lvec[i] <- log(fC) + log(abs(1 - dcop2(u, v, cop_para)))
    } else if (Delta[i] == 0 && Xi[i] == 0) {
      Lvec[i] <- log(max(1 - u - v + cop(u, v), 1e-20))
    }
  }
  
  return(-sum(Lvec))  # negative log-likelihood
}


####################################################################
###################################################################
# Log-likelihood function under a simplified "copula-single-Cox"  model.
# Assumes both the shape parameters are covariate-independent.



likelihood_singlecox <- function(y, D, p, q, dist, copula_name) {
  allpara <- y
  
  # Extract observed data
  Z <- D$time
  Z0 <- D$trunc
  Delta <- D$Delta
  Xi <- D$Xi
  X <- matrix(unlist(D[, 3:(p + 2)]), ncol = p)
  W <- matrix(unlist(D[, (p + 3):(p + q + 2)]), ncol = q)
  n <- length(Z)
  
  # Extract model parameters for T and C
  paraT <- allpara[1:(p + 2)]
  paraC <- allpara[(p + 3):(p + q + 4)]
  cop_para0 <- allpara[length(allpara)]
  
  # Transform copula parameter
  copula_type <- copula_name[2]
  cop_para <- switch(copula_type,
                     "frank" = cop_para0,
                     "gumbel" = exp(cop_para0) + 1,
                     "gaussian" = tanh(cop_para0),
                     "clayton" = exp(cop_para0),
                     stop("Unsupported copula type")
  )
  
  # Copula functions
  cop <- function(u, v) copulacdf2(copula_type, cbind(u, v), cop_para)
  partials <- func_partialcop(copula_type)
  dcop1 <- partials$dcop1
  dcop2 <- partials$dcop2
  
  # Compute log-likelihood
  Lvec <- numeric(n)
  Lvec[1] <- 0
  
  for (i in 1:n) {
    z <- Z[i]
    
    margT <- func_margsingle(z, X[i, ], p, dist, paraT)
    margC <- func_margsingle(z, W[i, ], q, dist, paraC)
    fT <- margT$densf
    fC <- margC$densf
    u <- margT$distF
    v <- margC$distF
    
    if (Delta[i] == 1 && Xi[i] == 0) {
      Lvec[i] <- log(fT) + log(abs(1 - dcop1(u, v, cop_para)))
    } else if (Delta[i] == 0 && Xi[i] == 1) {
      Lvec[i] <- log(fC) + log(abs(1 - dcop2(u, v, cop_para)))
    } else if (Delta[i] == 0 && Xi[i] == 0) {
      Lvec[i] <- log(max(1 - u - v + cop(u, v), 1e-20))
    }
  }
  
  return(-sum(Lvec))
}


####################################################################
###################################################################
# Log-likelihood for the independent double-Cox model
# It supports Weibull and Generalized Exponential distributions.

likelihood_indep <- function(y, D, p, q, dist, copula_name) {
  # --- Extract parameter vector ---
  allpara <- y
  
  #  Extract observed data 
  Z     <- D$time     # Event or censoring time
  Z0    <- D$trunc    # Truncation (entry) time
  Delta <- D$Delta    # Event indicator for T
  Xi    <- D$Xi       # Event indicator for C
  
  #  Covariates for T and C
  X <- cbind(D$covaX1, D$covaX2)
  W <- cbind(D$covaW1, D$covaW2)
  n <- length(Z)
  
  # Parse distribution-specific parameters 
  if (dist == "Weibull") {
    paraT <- allpara[1:(p * 2 + 2)]
    paraC <- allpara[(p * 2 + 3):((p + q) * 2 + 4)]
  } else if (dist == "GenGompertz") {
    paraT <- allpara[1:(p * 2 + 3)]
    paraC <- allpara[(p * 2 + 4):((p + q) * 2 + 6)]
  } else if (dist == "GenExp") {
    paraT <- allpara[1:(p * 2 + 2)]
    paraC <- allpara[(p * 2 + 3):((p + q) * 2 + 4)]
  } else {
    stop("Unsupported distribution type.")
  }
  
  # Initialize log-likelihood components
  Lvec   <- numeric(n)
  check1 <- numeric(n)
  check2 <- numeric(n)
  check3 <- numeric(n)
  Lvec[1] <- 0  # The first sample is often used for baseline initialization
  
  # Loop over each observation to compute contribution to log-likelihood
  for (i in 1:n) {
    z1 <- Z[i]
    
    # Compute marginal density and distribution for T
    margT <- func_marg(z1, X[i, ], p, dist, paraT)
    fT    <- margT$densf
    u     <- margT$distF
    
    # Compute marginal density and distribution for C
    margC <- func_marg(z1, W[i, ], q, dist, paraC)
    fC    <- margC$densf
    v     <- margC$distF
    
    # Initialize log-likelihood contribution
    L <- 0
    
    # Case 1: Observed failure time for T, censored for C
    if (Delta[i] == 1 && Xi[i] == 0) {
      L1         <- log(fT) + log(1 - v)
      L          <- L1
      check1[i]  <- L1
      
      # Case 2: Observed failure time for C, censored for T
    } else if (Delta[i] == 0 && Xi[i] == 1) {
      L2         <- log(fC) + log(1 - u)
      L          <- L2
      check2[i]  <- L2
      
      # Case 3: Both are censored (non-informative)
    } else if (Delta[i] == 0 && Xi[i] == 0) {
      L3         <- log(1 - u - v + u * v)  # independence assumption: P(T>z,C>z) = (1-F_T)(1-F_C)
      L          <- L3
      check3[i]  <- L3
    }
    
    Lvec[i] <- L
  }
  
  # --- Return negative log-likelihood ---
  return(-sum(Lvec))
}
