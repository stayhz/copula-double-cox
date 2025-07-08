
# Generate data from copula-based log-normal model

gendata_lognorm <- function(numseed, num1, n, p, q, truepara, copula_name, dist, pA) {
  
  set.seed(numseed)  # Set random seed
  
  # Extract copula type
  copula_name1 <- copula_name[1]
  
  # Extract parameter vectors
  paraT     <- truepara[1:(p * 2 + 2)]
  paraC     <- truepara[(p * 2 + 3):((p + q) * 2 + 4)]
  cop_para0 <- truepara[length(truepara)]
  
  a1 <- paraT[1]
  a2 <- paraT[2]
  betascale <- paraT[3:(p + 2)]
  betashape <- paraT[(p + 3):(p * 2 + 2)]
  
  b1 <- paraC[1]
  b2 <- paraC[2]
  etascale <- paraC[3:(q + 2)]
  etashape <- paraC[(q + 3):(q * 2 + 2)]
  
  # Generate covariates
  X2 <- rbinom(n, size = 1, prob = 0.5)
  X1 <- rnorm(n, mean = 0, sd = sqrt(0.2))
  X0 <- cbind(X1, X2)
  
  W2 <- runif(n, min = -1, max = 1)
  W0 <- cbind(X1, W2)
  
  # Construct copula object
  copulaspe <- switch(copula_name1,
                      "fgm"      = fgmCopula(param = cop_para0, dim = 2),
                      "clayton"  = claytonCopula(param = cop_para0, dim = 2),
                      "gumbel"   = gumbelCopula(param = cop_para0, dim = 2),
                      "frank"    = frankCopula(param = cop_para0, dim = 2),
                      "gaussian" = normalCopula(param = cop_para0, dim = 2),
                      stop("Unsupported copula type"))
  
  # Generate dependent uniform variables
  U <- rCopula(n, copulaspe)
  U1 <- U[, 1]
  U2 <- U[, 2]
  
  # Calculate scale and shape for log-normal model
  coxscaleT <- exp(X0 %*% betascale)
  coxscaleC <- exp(W0 %*% etascale)
  
  # Fixed log-normal parameters (can be changed for sensitivity)
  muT <- 0.8; sigmaT <- 0.7
  muC <- 1.5; sigmaC <- 1.0
  
  # Generate log-normal failure and censoring times
  T <- exp(muT + sigmaT * qnorm(1 - (1 - U1)^(1 / coxscaleT)))
  C <- exp(muC + sigmaC * qnorm(1 - (1 - U2)^(1 / coxscaleC)))
  
  # Administrative censoring
  A <- runif(n, 0, pA)
  Z <- pmin(T, C, A)
  Delta0 <- as.numeric(Z == T)
  Xi0    <- as.numeric(Z == C)
  
  # Combine and sort by observed time
  obs    <- cbind(Z, X0, W0, Delta0, Xi0)
  obsnew <- obs[order(obs[, 1]), ]
  
  # Extract variables
  time  <- obsnew[, 1]
  covaX <- obsnew[, 2:(1 + p)]
  covaW <- obsnew[, (2 + p):(p + q + 1)]
  Delta <- obsnew[, (p + q + 2)]
  Xi    <- obsnew[, (p + q + 3)]
  trunc <- rep(0, n)
  
  # Return generated data
  return(list(
    time  = time,
    covaX = covaX,
    covaW = covaW,
    Delta = Delta,
    Xi    = Xi,
    trunc = trunc
  ))
}
