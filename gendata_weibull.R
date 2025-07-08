
####################################################################################
####################################################################################

# Generate dependent censoring data under a Weibull distribution.

gendata_weibull <- function(numseed, num1, n, p, q, truepara, copula_name, dist, pA) {
  
  set.seed(numseed)
  
  copula_name1 <- copula_name[1]
  
  # Extract parameters
  paraT <- truepara[1:(p * 2 + 2)]
  paraC <- truepara[(p * 2 + 3):((p + q) * 2 + 4)]
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
  
  # Construct copula
  copulaspe <- switch(copula_name1,
                      fgm = fgmCopula(param = cop_para0, dim = 2),
                      clayton = claytonCopula(param = cop_para0, dim = 2),
                      gumbel = gumbelCopula(param = cop_para0, dim = 2),
                      frank = frankCopula(param = cop_para0, dim = 2),
                      gaussian = normalCopula(param = cop_para0, dim = 2),
                      stop("Unsupported copula type"))
  
  # Generate joint uniform samples
  U <- rCopula(n, copulaspe)
  U1 <- U[, 1]
  U2 <- U[, 2]
  
  HT <- -log(1 - U1)
  HC <- -log(1 - U2)
  
  # Compute  scale and shape
  coxscaleT <- exp(X0 %*% betascale)
  coxshapeT <- exp(X0 %*% betashape)
  coxscaleC <- exp(W0 %*% etascale)
  coxshapeC <- exp(W0 %*% etashape)
  
  # Weibull
  T <- a1 * exp((log(HT) - log(coxscaleT)) / (a2 * coxshapeT))
  C <- b1 * exp((log(HC) - log(coxscaleC)) / (b2 * coxshapeC))
  
  # Apply administrative censoring
  A <- runif(n, 0, pA)
  Z <- pmin(T, C, A)
  Delta0 <- as.numeric(Z == T)
  Xi0 <- as.numeric(Z == C)
  
  # Sort data
  obs <- cbind(Z, X0, W0, Delta0, Xi0)
  obsnew <- obs[order(obs[, 1]), ]
  
  time <- obsnew[, 1]
  covaX <- obsnew[, 2:(1 + p)]
  covaW <- obsnew[, (2 + p):(p + q + 1)]
  Delta <- obsnew[, (p + q + 2)]
  Xi <- obsnew[, (p + q + 3)]
  trunc <- rep(0, n)
  
  return(list(
    time = time,
    covaX = covaX,
    covaW = covaW,
    Delta = Delta,
    Xi = Xi,
    trunc = trunc
  ))
}


####################################################################################
####################################################################################
# Generate independent censoring data under a Weibull distribution.

gendata_weibull_indep <- function(numseed, num1, n, p, q, truepara, copula_name, dist, pA) {
  
  set.seed(numseed)
  
  # Extract parameters
  paraT <- truepara[1:(p * 2 + 2)]
  paraC <- truepara[(p * 2 + 3):((p + q) * 2 + 4)]
  
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
  
  # Generate independent uniforms
  U1 <- runif(n)
  U2 <- runif(n)
  
  HT <- -log(1 - U1)
  HC <- -log(1 - U2)
  
  # Compute scale and shape
  coxscaleT <- exp(X0 %*% betascale)
  coxshapeT <- exp(X0 %*% betashape)
  coxscaleC <- exp(W0 %*% etascale)
  coxshapeC <- exp(W0 %*% etashape)
  
  # Weibull 
  T <- a1 * exp((log(HT) - log(coxscaleT)) / (a2 * coxshapeT))
  C <- b1 * exp((log(HC) - log(coxscaleC)) / (b2 * coxshapeC))
  
  # Apply administrative censoring
  A <- runif(n, 0, pA)
  Z <- pmin(T, C, A)
  Delta0 <- as.numeric(Z == T)
  Xi0 <- as.numeric(Z == C)
  
  # Sort data
  obs <- cbind(Z, X0, W0, Delta0, Xi0)
  obsnew <- obs[order(obs[, 1]), ]
  
  time <- obsnew[, 1]
  covaX <- obsnew[, 2:(1 + p)]
  covaW <- obsnew[, (2 + p):(p + q + 1)]
  Delta <- obsnew[, (p + q + 2)]
  Xi <- obsnew[, (p + q + 3)]
  trunc <- rep(0, n)
  
  return(list(
    time = time,
    covaX = covaX,
    covaW = covaW,
    Delta = Delta,
    Xi = Xi,
    trunc = trunc
  ))
}



