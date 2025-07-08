
####################################################################
###################################################################
# Compute marginal distribution for dependent censoring data under the simplified double-Cox model.
# Assumes both margins follow Weibull or GenExp distributions without covariate effects on parameters.

func_margeasy <- function(Z, X, dist, para) {
  if (dist == "Weibull") {
    a1 <- exp(para[1])  # scale
    b1 <- exp(para[2])  # shape
    
    cumh <- (Z / a1)^b1
    survS <- exp(-cumh)
    distF <- 1 - survS
    densf <- b1 * (Z^(b1 - 1)) / (a1^b1) * survS
    mufull <- densf / survS
    
  } else if (dist == "GenExp") {
    a1 <- exp(para[1])
    b1 <- exp(para[2])
    
    termexp <- exp(-a1 * Z)
    distF <- (1 - termexp)^b1
    survS <- 1 - distF
    cumh <- -log(survS)
    densf <- a1 * b1 * termexp * (1 - termexp)^(b1 - 1)
    mufull <- densf / survS
  }
  
  return(list(distF = distF, cumh = cumh, survS = survS, densf = densf, mufull = mufull))
}



####################################################################
###################################################################
# Log-likelihood for dependent censoring data under the simplified double-Cox model.

likelihoodeasy <- function(y, D, p, q, dist, copula_name) {
  allpara <- y
  
  Z <- D$time
  Delta <- D$Delta
  Xi <- D$Xi
  X <- as.matrix(D[, 3:(p + 2)])
  W <- as.matrix(D[, (p + 3):(p + q + 2)])
  n <- length(Z)
  
  paraT <- allpara[1:2]
  paraC <- allpara[3:4]
  cop_para0 <- allpara[length(allpara)]
  
  # Copula parameter transformation
  copula_type <- copula_name[2]
  cop_para <- switch(copula_type,
                     frank = cop_para0,
                     gumbel = exp(cop_para0) + 1,
                     gaussian = tanh(cop_para0),
                     clayton = exp(cop_para0),
                     stop("Unsupported copula")
  )
  
  # Load copula CDF and partial derivative functions
  cop <- function(u, v) copulacdf2(copula_type, cbind(u, v), cop_para)
  partials <- func_partialcop(copula_type)
  dcop1 <- partials$dcop1
  dcop2 <- partials$dcop2
  
  logL <- numeric(n)
  
  for (i in 1:n) {
    z <- Z[i]
    margT <- func_margeasy(z, X[i, ], dist, paraT)
    margC <- func_margeasy(z, W[i, ], dist, paraC)
    fT <- margT$densf
    fC <- margC$densf
    u <- margT$distF
    v <- margC$distF
    
    logLi <- NA
    
    if (Delta[i] == 1 && Xi[i] == 0) {
      logLi <- log(fT) + log(abs(1 - dcop1(u, v, cop_para)))
    } else if (Delta[i] == 0 && Xi[i] == 1) {
      logLi <- log(fC) + log(abs(1 - dcop2(u, v, cop_para)))
    } else if (Delta[i] == 0 && Xi[i] == 0) {
      logLi <- log(max(1 - u - v + cop(u, v), 1e-20))
    }
    
    logL[i] <- logLi #log-likelihood
  }
  
  return(-sum(logL))
}


####################################################################
###################################################################

# Log-likelihood under the assumption of independent censoring.
# Margins still follow Weibull or GenExp distributions.

likelihoodeasy_indep <- function(y, D, p, q, dist, copula_name) {
  allpara <- y
  
  Z <- D$time
  Delta <- D$Delta
  Xi <- D$Xi
  X <- cbind(D$covaX1, D$covaX2)
  W <- cbind(D$covaW1, D$covaW2)
  n <- length(Z)
  
  paraT <- allpara[1:2]
  paraC <- allpara[3:4]
  
  logL <- numeric(n)
  
  for (i in 1:n) {
    z <- Z[i]
    margT <- func_margeasy(z, X[i, ], dist, paraT)
    margC <- func_margeasy(z, W[i, ], dist, paraC)
    fT <- margT$densf
    fC <- margC$densf
    u <- margT$distF
    v <- margC$distF
    
    logLi <- NA
    
    if (Delta[i] == 1 && Xi[i] == 0) {
      logLi <- log(fT) + log(1 - v)
    } else if (Delta[i] == 0 && Xi[i] == 1) {
      logLi <- log(fC) + log(1 - u)
    } else if (Delta[i] == 0 && Xi[i] == 0) {
      logLi <- log(max(1 - u - v + u * v, 1e-20))
    }
    
    logL[i] <- logLi
  }
  
  return(-sum(logL))
}


