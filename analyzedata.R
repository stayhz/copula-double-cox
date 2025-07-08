# Analyze estimation results under the dependence model
# Return Bias, SEE, ASE and CP.


analyzedata <- function(results, truepara, copula_name) {
  d <- length(truepara)
  copula_name1 <- copula_name[1]
  copula_name2 <- copula_name[2]
  
  # Initialize storage matrices
  matresttheta <- matrix(0, num1, d)
  matrse <- matrix(0, nrow = num1, ncol = d)
  matrZ <- matrX1 <- matrX2 <- matrDelta <- matrXi <- matrretry <- matretime  <- matrix(0, num1, n)
  matrlik <- matrix(0, num1, 1)
  tauvec <- matrix(0, num1, 1)
  CPvec <- matrix(0, 1, d)
  
  # Extract simulation results
  for (i in 1:num1) {
    res <- results[[i]]
    matresttheta[i, ] <- res$params
    matrZ[i, ]        <- res$Z
    matrX1[i, ]       <- res$X1
    matrX2[i, ]       <- res$X2
    matrDelta[i, ]    <- res$Delta
    matrXi[i, ]       <- res$Xi
    matrretry[i, ]    <- res$retry_count
    matretime[i, ]    <- res$etime
    matrse[i, ]       <- res$SE
    # matrlik[i, ]      <- res$likvalue
  }
  
  # Remove rows with NA
  complete_rows <- complete.cases(matresttheta) & complete.cases(matrse)
  matresttheta <- matresttheta[complete_rows, ]
  matrse <- matrse[complete_rows, ]
  num1 <- nrow(matresttheta)
  
  # Function for copula parameter transformation
  transform_copula_param <- function(x, name) {
    switch(name,
           "frank" = x,
           "gumbel" = exp(x) + 1,
           "gaussian" = tanh(x),
           "clayton" = exp(x),
           x)
  }
  
  # Function to compute tau and SE
  compute_tau_and_se <- function(copula_name, rstar, SEi) {
    if (copula_name == "frank") {
      r <- max(rstar, 1e-10)
      tau <- copulastat("frank", r)
      SE_tau <- SEi * abs(4 / r^2 - (tau - 1 + 4 / r) * 2 / r + 4 / (r * (exp(r) - 1)))
    } else if (copula_name == "gumbel") {
      r <- exp(rstar) + 1
      tau <- 1 - 1 / r
      SE_tau <- SEi * exp(rstar) / r^2
    } else if (copula_name == "gaussian") {
      r <- tanh(rstar)
      tau <- 2 * asin(r) / pi
      SE_tau <- (2 / pi) * SEi * (1 - r^2)
    } else if (copula_name == "clayton") {
      r <- exp(rstar)
      tau <- r / (r + 2)
      SE_tau <- SEi * 2 * r / (r + 2)^2
    } else {
      tau <- 0
      SE_tau <- 0
    }
    return(list(tau = tau, SE = SE_tau))
  }
  
  # Transform parameters
  exp_idx <- c(1, 2, 7, 8)
  matresttheta2 <- matresttheta
  matresttheta2[, exp_idx] <- exp(matresttheta[, exp_idx])
  matresttheta2[, d] <- transform_copula_param(matresttheta[, d], copula_name2)
  truestar <- c(log(truepara[1:2]), truepara[3:6], log(truepara[7:8]), truepara[9:13])
  
  # Compute tau and SE
  matrse_tau <- numeric(num1)
  for (j in 1:num1) {
    tau_se <- compute_tau_and_se(copula_name2, matresttheta[j, d], matrse[j, d])
    tauvec[j] <- tau_se$tau
    matrse_tau[j] <- tau_se$SE
  }
  
  truetau <- ifelse(truepara[d] == 0, 0, copulastat(copula_name1, truepara[d]))
  bias_tau <- mean(tauvec) - truetau
  propbias_tau <- bias_tau / truetau
  
  # Compute mean, bias, and proportional bias
  meanestpara <- colMeans(matresttheta)
  bias <- meanestpara - truestar
  bias_exp <- colMeans(matresttheta2) - truepara
  propbias <- bias / truestar
  propbias_exp <- bias_exp / truepara
  
  # Compute coverage probability
  for (i in 1:d) {
    count <- sum(sapply(1:num1, function(j) {
      if (i %in% exp_idx) {
        SEi <- matrse[j, i] * exp(matresttheta[j, i])
        CI <- exp(matresttheta[j, i]) + c(-1.96, 1.96) * SEi
        true_val <- truepara[i]
      } else if (i == d) {
        CI <- tauvec[j] + c(-1.96, 1.96) * matrse_tau[j]
        true_val <- truetau
      } else {
        SEi <- matrse[j, i]
        CI <- matresttheta[j, i] + c(-1.96, 1.96) * SEi
        true_val <- truepara[i]
      }
      true_val >= CI[1] && true_val <= CI[2]
    }))
    CPvec[i] <- count / num1
  }
  
  sdvec <- apply(matresttheta, 2, sd)
  sdvec[d] <- sd(tauvec)
  ASE <- colMeans(matrse)
  ASE[d] <- mean(matrse_tau)
  
  return(list(
    bias = bias,
    bias_exp = bias_exp,
    propbias = propbias,
    propbias_exp = propbias_exp,
    truestar = truestar,
    meanestpara = meanestpara,
    meanretry = mean(matrretry),
    meanse = colMeans(matrse),
    CPvec = t(CPvec),
    bias_tau = bias_tau,
    propbias_tau = propbias_tau,
    sdvec = sdvec,
    ASE = ASE
  ))
}
