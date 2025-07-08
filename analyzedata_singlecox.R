# Analyze estimation results under the single-Cox model
# Return Bias, SEE, ASE and CP.


analyzedata_singlecox <- function(results, truepara, copula_name) {
  d <- length(truepara)
  copula_name1 <- copula_name[1]
  copula_name2 <- copula_name[2]
  
  # Initialize storage
  matresttheta <- matrix(0, num1, d - 4)
  matrZ <- matrX1 <- matrX2 <- matrDelta <- matrXi <- matrretry <- matretime <- matrix(0, num1, n)
  matrse <- matrix(0, num1, d - 4)
  CPvec <- matrix(0, 1, d)
  tauvec <- matrse_tau <- numeric(num1)
  
  for (i in 1:num1) {
    res <- results[[i]]
    matresttheta[i, ] <- res$params
    matrZ[i, ] <- res$Z
    matrX1[i, ] <- res$X1
    matrX2[i, ] <- res$X2
    matrDelta[i, ] <- res$Delta
    matrXi[i, ] <- res$Xi
    matrretry[i, ] <- res$retry_count
    matretime[i, ] <- res$etime
    matrse[i, ] <- res$SE
  }
  
  # Remove incomplete rows
  complete_rows <- complete.cases(matresttheta) & complete.cases(matrse)
  matresttheta <- matresttheta[complete_rows, ]
  matrse <- matrse[complete_rows, ]
  num1 <- nrow(matresttheta)
  
  # Expand to full d columns with inserted zeros (for removed transformation params)
  insert_zero <- function(mat, idx) {
    out <- matrix(0, nrow = nrow(mat), ncol = d)
    out[, idx] <- mat
    return(out)
  }
  keepidx <- c(1:4, 5:8, 9)
  insertidx <- c(1:4, 7:10, 13)
  matresttheta <- insert_zero(matresttheta, insertidx)
  matrse <- insert_zero(matrse, insertidx)
  
  # Transform copula parameter to original scale
  transform_copula_param <- function(x, name) {
    switch(name,
           "frank" = x,
           "gumbel" = exp(x) + 1,
           "gaussian" = tanh(x),
           "clayton" = exp(x),
           x)
  }
  matresttheta2 <- matresttheta
  matresttheta2[, d] <- transform_copula_param(matresttheta[, d], copula_name2)
  
  # Determine selvec and truestar based on dist
  selvec <- c(1, 2, 7, 8)
  truestar <- c(log(truepara[1:2]), truepara[3:6], log(truepara[7:8]), truepara[9:13])
  matresttheta2[, selvec] <- exp(matresttheta[, selvec])
  
  # Compute tau and SE
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
  for (j in 1:num1) {
    tau_se <- compute_tau_and_se(copula_name2, matresttheta[j, d], matrse[j, d])
    tauvec[j] <- tau_se$tau
    matrse_tau[j] <- tau_se$SE
  }
  
  # True tau
  truetau <- ifelse(truepara[d] == 0, 0, copulastat(copula_name1, truepara[d]))
  
  # Bias and proportion bias
  meanestpara <- colMeans(matresttheta)
  bias <- meanestpara - truestar
  bias_exp <- colMeans(matresttheta2) - truepara
  propbias <- bias / truestar
  propbias_exp <- bias_exp / truepara
  bias_tau <- mean(tauvec) - truetau
  propbias_tau <- bias_tau / truetau
  bias_exp[d] <- bias_tau
  
  # CP calculation
  for (i in 1:d) {
    CPvec[i] <- mean(sapply(1:num1, function(j) {
      if (i %in% selvec) {
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
