

copulacdf2 <- function(copula_name, vec, copula_para) {
  r <- copula_para
  u <- vec[1]
  v <- vec[2]
  
  
  if (copula_name == "fgm") {
    cdfvalue <- u * v + copula_para * u * v * (1 - u) * (1 - v)
  } else if (copula_name == "frank") {
    cdfvalue <- pCopula(c(u, v), copula = frankCopula(param = copula_para))
  } else if (copula_name == "gaussian") {
    cdfvalue <- pCopula(c(u, v), copula = normalCopula(param = copula_para))
  } else if (copula_name == "clayton") {
    cdfvalue <- pCopula(c(u, v), copula = claytonCopula(param = copula_para))
  } else if (copula_name == "gumbel") {
    cdfvalue <- pCopula(c(u, v), copula = gumbelCopula(param = copula_para))
  } else {
    stop("Unsupported copula name.")
  }

  
  return(cdfvalue)
}






