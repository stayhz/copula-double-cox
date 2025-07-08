

copulastat <- function(cop_name, paramvec) {
 
  n=length(paramvec)
  tauvec <- rep(0,n)
  
  for (i in 1:n){
    param <- paramvec[i]
  # Kendall tau
  if (cop_name == "frank") {
    frank_copula <- frankCopula(param)
    kendall_tau <- tau(frank_copula)
  } else if (cop_name == "clayton") {
    clayton_copula <- claytonCopula(param)
    kendall_tau <- tau(clayton_copula)
  } else if (cop_name == "gumbel") {
    gumbel_copula <- gumbelCopula(param)
    kendall_tau <- tau(gumbel_copula)
  } else if (cop_name == "gaussian") {
    gaussian_copula <- normalCopula(param)
    kendall_tau <- tau(gaussian_copula)
    
  } else if (cop_name == "fgm") {
    kendall_tau <- 2 * param / 9
  } else {
    stop("Unknown copula name")
  }
  
    tauvec[i] <- kendall_tau
  }
  
  
  return(tauvec)
}
