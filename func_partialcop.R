#' @title func_partialcop   # 标题，描述函数的名称或作用，在帮助文档中显示为函数的简短标题。
#'
#' @description This function calculates the first and second partial derivatives 
#' of various copulas, including Frank, Clayton, and Gumbel. 
#' This function is not applicable for Gaussian copulas due to their complexity.
#' 
#' @param copula_name A character string indicating the name of the copula for which 
#' the derivatives are to be calculated. Acceptable values are "Frank", "Clayton", 
#' "FGM", and "Gumbel".
#' 
#' @return A list containing the first and second partial derivatives of the specified copula:
#' \item{dcop1}{The first partial derivative with respect to the first variable.}
#' \item{dcop2}{The first partial derivative with respect to the second variable.}
#' \item{dcop11}{The second partial derivative with respect to the first variable.}
#' \item{dcop22}{The second partial derivative with respect to the second variable.}
#'
#' @examples   # 示例，提供了函数使用的示例代码，帮助用户理解如何调用这个函数。
#' \dontrun{   # 表示该代码块不会在帮助文档中实际运行，但可以作为调用示例提供给用户。
#' derivatives <- func_partialcop("Frank")
#' }
#'
#' @export # 指示该函数应被导出，使其成为包的公共接口的一部分，从而用户可以使用该函数。
#'
#' 该函数用于计算给定 copula 的一阶和二阶偏导数，适用于 Frank、Clayton 和 Gumbel copulas。
#' 输入参数为 copula 名称，返回相应的偏导数值。


func_partialcop <- function(copula_name) {
  if (copula_name == "frank") {
    # Frank copula偏导
    dcop1 <- function(u, v, r) {
      -(1/r)*(1/ (1 + (exp(-r*u)-1)*(exp(-r*v)-1)/(exp(-r)-1))) * (-r)*exp(-r*u)*(exp(-r*v)-1)/(exp(-r)-1)
    }
    
    dcop2 <- function(u, v, r) {
      -(1/r)*(1/ (1 + (exp(-r*u)-1)*(exp(-r*v)-1)/(exp(-r)-1))) * (-r)*exp(-r*v)*(exp(-r*u)-1)/(exp(-r)-1)
    }
    
    dcop11 <- function(u, v, r) {
      (1/r)*(1/ (1 + (exp(-r*u)-1)*(exp(-r*v)-1)/(exp(-r)-1)))^2 * ((-r)*exp(-r*u)*(exp(-r*v)-1)/(exp(-r)-1))^2 - 
      (1/r)* (1/ (1 + (exp(-r*u)-1)*(exp(-r*v)-1)/(exp(-r)-1)))  * (r^2)*exp(-r*u)*(exp(-r*v)-1)/(exp(-r)-1)
    }
    
    dcop22 <- function(u, v, r) {
      (1/r)*(1/ (1 + (exp(-r*u)-1)*(exp(-r*v)-1)/(exp(-r)-1)))^2 * ((-r)*exp(-r*v)*(exp(-r*u)-1)/(exp(-r)-1))^2 - 
      (1/r)* (1/ (1 + (exp(-r*u)-1)*(exp(-r*v)-1)/(exp(-r)-1)))  * (r^2)*exp(-r*v)*(exp(-r*u)-1)/(exp(-r)-1)
    }
    
  } else if (copula_name == "clayton") {
    # Clayton copula偏导
    dcop1 <- function(x1, x2, r) {
      ( x1^(-r)+x2^(-r)-1 )^((-1/r) -1) *x1^(-r-1)
    }
    
    dcop2 <- function(x1, x2, r) {
      ( x1^(-r)+x2^(-r)-1 )^((-1/r) -1) *x2^(-r-1)
    }
    
    dcop11 <- function(x1, x2, rr) {
      (1+rr)*(x1^(-rr)+x2^(-rr)-1)^(-1/rr-2)*x1^(-2*rr-2) - (x1^(-rr)+x2^(-rr)-1)^(-1/rr-1)  * (rr+1)*x1^(-rr-2) 
    }
    
    dcop22 <- function(x1, x2, rr) {
      (1+rr)*(x1^(-rr)+x2^(-rr)-1)^(-1/rr-2)*x2^(-2*rr-2) - (x1^(-rr)+x2^(-rr)-1)^(-1/rr-1)  * (rr+1)*x2^(-rr-2) 
    }
    
  } else if (copula_name == "fgm") {
    # FGM copula偏导
    dcop1 <- function(x1, x2, r) {
      x2 + r*x2*(1-x1)*(1-x2) -  r*x1*x2*(1-x2)
    }
    
    dcop2 <- function(x1, x2, r) {
      x1 + r*x1*(1-x1)*(1-x2) -  r*x1*x2*(1-x1)
    }
    
    dcop11 <- function(x1, x2, rr) {
      -2*rr*x2*(1-x2) 
    }
    
    dcop22 <- function(x1, x2, rr) {
      -2*rr*x1*(1-x1) 
    }
    
  } else if (copula_name == "gumbel") {
    # Gumbel copula偏导
    dcop1 <- function(u, v, r) {
      -exp(-( ( (-log(u))^r + (-log(v))^r )  )^(1/r) ) * ( ( (-log(u))^r + (-log(v))^r )  )^(1/r -1) * (-log(u))^(r-1) *(-1/u)  
    }
    
    dcop2 <- function(u, v, r) {
      -exp(-( ( (-log(u))^r + (-log(v))^r )  )^(1/r) ) * ( ( (-log(u))^r + (-log(v))^r )  )^(1/r -1) * (-log(v))^(r-1) *(-1/v)
    }
    
    dcop11 <- function(u, v, r) {
      exp(-( ( (-log(u))^r + (-log(v))^r )  )^(1/r) ) * (( ( (-log(u))^r + (-log(v))^r )  )^(1/r -1) * (-log(u))^(r-1) *(-1/u))^2 -
      exp(-( ( (-log(u))^r + (-log(v))^r )  )^(1/r) ) * (1/r -1) *( ( (-log(u))^r + (-log(v))^r )  )^(1/r -2) * r * ((-log(u))^(r-1) *(-1/u))^2 - 
      exp(-( ( (-log(u))^r + (-log(v))^r )  )^(1/r) ) * ( ( (-log(u))^r + (-log(v))^r )  )^(1/r -1) * (r-1) * (-log(u))^(r-2) *(-1/u)^2 -
      exp(-( ( (-log(u))^r + (-log(v))^r )  )^(1/r) ) * ( ( (-log(u))^r + (-log(v))^r )  )^(1/r -1) * (-log(u))^(r-1) *(-1/u)^2
    }
    
    dcop22 <- function(u, v, r) {
      exp(-( ( (-log(u))^r + (-log(v))^r )  )^(1/r) ) * (( ( (-log(u))^r + (-log(v))^r )  )^(1/r -1) * (-log(v))^(r-1) *(-1/v))^2 -
      exp(-( ( (-log(u))^r + (-log(v))^r )  )^(1/r) ) * (1/r -1) *( ( (-log(u))^r + (-log(v))^r )  )^(1/r -2) * r * ((-log(v))^(r-1) *(-1/v))^2 -
      exp(-( ( (-log(u))^r + (-log(v))^r )  )^(1/r) ) * ( ( (-log(u))^r + (-log(v))^r )  )^(1/r -1) * (r-1) * (-log(v))^(r-2) *(-1/v)^2 -
      exp(-( ( (-log(u))^r + (-log(v))^r )  )^(1/r) ) * ( ( (-log(u))^r + (-log(v))^r )  )^(1/r -1) * (-log(v))^(r-1) *(-1/v)^2
    }
  } else if (copula_name == "gaussian") {
    # Gaussian copula偏导
  dcop1 <- function(u, v, r) {
    pnorm((qnorm(v) - r * qnorm(u)) / sqrt(1 - r^2))
  }
  
  dcop2 <- function(u, v, r) {
    pnorm((qnorm(u) - r * qnorm(v)) / sqrt(1 - r^2))
  }
  
  dcop11 <- function(u, v, r) {
    (sqrt(1 / (1 - r^2))) * (-r) * dnorm((qnorm(v) - r * qnorm(u)) / sqrt(1 - r^2)) / dnorm(qnorm(u))
  }
  
  dcop22 <- function(u, v, r) {
    (sqrt(1 / (1 - r^2))) * (-r) * dnorm((qnorm(u) - r * qnorm(v)) / sqrt(1 - r^2)) / dnorm(qnorm(v))
  }
  
  }
  
  
  return(list(dcop1 = dcop1, dcop2 = dcop2, dcop11 = dcop11, dcop22 = dcop22))

  }


