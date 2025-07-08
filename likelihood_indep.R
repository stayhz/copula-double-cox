#' @title Likelihood_indep
#'
#' @description   This function calculates the neg-loglikelihood
#' for the parametric gamma-frailty model with non-proportional hazard functions
#'
#' @param y Vector of parameters in the form
#' \deqn{y = (ln(a), ln(b), \beta _{shape}, \beta _{scale}, ln(\sigma ^2))}
#' for Weibull hazard function and
#' \deqn{y = (ln(10^3a), ln(10^2b), \beta _{shape}, \beta _{scale}, ln(\sigma ^2))}
#' for Gompertz hazard function, where a and b are slope and shape parameters,
#' \eqn{\beta _{shape}} and \eqn{\beta _{scale}} are the Cox-regression parameters
#'  for shape and scale, respectively, and \eqn{\sigma ^2} is the variance of frailty.
#'  This vector must include at least two parameters, \eqn{ln(a)} and \eqn{ln(b)}.
#' @param D  A data.frame in which to interpret the variables named in the formula.
#' The data set includes the following fields:
#' \enumerate{
#' \item time-to-failure and censoring in the case without left truncation
#' or time-of-start, time-of-failure, and censoring in the case with left truncation at the time of begin
#' (censoring must be either 0 for no event or 1 for event);
#' \item Covariates (continuous or categorical) used in a study (can be empty set).
#' }
#' @param nf The number of continuous and binary factors in the data set D corresponding to the covariates
#' used in the Cox-regression for proportional hazard term.
#' @param nk  The number of continuous and binary factors in the data set D corresponding
#' to the covariates used in the Cox-regression for shape b.
#' @param ncl The number of clusters in the data set D corresponding to the cluster covariate.
#' Is equal to 0 for the fixed-effect model.
#' @param dist Baseline hazard function ('Weibull' or 'Gompertz').
#'
#' @return Neg-loglikelihood
#'
#' @examples
#' \dontrun{
#' LikGenNPH(y, D, nf, nk, ncl, dist)
#' }
#'
#' @export
#'
#'LikGenNPH 是用来计算对数似然函数的，
#'具体是用于两种生存分析模型：Weibull 和 Gompertz 模型。
#'y #参数向量，包含了(log(lambda0),log(k0),beta_shape,beta_scale)
#'D是数据框，包含时间、截尾、截断和cluster信息等。
#'n1: scale部分的协变量数量。
#'n2: shape部分的协变量数量。
#'#'ncl: cluster数量，定义集群效应。
#'dist: 分布类型，取值为 'Weibull' 或 'Gompertz'。
likelihood_indep=function(y,D,p,q,dist,copula_name){
  allpara=y #参数向量
  
  Z=D$time  #生存时间
  Z0=D$trunc #截断，或者说起始时间
  Delta=D$Delta
  Xi=D$Xi
  X=cbind(D$covaX1,D$covaX2)
  W=cbind(D$covaW1,D$covaW2)
  n=length(Z)
  
  ############################# 下面计算似然函数 ########################
  # 第一步，准备好所需要的函数和参数
  if (dist=='Weibull'){
    paraT=allpara[1:(p*2+2)]
    paraC=allpara[(p*2+3):((p+q)*2+4)]
    cop_para=allpara[length(allpara)]
    } else if (dist=='GenGompertz') { 
      paraT=allpara[1:(p*2+3)]
      paraC=allpara[(p*2+4):((p+q)*2+6)]
      cop_para=allpara[length(allpara)]
    } else if (dist=='GenExp') { 
      paraT=allpara[1:(p*2+2)]
      paraC=allpara[(p*2+3):((p+q)*2+4)]
      cop_para=allpara[length(allpara)]
    }
  
  
  
  # cop <- function(x1, x2) {
  #   copulacdf2(copula_name, cbind(x1, x2), cop_para)  # 使用自定义的 copulacdf2 函数
  # }
  
  
  # 调用 dcop1,dcop2 函数
  # partial_derivatives <- func_partialcop(copula_name)
  # dcop1 <- partial_derivatives$dcop1
  # dcop2 <- partial_derivatives$dcop2
  
  

  # 第二步，按照三种case计算对数似然
  Lvec <- numeric(n)  # 记录每个样本对应的概率（似然）
  check1 <- numeric(n)
  check2 <- numeric(n)
  check3 <- numeric(n)
  
  Lvec[1] <- 0  # 第一个 dH 是等于 0 的
  
  
  
  for (i in 1:n) {
    
    z1=Z[i]
    
    margT <- func_marg(z1,X[i,],p,dist,paraT)
    fT <- margT$densf
    u <- margT$distF
    
    margC <- func_marg(z1,W[i,],q,dist,paraC)
    fC <- margC$densf
    v <- margC$distF
    
    
    # fT <- func_density(D,p,p,dist,paraT,ID)
    # fC <- func_density(D,q,q,dist,paraC,ID)
    # u <- func__distribution(D,n1,n2,dist,paraT,ID) #F_T(z)
    # v <- func__distribution(D,n1,n2,dist,paraC,ID) #F_C(z)
    
    # C1 <- dcop1(u, v, cop_para)   # copula的关于第一个分量的一阶偏导
    # C2 <- dcop2(u, v, cop_para)   # copula的关于第二个分量的一阶偏导
    
    # 第一种似然
    L <- 0
    if (Delta[i] == 1 && Xi[i] == 0) {
      L1 <- log(fT) + log(1-v)
      L <- L1   # 第一种似然
      check1[i] <- L1
      
      # 第二种似然
    } else if (Delta[i] == 0 && Xi[i] == 1) {
      L2 <- log(fC) + log(1-u)
      L <- L2   # 第二种似然
      check2[i] <- L2
      
      # 第三种似然
    } else if (Delta[i] == 0 && Xi[i] == 0) {
      # L3 <- log(max(1 - u - v + u*v, 1e-20))
      L3 <- log(1 - u - v + u*v)
      L <- L3  
      check3[i] <- L3
    }
    
    Lvec[i] <- L
  }
  
  Lik <- sum(Lvec)
  
  
  
  Lik = -Lik
  return(Lik)
}




