#' @title func_marg
#'
#' @description   This function calculates marginal survival
#'
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
#' to the covariates used in the Cox-regression for shape \eqn{b}.
#' @param ncl The number of clusters in the data set D corresponding to the cluster covariate.
#' Is equal to 0 for the fixed-effect model).
#' @param dist Baseline hazard function ('Weibull' or 'Gompertz').
#' @param vpar The sample of the vector-row-parameters in matrix form.
#' @param ID The order number of an object in the data set.
#'
#' @return Marginal survival for the object with order number ID
#'
#' @examples
#' \dontrun{
#' msurv(D,nf,nk,ncl,dist,vpar,ID)
#' }
#'
#' @export
#'这段代码定义了一个名为 msurv 的函数，用于计算生存概率。
#'函数的输入包括数据集、不同参数的数量、分布类型、参数向量以及特定行的ID。
#'D: 数据集，通常是一个数据框。
#'nf: 固定效应协变量的数量。
#'nk: 交互效应协变量的数量。
#'ncl: 聚类效应的数量。
#'dist: 分布类型，可以是 'Weibull' 或 'Gompertz'。
#'para: 参数向量，包含了模型的参数。
#'ID: 指定行的索引，用于从数据集中提取特定行的数据。
#'
func_marg=function(Z,X,p,dist,para){
  # 对于Weibull分布，a1和 a2 分别是scale参数和shape参数的baseline值。
  # #对于Gompertz分布，参数需要额外的比例因子进行缩放。
  if (dist=='Weibull'){
    a1=exp(para[1]) # 可保证baseline是正的
    a2=exp(para[2])} else if (dist=='GenGompertz') {
      a1=exp(para[1])
      a2=exp(para[2])
      parav=1e-1*exp(para[length(para)])
    }else if (dist=='GenExp'){
      a1=exp(para[1]) # 可保证baseline是正的
      a2=exp(para[2])} 
  

  
  n1=length(X)
  n2=n1
  
  
  #计算Cox回归模型中的比例风险部分
  #Coxshape 和 Coxscale 分别用于存储比例风险部分和固定效应部分的值。
  Coxshape=0
  Coxscale=0
  if (p>0){
    for (i in 1:n1){
      Coxscale=Coxscale+X[i]*para[2+i]}
    for (i in 1:n2){
      Coxshape=Coxshape+X[i]*para[2+n1+i]}
  }
  Coxshape=exp(Coxshape)
  Coxscale=exp(Coxscale)
  
  
  #调整shape参数和提取观测时间
  parascale=a1*Coxscale
  parashape=a2*Coxshape

  #计算 Hfull1 和 Hfull0，分别对应生存时间和截断时间。
  # 第一种weibull分布：基于S=exp(-(t/lambda)^v)
  if (dist=='Weibull'){
    # Hfull1= Coxscale*(Z/a1)^parashape
    # Hfull0= Coxscale*(z0/a1)^parashape
    # cumh=Hfull1-Hfull0
    cumh <- Coxscale*((Z/a1)^parashape)
    mufull <- Coxscale*parashape*(Z^(parashape-1))/(a1^parashape) #hazard function
    survS <- exp(-cumh)
    distF <- (1-exp(-cumh))
    densf <- mufull*survS
  }
  
  
  # # 第二种weibull分布：基于S=exp(-(t^v)/lambda)
  # if (dist=='Weibull'){
  #   # Hfull1= Coxscale*(Z/a1)^parashape
  #   # Hfull0= Coxscale*(z0/a1)^parashape
  #   # cumh=Hfull1-Hfull0
  #   cumh <- (Z^parashape)/parascale
  #   mufull <- parashape*(Z^(parashape-1))/parascale #hazard function
  #   survS <- exp(-cumh)
  #   distF <- (1-exp(-cumh))
  #   densf <- mufull*survS
  # }
  
  
  # if (dist=='Gompertz'){
  #   # Hfull1= (Cox*lambda0/k0)*(exp(k0*Z)-1)
  #   # Hfull0= (Cox*lambda0/k0)*(exp(k0*z0)-1)
  #   # mufull1= (Cox*lambda0)*exp(k0*Z)
  #   cumh <- (Cox*lambda0/k0)*(exp(k0*Z)-1)
  #   mufull1 <- (Coxscale*a1)*exp(parashape*x1)
  #   survS <- exp(-cumh)
  #   distF <- (1-exp(-cumh))
  #   densf=mufull1*msurv
  # }
  
  
  # if (dist=='GenGompertz'){
  #   distF <- (1-exp(parascale*(1-exp(parav*Z))/parav))^parashape
  #   cumh <- -log(1- distF)
  #   survS <- 1- distF
  #   termexp <- exp(parascale*(1-exp(parav*Z))/parav)
  #   densf <- parashape*parascale*exp(parav*Z)*termexp*(1-termexp)^(parashape-1)
  #   mufull <- densf/survS    #hazard function
  # }
  
  if (dist=='GenGompertz'){
  termexp <- exp(parascale*(1-exp(parav*Z))/parav)
  distF <- (1-termexp)^parashape
  survS <- 1-distF
  cumh <- -log(survS)
  densf <- parascale *parashape *exp(parav*Z) *termexp *(1-termexp)^(parashape-1)
  mufull <- densf/survS
  }
  
  
  if (dist=='GenExp'){
    termexp <- exp(-parascale*Z)
    distF <- (1-termexp)^parashape
    survS <- 1-distF
    cumh <- -log(survS)
    densf <- parascale *parashape *termexp *(1-termexp)^(parashape-1)
    mufull <- densf/survS
  }
  
  
  #返回概率密度：# return(densf)
  return(list(distF = distF, cumh = cumh, survS = survS, densf = densf,mufull=mufull ))
  
}

#该函数 msurv 用于计算指定数据行的生存概率，基于 Weibull 或 Gompertz 分布模型，
#考虑固定效应、交互效应和聚类效应。
#函数的核心是根据不同的生存分布模型计算累积风险函数，然后通过累积风险函数计算生存概率。






func_margsingle=function(Z,X,p,dist,para){
  # 对于Weibull分布，a1和 a2 分别是scale参数和shape参数的baseline值。
  # #对于Gompertz分布，参数需要额外的比例因子进行缩放。
  if (dist=='Weibull'){
    a1=exp(para[1]) # 可保证baseline是正的
    a2=exp(para[2])} else if (dist=='GenGompertz') {
      a1=exp(para[1])
      a2=exp(para[2])
      parav=1e-1*exp(para[length(para)])
    }else if (dist=='GenExp'){
      a1=exp(para[1]) # 可保证baseline是正的
      a2=exp(para[2])} 
  
  
  
  n1=length(X)
  n2=n1
  
  
  #计算Cox回归模型中的比例风险部分
  #Coxshape 和 Coxscale 分别用于存储比例风险部分和固定效应部分的值。
  # Coxshape=0
  Coxscale=0
  if (p>0){
    for (i in 1:n1){
      Coxscale=Coxscale+X[i]*para[2+i]}
  #   for (i in 1:n2){
  #     Coxshape=Coxshape+X[i]*para[2+n1+i]}
  }
  # Coxshape=exp(Coxshape)
  Coxscale=exp(Coxscale)
  
  
  
  
  
  
  #调整shape参数和提取观测时间
  parascale=a1*Coxscale
  # parashape=a2*Coxshape
  parashape=a2
  
  #计算 Hfull1 和 Hfull0，分别对应生存时间和截断时间。
  # 第一种weibull分布：基于S=exp(-(t/lambda)^v)
  if (dist=='Weibull'){
    # Hfull1= Coxscale*(Z/a1)^parashape
    # Hfull0= Coxscale*(z0/a1)^parashape
    # cumh=Hfull1-Hfull0
    cumh <- Coxscale*((Z/a1)^parashape)
    mufull <- Coxscale*parashape*(Z^(parashape-1))/(a1^parashape) #hazard function
    survS <- exp(-cumh)
    distF <- (1-exp(-cumh))
    densf <- mufull*survS
  }
  
  
  # # 第二种weibull分布：基于S=exp(-(t^v)/lambda)
  # if (dist=='Weibull'){
  #   # Hfull1= Coxscale*(Z/a1)^parashape
  #   # Hfull0= Coxscale*(z0/a1)^parashape
  #   # cumh=Hfull1-Hfull0
  #   cumh <- (Z^parashape)/parascale
  #   mufull <- parashape*(Z^(parashape-1))/parascale #hazard function
  #   survS <- exp(-cumh)
  #   distF <- (1-exp(-cumh))
  #   densf <- mufull*survS
  # }
  
  
  # if (dist=='Gompertz'){
  #   # Hfull1= (Cox*lambda0/k0)*(exp(k0*Z)-1)
  #   # Hfull0= (Cox*lambda0/k0)*(exp(k0*z0)-1)
  #   # mufull1= (Cox*lambda0)*exp(k0*Z)
  #   cumh <- (Cox*lambda0/k0)*(exp(k0*Z)-1)
  #   mufull1 <- (Coxscale*a1)*exp(parashape*x1)
  #   survS <- exp(-cumh)
  #   distF <- (1-exp(-cumh))
  #   densf=mufull1*msurv
  # }
  
  
  # if (dist=='GenGompertz'){
  #   distF <- (1-exp(parascale*(1-exp(parav*Z))/parav))^parashape
  #   cumh <- -log(1- distF)
  #   survS <- 1- distF
  #   termexp <- exp(parascale*(1-exp(parav*Z))/parav)
  #   densf <- parashape*parascale*exp(parav*Z)*termexp*(1-termexp)^(parashape-1)
  #   mufull <- densf/survS    #hazard function
  # }
  
  if (dist=='GenGompertz'){
    termexp <- exp(parascale*(1-exp(parav*Z))/parav)
    distF <- (1-termexp)^parashape
    survS <- 1-distF
    cumh <- -log(survS)
    densf <- parascale *parashape *exp(parav*Z) *termexp *(1-termexp)^(parashape-1)
    mufull <- densf/survS
  }
  
  
  if (dist=='GenExp'){
    termexp <- exp(-parascale*Z)
    distF <- (1-termexp)^parashape
    survS <- 1-distF
    cumh <- -log(survS)
    densf <- parascale *parashape *termexp *(1-termexp)^(parashape-1)
    mufull <- densf/survS
  }
  
  
  #返回概率密度：# return(densf)
  return(list(distF = distF, cumh = cumh, survS = survS, densf = densf,mufull=mufull ))
  
}




