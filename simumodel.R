####################################################################################
####################### All Simulation Code for Model Fitting ######################
####################################################################################

library(copula)
library(parallel)
library(pbapply)
library(mvtnorm)
library(ucminf)
library(MASS)       # For multivariate normal distributions
library(progress)
library(survival)
library(optimx)
library(writexl)
library(ggplot2)

# Load all required user-defined functions
source("func_marg.R"); source("func_partialcop.R"); source("likelihood.R")
source("copulastat.R"); source("copulacdf2.R"); source("func_revkendall.R")
source("easydoublecox.R"); source("gendata_weibull.R"); source("gendata_genexp.R")
source("method_doublecox.R"); source("method_indep.R"); source("method_singlecox.R")
source("simu_indep1.R"); source("simu_indep2.R"); source("simu_dep1.R"); source("simu_dep2.R")
source("simu_single.R"); source("analyzedata.R"); source("analyzedata_indep.R")
source("analyzedata_singlecox.R"); 
source("gendata_lognorm.R"); source("simu_lognorm.R"); source("simu_logncox.R");source("simu_indepsingle.R")

####################### Scenario 1:  Comparison with independent model #######################

## Data generated from independent or dependent models
## Compare proposed model with independent model

### Parameter setup: Weibull model
truepara <- c(5, 0.6, 0.5, -1.2, 0.1, -0.2, 11, 1.2, 0.6, -1.4, 0.1, 0.2, 5)
pA <- 40
n <- 500; num1 <- 500; p <- 2; q <- 2; dist <- 'Weibull'
copula_name <- c("frank", "frank")
numseed0 <- 3
numCores <- detectCores()
cl <- makeCluster(numCores)

### (Alternative) Parameter setup: GenExp model
truepara <- c(0.04, 0.4, 0.5, -1.2, 0.1, -0.2, 0.07, 0.8, 0.6, -1.4, 0.1, 0.2, 5)
dist <- 'GenExp'

##### Case 1: Data from independent model, estimated using copula-double-Cox
results_indep1 <- parLapply(cl, 1:num1, simu_indep1, n, p, q, truepara, copula_name, dist, numseed0, pA)
anares1 <- analyzedata(results_indep1, truepara, copula_name)
anares1
save(list = ls(), file = "simudata1.RData")

##### Case 2: Data from independent model, estimated using independent model
results_indep2 <- parLapply(cl, 1:num1, simu_indep2, n, p, q, truepara, copula_name, dist, numseed0, pA)
anares2 <- analyzedata(results_indep2, truepara, copula_name)
anares2

##### Case 3: Data from dependent model, estimated using copula-double-Cox
results_dep1 <- parLapply(cl, 1:num1, simu_dep1, n, p, q, truepara, copula_name, dist, numseed0, pA)
anares3 <- analyzedata(results_dep1, truepara, copula_name)
anares3

##### Case 4: Data from dependent model, estimated using independent model
results_dep2 <- parLapply(cl, 1:num1, simu_dep2, n, p, q, truepara, copula_name, dist, numseed0, pA)
anares4 <- analyzedata(results_dep2, truepara, copula_name)
anares4

### Copy results to clipboard for Excel or other tools
# anares <- anares4
# anares <- anares3
# anares <- anares2
anares <- anares1

combined_data <- data.frame(
  bias_exp = anares$bias_exp,
  sdvec    = anares$sdvec,
  ASE      = anares$ASE,
  CPvec    = anares$CPvec
)

write.table(
  combined_data,
  "clipboard",  # For Windows clipboard output
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

####################################################################################
########################### Scenario 2: Sensitivity test of Copula Structure ####################
## Data generated from Frank copula, estimated using Gumbel/Gaussian copulas

truepara <- c(5, 0.6, 0.5, -1.2, 0.1, -0.2, 11, 1.2, 0.6, -1.4, 0.1, 0.2, 5)
dist <- 'Weibull'
copula_namegum <- c("frank", "gumbel")
copula_namegau <- c("frank", "gaussian")

results_miscop1 <- parLapply(cl, 1:num1, simu_dep1, n, p, q, truepara, copula_namegum, dist, numseed0, pA)
anares_miscop1 <- analyzedata(results_miscop1, truepara, copula_namegum)
anares_miscop1

results_miscop2 <- parLapply(cl, 1:num1, simu_dep1, n, p, q, truepara, copula_namegau, dist, numseed0, pA)
anares_miscop2 <- analyzedata(results_miscop2, truepara, copula_namegau)
anares_miscop2

####################################################################################
########################### Scenario 3: Over-parameterization ######################
## Data generated from single-Cox model, estimated using both full and simplified models

# copula-double-Cox
truepara <- c(5, 0.6, 0.5, -1.2, 0, 0, 11, 1.2, 0.6, -1.4, 0, 0, 5)
dist <- 'Weibull'
results_over <- parLapply(cl, 1:num1, simu_dep1, n, p, q, truepara, copula_name, dist, numseed0, pA)
anares_over <- analyzedata(results_over, truepara, copula_name)
anares_over

# copula-single-Cox
results_single <- parLapply(cl, 1:num1, simu_single, n, p, q, truepara, copula_name, dist, numseed0, pA)
anares_single <- analyzedata_singlecox(results_single, truepara, copula_name)
anares_single

####################################################################################
#################### Scenario 4: Comparison with classical Cox model ################
## Data generated from independent single-Cox model

truepara <- c(5, 0.6, 0.5, -1.2, 0, 0, 11, 1.2, 0.6, -1.4, 0, 0, 0)

# copula-double-Cox
results_over <- parLapply(cl, 1:num1, simu_dep1, n, p, q, truepara, copula_name, dist, numseed0, pA)
anares_over <- analyzedata(results_over, truepara, copula_name)
anares_over

# Classical Cox model
results_cox <- parLapply(cl, 1:num1, simu_indepsingle, n, p, q, truepara, copula_name, dist, numseed0, pA)
anares_cox <- analyzedata_indep(results_cox, truepara, copula_name)
anares_cox

####################################################################################
########## Scenario 5: Sensitivity to misspecified marginal distributions ###########
## Data generated from log-normal but fitted with incorrect distribution

truepara <- c(5, 0.6, 0.5, -1.2, 0, 0, 11, 1.2, 0.6, -1.4, 0, 0, 0)
dist <- 'Weibull'

# copula-double-Cox
results_lognorm <- parLapply(cl, 1:num1, simu_lognorm, n, p, q, truepara, copula_name, dist, numseed0, pA)
anares_lognorm <- analyzedata(results_lognorm, truepara, copula_name)
anares_lognorm

# Cox model
results_logncox <- parLapply(cl, 1:num1, simu_logncox, n, p, q, truepara, copula_name, dist, numseed0, pA)
anares_logncox <- analyzedata_indep(results_logncox, truepara, copula_name)
anares_logncox

