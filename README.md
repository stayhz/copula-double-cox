
List and Description of Custom R Functions Provided

This file lists and briefly describes all custom R functions used in the simulation and estimation process associated with the manuscript. The project involves five main scenarios, each with its own data generation and model fitting structure.

1. Simulation Functions (by Scenario)

- simumodel.R
  Main script that coordinates the simulations across all scenarios.

- Scenario 1 (Dependent and Independent Data, Fitted by Double-Cox and Independent Models):
  - simu_dep1: Simulates data under dependence and fits the copula-double-Cox model.
  - simu_dep2: Simulates data under dependence and fits the independent model.
  - simu_indep1: Simulates independent data and fits the copula-double-Cox model.
  - simu_indep2: Simulates independent data and fits the independent model.

- Scenario 3 (Over-parameterized Model Test):
  - simu_single: Simulates data from a single-Cox model with copula and fits it using the copula-single-Cox model.

- Scenario 4 (Comparison with Cox Model):
  - simu_indepsingle: Simulates data from an independent single-Cox model.

- Scenario 5 (Misspecified Marginal Distributions):
  - simu_lognorm: Simulates data with log-normal margins and fits with copula-double-Cox model.
  - simu_logncox: Simulates data with log-normal margins and fits with the Cox model.

2. Model Estimation Functions

- method_doublecox.R
  Fits the proposed copula-double-Cox model.

- method_indep.R
  Fits the independent double-Cox model (i.e., no dependence structure).

- method_singlecox.R
  Fits the copula-single-Cox model.

3. Data Generation Functions

- gendata_weibull.R
  Generates survival data with Weibull marginal distributions.

- gendata_genexp.R
  Generates survival data with generalized exponential marginal distributions.

- gendata_lognorm.R
  Generates survival data with log-normal marginal distributions.

4. Likelihood Functions

- likelihood.R
  Provides the full likelihood function for the copula-double-Cox model.

- likelihood_indep.R
  Provides the likelihood function for the independent model (no copula).

5. Supporting Functions for Distribution and Copula Calculations

- func_marg.R
  Computes marginal probability density and cumulative distribution functions.

- func_partialcop.R
  Computes partial derivatives of the copula functions.

- copulastat.R
  Computes copula-related statistics such as Kendall's tau.

- copulacdf2.R
  Computes bivariate copula cumulative distribution functions.

- func_revkendall.R
  Approximates Kendall's tau based on copula parameter estimates.

6. Initial Value Estimation

- easydoublecox.R
  Provides initial parameter estimates for the copula-double-Cox model to accelerate convergence.

7. Simulation Result Analysis

- analyzedata.R
  Analyzes and summarizes results from the copula-double-Cox model.

- analyzedata_indep.R
  Analyzes results from the independent model.

- analyzedata_singlecox.R
  Analyzes results from the copula-single-Cox model.
