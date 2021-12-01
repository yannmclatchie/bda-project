# clean up environment
rm(list = ls())
gc(reset = TRUE)
options(warn=-1)

# set working directory
setwd("./R")

# install libraries
library(survival)
library(tidyverse)
library(rstan)
library(bayesplot)
library(loo)
library(reshape2)
# set number of cores
options(mc.cores = parallel::detectCores())

# read lung cancer data from "survival' library
data("cancer", package = "survival")

# omittimg NAs
data = cancer %>% na.omit()

# Censoring status is transformed to the column with 0-censored, 1-observed.
# The continuous variables are centered. Some values in institutional variable were changed to form the sequence.
X=data
X$status[X$status==1] = 0
X$status[X$status==2] = 1
X$age=X$age-mean(X$age)
X$meal.cal=X$meal.cal-mean(X$meal.cal)
X$wt.loss=X$wt.loss-mean(X$wt.loss)
X$inst[X$inst==21] = 8
X$inst[X$inst==22] = 9
X$inst[X$inst==26] = 14
X$inst[X$inst==32] = 17

Xcens=X[X$status==0,]
inst_cens=Xcens[1]
Xcens = Xcens[-c(1,2,3)]
ycens=X$time[X$status==0]

Xobs=X[X$status==1,]
inst_obs=Xobs[1]
Xobs = as.matrix(Xobs[-c(1,2,3)])
yobs=X$time[X$status==1]

# build data list for Stan model
data_model = list(
  yobs = yobs,
  Xobs = Xobs,
  N = nrow(Xobs),
  M = ncol(Xobs),
  ycen = ycens,
  Xcen = Xcens,
  Ncen = nrow(Xcens),
  inst_cens = inst_cens[,1],
  inst_obs = inst_obs[,1],
  J = 17
)

# compile and run separate censored model
wm = rstan::stan_model(file = "../stan/weibull_censored_hier.stan")

# print out Stan code
print(wm)

# learn the model with parameters 4 chains, 5000 iterations for each, 2500 iterations for warm-up
iters=5000
weibull_cens = rstan::sampling(wm, data = data_model, iter = iters)

# Checking convergence
# Rhat values are lower than 1.05. It means that the number of iterations is  enough and samples converged.
rstan::monitor(weibull_cens)
mcmc_rhat(rhat = rhat(weibull_cens)[1:25])
weibull_cens
# Effective sample size ratio - majority is good
# light: between 0.5 and 1 (high)
# mid: between 0.1 and 0.5 (good)
# dark: below 0.1 (low)
neff_ratio(weibull_cens) %>% mcmc_neff()

# Computation of the PSIS-LOO elpd values and the k-values
log_lik_hier = extract_log_lik(weibull_cens, merge_chains = FALSE)
r_eff_hier = relative_eff(exp(log_lik_hier), cores = 2)
loo_hier = loo(log_lik_hier, r_eff = r_eff_hier, cores = 2)
print(loo_hier)
psis_hier = psis(-log_lik_hier, r_eff = r_eff_hier)
plot(psis_hier)
