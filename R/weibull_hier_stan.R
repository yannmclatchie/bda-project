# clean up environment
rm(list = ls())
gc(reset = TRUE)

# set working directory
setwd("~/Desktop/Aalto/BDA/bda-project/R")

# install libraries
library(survival)
library(tidyverse)
library(rstan)
library(bayesplot)
# set number of cores
options(mc.cores = parallel::detectCores())

# read lung cancer data from `survival` library
data("cancer", package = "survival")
# build dataset from only those non-censored data points
data <- cancer %>%
  filter(status == 2) %>%
  filter(inst %in% c(1, 12, 13, 3, 11)) %>%
  drop_na()
# identify covariate labels and build design matrix
cov_labels <- data %>%
  dplyr::select(-status,-time,-inst) %>%
  colnames()
# build design matrix
X <- as.matrix(data[cov_labels])
# institution indicators as numeric factors
ll <- as.numeric(as.factor(data$inst))
# response variable
y <- data$time
# build data list for Stan model
hier_data = list(
  y = y,
  X = X,
  ll = ll,
  N = length(y),
  M = ncol(X),
  J = length(unique(ll))
)

# compile and run seperate model
hwm <- rstan::stan_model(file = "../stan/weibull_hier.stan")
# print out Stan code
print(hwm)
# learn the model parameters
hier_model <- rstan::sampling(hwm, data = hier_data, iter = 10000)
# verify convergence
mcmc_rhat(rhat = rhat(hier_model))
rstan::monitor(hier_model)

# perform approximate loo and psis-loo
hwm_log_lik <- extract_log_lik(hier_model, merge_chains = FALSE)
# estimate the PSIS effective sample size
hwm_r_eff <- relative_eff(exp(hwm_log_lik), cores = parallel::detectCores())
# compute loo
hwm_loo <- loo(hwm_log_lik, r_eff = hwm_r_eff, cores = parallel::detectCores())
# compute waic
hwm_waic <- waic(hwm_log_lik, cores = parallel::detectCores())
