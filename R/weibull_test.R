# install libraries
library(survival)
library(tidyverse)
library(rstan)
library(bayesplot)
library(loo)
# set number of cores
options(mc.cores = parallel::detectCores())

# read lung cancer data from `survival` library
data("cancer", package = "survival")
# build dataset
data <- cancer %>%
  drop_na()
# identify covariate labels and build design matrix
cov_labels <- data %>%
  dplyr::select(-status, -time) %>%
  colnames()

X <- as.data.frame(data[cov_labels])
# print(dim(X))
# [1] 167   8
y <- data$time
# build data list for Stan model
weibull_data = list(
  y = y, X = X, N = length(y), M = ncol(X)
)

# compile and run seperate model
wm <- rstan::stan_model(file = "../stan/weibull_test.stan")
# print out Stan code
print(wm)
# learn the model parameters
weibull_model <- rstan::sampling(wm, data = weibull_data)
# verify convergence
rstan::monitor(weibull_model)

# plot posterior distribution of regressors
posterior <- as.array(weibull_model)
mcmc_areas(
  posterior,
  pars = c(
    "beta[1]",
    "beta[2]",
    "beta[3]",
    "beta[4]",
    "beta[5]",
    "beta[6]",
    "beta[7]",
    "beta[8]"
  ),
  prob = 0.8,
  # 80% intervals
  prob_outer = 0.99,
  # 99%
  point_est = "mean"
)
