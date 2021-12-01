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
# build dataset without censored data
data <- cancer %>%
  filter(status == 2) %>%
  drop_na()
# identify covariate labels
xtags <- data %>%
  dplyr::select(-status, -time, -inst) %>%
  colnames()
# retrieve the number of covariates
D <- length(xtags)
# build model formula of the form
formula <- formula(paste("time ~ (",
                         paste0(xtags, collapse = " + "), 
                         "| inst)"))
# define priors over regressors and shape
prior <- c(
  prior_string("normal(0, 1)", class = "sd"),
  prior_string("cauchy(0, 5)", class = "shape")
)
# generate model stan code
make_stancode(formula, data, weibull, prior)
# fit hierarchical GLM model with BRMS, expliciting a non-exponential family
weibull_hier <- brm(
  formula,
  data = data,
  prior = prior,
  family = weibull(link = "log", link_shape = "log"),
  iter = 50000,
  cores = parallel::detectCores()
)
# have a look at some convergence diagnostics
summary(weibull_hier)

# plot rhat
mcmc_rhat(rhat = rhat(weibull_hier))

# plot posterior distribution of parameters
posterior <- as.array(weibull_hier)
bayesplot::mcmc_areas(
  posterior,
  prob = 0.8,
  # 80% intervals
  prob_outer = 0.99,
  # 99%
  point_est = "mean"
)

# perform approximate loo and psis-loo
log_lik <- extract_log_lik(weibull_hier, merge_chains = FALSE)
# estimate the PSIS effective sample size and Monte Carlo error
r_eff <- relative_eff(exp(weibull_hier), cores = parallel::detectCores())
# compute loo
loo(weibull_hier, cores = parallel::detectCores())
# compute waic
waic(weibull_hier, cores = parallel::detectCores())
