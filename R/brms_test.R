# install libraries
library(survival)
library(tidyverse)
library(brms)
library(loo)
# set number of cores
options(mc.cores = parallel::detectCores())

# read lung cancer data from `survival` library
data("cancer", package = "survival")
# build dataset
data <- cancer %>%
  drop_na()
# identify covariate labels
xtags <- data %>%
  dplyr::select(-status, -time) %>%
  colnames()
# retrieve the number of covariates
D <- length(xtags)
# build model formula of the form
# p(time | status) = sum(covariates)
formula <- formula(paste("time | cens(1 - status) ~",
                         paste0(xtags, collapse = " + ")))
# define prior assumption of variate
prior <- set_prior("normal(0, 1)", class = "b")
# generate model stan code
make_stancode(formula, data, weibull, prior)
# fit covariate GLM model with BRMS, expliciting a non-exponential family
fit <- brm(
  formula,
  data = data,
  prior = prior,
  family = weibull
)
# have a look at some convergence diagnostics
summary(fit)
# plot posterior distribution of parameters
posterior <- as.array(fit)
mcmc_areas(
  posterior,
  prob = 0.8,
  # 80% intervals
  prob_outer = 0.99,
  # 99%
  point_est = "mean"
)
