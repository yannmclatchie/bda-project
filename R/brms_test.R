# install libraries
library(survival)
library(tidyverse)
library(brms)
library(loo)
# set number of cores
options(mc.cores = detectCores())

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
# fit covariate GLM model with BRMS, expliciting a non-exponential family
fit <- brm(
  formula,
  data = data,
  prior = prior,
  family = weibull,
  inits = 0,
  chains = 4,
  iter = 200
)
# have a look at some convergence diagnostics
summary(fit)
