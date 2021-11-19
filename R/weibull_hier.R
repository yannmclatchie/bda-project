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
# build dataset
data <- cancer %>%
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
# define prior assumption of variate
prior <- get_prior(formula,
          data,
          family = weibull)
# generate model stan code
make_stancode(formula, data, weibull, prior)
# fit hierarchical GLM model with BRMS, expliciting a non-exponential family
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
bayesplot::mcmc_areas(
  posterior,
  prob = 0.8,
  # 80% intervals
  prob_outer = 0.99,
  # 99%
  point_est = "mean"
)
