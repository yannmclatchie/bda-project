# install libraries
library(survival)
library(tidyverse)
library(rstan)
library(bayesplot)
library(loo)
library(splines2)
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
formula <- formula(paste("time ~ ",
                         paste0(xtags, collapse = " + ")))
# define priors over regressors and shape
prior <- get_prior(
  formula,
  data,
  cox
)
prior <- c(
  prior_string("normal(0, 1)", class = "b")
)
# generate model stan code
make_stancode(formula, data, cox, prior)
# fit hierarchical GLM model with BRMS, expliciting a non-exponential family
cox_model <- brm(
  formula,
  data = data,
  prior = prior,
  family = cox,
  iter = 50000,
  cores = parallel::detectCores()
)
# have a look at some convergence diagnostics
summary(cox_model)

# plot rhat
mcmc_rhat(rhat = rhat(cox_model))

# plot posterior distribution of parameters
posterior <- as.array(cox_model)
bayesplot::mcmc_areas(
  posterior,
  prob = 0.8,
  # 80% intervals
  prob_outer = 0.99,
  # 99%
  point_est = "mean"
)

# perform approximate loo and psis-loo
log_lik <- extract_log_lik(cox_model, merge_chains = FALSE)
# estimate the PSIS effective sample size and Monte Carlo error
r_eff <- relative_eff(exp(cox_model), cores = parallel::detectCores())
# compute loo
loo(cox_model, cores = parallel::detectCores())
# compute waic
waic(cox_model, cores = parallel::detectCores())
