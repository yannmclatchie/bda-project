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
# p(time) = sum(covariates)
formula <- formula(paste("time ~",
                         paste0(xtags, collapse = " + ")))
# define prior assumption of variate
prior <- set_prior("normal(0, 10)", class = "b")
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

# import latent projpred helper functions
setwd("/Users/yannmclatchie/Desktop/Aalto")
source("latent_funs.R")
# extract posterior draws of latent predictor
eta_post_draws <- extract_eta(fit, data)

# define new reference model for latent predictor
latent_ref <- fit_latent(fit,
                         eta_post_draws,
                         data)

# perform latent projection predictive variable selection on model
vs <- varsel(
  latent_ref,
  nterms_max = D,
  method = "l1",
  ndraws_pred = 100
)
# selection order of the variables
solution_terms(vs)
# plot varsel results
plot(vs, stats = c('elpd', 'rmse'))
plot(vs, stats = c('elpd'))
plot(vs, stats = c('elpd', 'rmse'), deltas = TRUE)
# compare predictive power
pred <- proj_linpred(vs, newdata = data, nterms = 4)
ggplot() +
  geom_point(aes(x = colMeans(pred$pred), y = data$time)) +
  geom_abline(slope = 1, color = "red") +
  labs(x = "prediction", y = "y")
