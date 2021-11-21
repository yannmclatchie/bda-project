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
library(loo)
library(reshape2)
# set number of cores
options(mc.cores = parallel::detectCores())

# read lung cancer data from `survival` library
data("cancer", package = "survival")
# build dataset from only those non-censored data points
data <- cancer %>%
  filter(status == 2) %>%
  drop_na()
# identify covariate labels and build design matrix
cov_labels <- data %>%
  dplyr::select(-status,-time,-inst) %>%
  colnames()
# build design matrix
X <- as.matrix(data[cov_labels])
# print(dim(X))
# [1] 120   7
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
weibull_model <-
  rstan::sampling(wm, data = weibull_data, iter = 10000)
# verify convergence
rstan::monitor(weibull_model)
mcmc_rhat(rhat = rhat(weibull_model)[1:8])
# take note of important parameters
regression_pars <- c(
  "beta[1]",
  "beta[2]",
  "beta[3]",
  "beta[4]",
  "beta[5]",
  "beta[6]",
  "beta[7]",
  "alpha"
)

# investigate trace plots
color_scheme_set("mix-blue-red")
posterior <- as.array(weibull_model)
mcmc_trace(posterior,
           pars = regression_pars,
           facet_args = list(nrow = 4, labeller = label_parsed))

# plot posterior distribution of regressors
mcmc_areas(
  posterior,
  pars = regression_pars,
  prob = 0.8,
  # 80% intervals
  prob_outer = 0.99,
  # 99%
  point_est = "mean"
)

# plot effective sample size ratio
neff_ratio(weibull_model) %>% 
  mcmc_neff()

# define the inverse link function between the latent predictor and 
# the weibull scale parameter
inv_link <- function(eta) {
  sigma <- exp(-eta / alpha)
  return(sigma)
}

## cross-validation

# perform approximate loo and psis-loo
log_lik <- extract_log_lik(weibull_model, merge_chains = FALSE)
# estimate the PSIS effective sample size and Monte Carlo error
r_eff <- relative_eff(exp(log_lik), cores = parallel::detectCores())
# compute loo
loo(log_lik, r_eff = r_eff, cores = parallel::detectCores())
# compute waic
waic(log_lik)


# we now predict the survival time of the training data and compare 
# the results to the true values
beta <- as.matrix(summary(weibull_model)$summary[, "mean"][1:7])
alpha <-
  as.matrix(summary(weibull_model)$summary[, "mean"]["alpha"])[1]
eta <- as.matrix(X) %*% beta
sigma <- inv_link(eta)
# compute point estimate of y given shape and scale point estimates
mu <- sigma * gamma(1 + 1 / alpha)
plot_df <- melt(data.frame(mu, data$time))
plot_df["patient"] <-
  rep(as.numeric(rownames(data.frame(mu, data$time))), 2)
plot_df %>%
  ggplot(aes(x = patient, y = value, color = variable)) +
  geom_point(aes(color = variable))

