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
# The continuous variables are centered.
X=data
X$status[X$status==1] = 0
X$status[X$status==2] = 1
#X$male = ifelse(X$sex==1,1,0)
#X$female = ifelse(X$sex==2,1,0)
X$age=X$age-mean(X$age)
X$meal.cal=X$meal.cal-mean(X$meal.cal)
X$wt.loss=X$wt.loss-mean(X$wt.loss)

Xcens=X[X$status==0,]
Xcens = Xcens[-c(1,2,3)]
ycens=X$time[X$status==0]

Xobs=X[X$status==1,]
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
  Ncen = nrow(Xcens)
)

# compile and run separate censored model
wm = rstan::stan_model(file = "../stan/weibull_censored.stan")
# print out Stan code
print(wm)
# learn the model with parameters 4 chains, 5000 iterations for each, 2500 iterations for warm-up
weibull_cens = rstan::sampling(wm, data = data_model, iter = 5000)

# Checking convergence
# Rhat values are lower than 1.05. It means that the number of iterations is  enough and samples converged.
rstan::monitor(weibull_cens)
mcmc_rhat(rhat = rhat(weibull_cens)[1:8])

# All chains for $\alpha$ in a single line-plot (does not include warm-up).
betas=extract(weibull_cens)$beta
alphas=extract(weibull_cens)$alpha

plot(alphas[1:2500],type="l",col="red",ylab="alphas")
lines(alphas[2501:5000],col="yellow")
lines(alphas[5001:7500],col="blue")
lines(alphas[7501:10000],col="purple")

# All chains for, for instance, $\beta_{1}$ in a single line-plot (does not include warm-up).

plot(betas[,1][1:2500],type="l",col="red", ylab="betas1")
lines(betas[,1][2501:5000],col="yellow")
lines(betas[,1][5001:7500],col="blue")
lines(betas[,1][7501:10000],col="purple")

# The sequences are mixed, they are not separated from each other, 
# and two variance components (the variance within each sequence and the variance between sequences) don’t differ much. The chains have converged.

# Model distributions and characteristics

# Posterior distribution of alpha-shape parameter
hist(extract(weibull_cens)$alpha, main="Posterior distribution
     of alpha-shape parameter", xlab="",breaks = 50)


# 95% credible intervals of sampled alpha. 
bayesplot::mcmc_areas(as.matrix(weibull_cens), 
                      pars = c("alpha"), prob = 0.95)

# 95% credible intervals of sampled betas. 
bayesplot::mcmc_areas(as.matrix(weibull_cens), 
                      pars = c("beta[1]","beta[4]","beta[5]","beta[7]"), prob = 0.95)

# 95% credible intervals of sampled betas. 
par(mfrow=c(2,2))
bayesplot::mcmc_areas(as.matrix(weibull_cens), 
                      pars = c("beta[2]","beta[3]"), prob = 0.95)
bayesplot::mcmc_areas(as.matrix(weibull_cens), 
                      pars = c("beta[6]"), prob = 0.95)

# Effective sample size ratio - majority is good
# light: between 0.5 and 1 (high)
# mid: between 0.1 and 0.5 (good)
# dark: below 0.1 (low)

neff_ratio(weibull_cens) %>% mcmc_neff()


# Posterior predictive checks
color_scheme_set("brightblue")
yrep=extract(weibull_cens)$ypred
bayesplot::ppc_dens_overlay(yobs, yrep[1:50, ])
bayesplot::ppc_hist(yobs, yrep[1:8, ])
bayesplot::ppc_dens_overlay_grouped(yobs, yrep[1:50, ], group = Xobs[,2])
bayesplot::ppc_freqpoly_grouped(yobs, yrep[1:3,], Xobs[,2]) + yaxis_text()

# Computation of the PSIS-LOO elpd values and the k-values
log_lik_sep = extract_log_lik(weibull_cens, merge_chains = FALSE)
r_eff_sep = relative_eff(exp(log_lik_sep), cores = 2)
loo_sep = loo(log_lik_sep, r_eff = r_eff_sep, cores = 2)
print(loo_sep)
psis_sep = psis(-log_lik_sep, r_eff = r_eff_sep)
plot(psis_sep)
