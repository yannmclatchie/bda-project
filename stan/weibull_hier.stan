data {
  int<lower=0> N;               // number of object
  int<lower=0> M;               // number of features
  int<lower =0> J;              // number of institutions
  vector<lower=0>[N] y;         // target-survival time
  matrix[N, M] X;               // covariates
  int<lower=1, upper=J> ll[N];  // instution labels
}

transformed data {
  matrix[N, M] Xc;  // centered version of X without an intercept
  vector[M] means_X;  // column means of X before centering
  
  // column-center the design matrix for fitting the model
  for (m in 1:M) {
    means_X[m] = mean(X[, m]);
    Xc[, m] = X[, m] - means_X[m];
  }
}

parameters {
  // hyperpriors
  real mu_beta;
  real sigma_beta;
  
  // regressors
  matrix[M, J] beta;            // regressors weights for different institutions
  real<lower=0> alpha;          // shape parameter
}

transformed parameters {
  vector[N] eta;
  vector<lower=0>[N] sigma;
  // compute latent predictor term
  for (n in 1:N) {
    eta[n] = Xc[ll[n], ] * beta[, ll[n]];
  }
  // apply the log inverse link function
  sigma = exp(-eta / alpha);
}

model {
  // hyperpriors
  mu_beta ~ normal(0, 1);
  sigma_beta ~ gamma(1, 1);
  
  // prior over regressor and shape parameters
  for (j in 1:J) {
    beta[, j] ~ normal(mu_beta, sigma_beta);
  }
  alpha ~ cauchy(0, 5);
  
  // fitting model
  for (n in 1:N) {
    y[n] ~ weibull(alpha, sigma[n]);
  }
}

generated quantities {
  // define quantities
  real ypred[N];
  vector[N] log_lik;
  
  // compute predictive distribution for survival time
  for (n in 1:N) {
    ypred[n] = weibull_rng(alpha, sigma[n]);
  }
  
  // log-likelihood
  for (n in 1:N) {
    log_lik[n] = weibull_lpdf(y[n] | alpha, sigma[n]);
  }
}
