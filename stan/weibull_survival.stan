data {
  int<lower=0> N; // number of data realisations
  int<lower=0> M; // feature dimensionality
  vector<lower=0>[N] y; // survival time
  matrix[N, M] X; // design matrix
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
  // GLM parameters
  vector[M] beta; // regressors
  real<lower=0> alpha;  // shape parameter
}

transformed parameters {
  // compute latent predictor term
  vector[N] eta = Xc * beta;
  // apply the log inverse link function
  vector<lower=0>[N] sigma = exp(-eta / alpha);
}

model {
  // prior over regressor and shape parameters
  beta ~ normal(0, 10);
  alpha ~ cauchy(0, 10);

  // fit model
  y ~ weibull(alpha, sigma);
}

generated quantities {
  // compute predictive distribution for survival time
  real ypred[N] = weibull_rng(alpha, sigma);
  
  // log-likelihood
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = weibull_lpdf(y[n] | alpha, sigma[n]);
  }
}
