data {
  int<lower=0> N; // number of data realisations
  int<lower=0> M; // feature dimensionality
  vector<lower=0>[N] y; // variate
  matrix[N, M] X; // design matrix
}

transformed data {
  // column-center the design matrix for fitting the model
  matrix[N, M] Xc;  // centered version of X without an intercept
  vector[M] means_X;  // column means of X before centering
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
  // compute latent predictor term and Weibull scale parameter
  vector[N] eta = Xc * beta;
  vector[N] sigma;
  for (n in 1:N) {
    // apply the log inverse link function
    sigma[n] = exp(eta[n]);
  }
}

model {
  // prior over regressor and shape parameters
  beta ~ normal(0, 1);
  alpha ~ gamma(1, 1);

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
