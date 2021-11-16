data {
  int<lower=0> N; // number of data realisations
  int<lower=0> M; // feature dimensionality
  int<lower=0> y[N]; // variate observations
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
  // Cox regression model parameters
  vector[M] beta; // regressors
  real log_lambda0; // intercept for Cox model
}

transformed parameters {
  // parameter for Poisson distribution
  vector<lower=0>[N] lambda;
  lambda = exp(log_lambda0 + Xc * beta);
}

model {
  // priors
  beta ~ normal(0, 100);
  log_lambda0 ~ normal(0, 100);

  // likelihood
  y ~ poisson(lambda);
}

generated quantities {
  // log-likelihood
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = poisson_lpmf(y[n] | lambda[n]);
  }
}
