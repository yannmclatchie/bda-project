data {
  int<lower=0> N; // number of data realisations
  int<lower=0> M; // feature dimensionality
  vector<lower=0>[N] y; // variate observations
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
  real intercept;  // temporary intercept for centered predictors
  real<lower=0> shape;  // shape parameter
}

model {
  // TODO: add prior over beta and shape parameters
  
  // initialize linear predictor term
  vector[N] mu = intercept + Xc * beta;
  for (n in 1:N) {
    // apply the inverse link function
    mu[n] = exp(mu[n]) / tgamma(1 + 1 / shape);
  }
  // fit model
  y ~ weibull(shape, mu);
}

generated quantities {
  // actual population-level intercept
  real b_intercept = intercept - dot_product(means_X, beta);
}