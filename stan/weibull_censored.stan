data {
  int<lower=0> N;    // number of object
  int<lower=0> M;    // number of features
  int<lower=0> Ncen;    // number of object
  vector<lower=0>[N] yobs; // target-survival time
  matrix[N, M] Xobs;    // covariates
  vector<lower=0>[Ncen] ycen; // target-survival time
  matrix[Ncen, M] Xcen;    // covariates
}

parameters {
  vector[M] beta;       // regressors
  real<lower=0> alpha;  // shape parameter
}

transformed parameters {
  // Log inverse link function
  vector<lower=0>[N] sigma = exp(-Xobs*beta / alpha);
}

model {
  // priors
  beta ~ normal(0, 10);
  alpha ~ gamma(1,1);

  // fitting model
  yobs ~ weibull(alpha, sigma);
  
  // Increment log-density with Survival Function
  target += weibull_lccdf(ycen | alpha, exp(-Xcen*beta / alpha));
}

generated quantities {
  // compute predictive distribution for survival time
  real ypred[N] = weibull_rng(alpha, sigma);
  
  // log-likelihood
  vector[N+Ncen] log_lik;
  for (i in 1:N) {
    log_lik[i] = weibull_lpdf(yobs[i] | alpha, sigma[i]);
  }
  // Survival function
  for (j in 1:Ncen){
    log_lik[N+j] = weibull_lccdf(ycen[j] | alpha, exp(-Xcen[j,]*beta / alpha)); 
  }
}
