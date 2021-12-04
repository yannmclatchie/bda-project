data {
  int<lower=0> N;    // number of object
  int<lower=0> M;    // number of features
  vector<lower=0>[N] y; // target-survival time
  matrix[N, M] X;    // covariates
  int <lower =0> J; // number of institutions
  int<lower=1,upper=J> inst[N]; // inst labels
}

parameters {
  vector[M] mu; // mean hyperior over regressors
  vector<lower=0>[M] sigma; // variance hyperior over regressors
  matrix[M,J] beta;       // regressors weights for different institutions
  real<lower=0> alpha;  // shape parameter
}

transformed parameters {
  // compute latent predictor term
  vector[N] eta = Xc * beta;
  // apply the log inverse link function
  vector<lower=0>[N] sigma = exp(-eta / alpha);
}

transformed parameters {
  // Log inverse link function
  vector<lower=0>[N] sigma;
  for (i in 1:N){
    sigma[i] = exp(dot_product(-Xobs[i,], beta[,inst_obs[i]]) / alpha);
  }

}

model {
  // priors
  for(i in 1:M){
    mu0 ~ normal(0, 1);
  }
  for(i in 1:M){
    sigma0 ~ normal(0, 1);
  }
  
  for(i in 1:M){
    for(j in 1:J){
      beta[i,j] ~ normal(mu0[i], sigma0[i]);
    }
  }
  alpha ~ cauchy(0,5);

  // fitting model
  yobs ~ weibull(alpha, sigma);
  
  // Increment log-density with Survival Function
  for (i in 1:Ncen){
    target += weibull_lccdf(ycen[i] | alpha, exp(-Xcen[i,]*beta[,inst_cens[i]]/alpha));
  }
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
    log_lik[N+j] = weibull_lccdf(ycen[j] | alpha, exp(-Xcen[j,]*beta[,inst_cens[j]] / alpha)); 
  }
}
