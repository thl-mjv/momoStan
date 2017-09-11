// FluMOMO algorithm

data {
  int<lower=1> N ; // observations per time series
  int<lower=1> M ; // number of age- and region groups
  int<lower=1> P ; // covariates including intercept
  int<lower=1> Q ; // number of covariates
  int y  [N,M]; // observations
  matrix [N,M] pop; 
  matrix [N,P] x;
  real z [N,M,Q];
}
transformed data {
  int yy [N,M];
  yy=y;
}
parameters {
  matrix[P,M] alpha_baseline;
  matrix[P,M] alpha_baseline_null;  
  matrix<lower=0>[Q,M] alpha_covariates; 
}
transformed parameters {
  matrix<lower=-10,upper=10>[N,M] baseline;
  matrix<lower=-10,upper=10>[N,M] baseline_null;
  matrix<lower=-10,upper=10>[N,M] covariate_effects;
  
  for(m in 1:M) {
    baseline         [,m] = x * alpha_baseline     [,m]+pop[,m];
    baseline_null    [,m] = x * alpha_baseline_null[,m]+pop[,m];
    covariate_effects[,m] = to_matrix(z[,m,])*alpha_covariates[,m];
  }
  
}
model {
  for(m in 1:M) {
    alpha_baseline     [,m] ~ normal(0.0,1.0);
    alpha_baseline_null[,m] ~ normal(0.0,1.0);
    alpha_covariates   [,m] ~ normal(0.0,1.0);
    y    [,m] ~ poisson_log(baseline     [1:N ,m]+covariate_effects[1:N,m]); 
    yy   [,m] ~ poisson_log(baseline_null[1:N ,m]); 
  }
}
generated quantities {
  matrix[N,M] pred_baseline;
  matrix[N,M] pred_baseline_null;
  real pred_effects [N,M,Q] ;
  matrix[N,M] pred_full;
  for(n in 1:N)
    for(m in 1:M) {
      pred_baseline     [n,m] = poisson_log_rng(baseline     [n,m]); 
      pred_baseline_null[n,m] = poisson_log_rng(baseline_null[n,m]);
      pred_full         [n,m] = poisson_log_rng(baseline     [n,m]+covariate_effects[n,m]); 
      for(q in 1:Q) {
	pred_effects[n,m,q] = poisson_log_rng(baseline[n,m]+z[n,m,q]*alpha_covariates[q,m]); // OK?
      }
    }
}
