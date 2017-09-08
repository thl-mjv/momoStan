// A-MOMO algorithm

data {
  int<lower=1> N ; // observations per time series
  int<lower=1> M ; // number of age- and region groups
  int<lower=1> P ; // covariates including intercept
  int<lower=1> NP; // number of predictions
  int<lower=1> NT; // N+NP
  int y  [N,M];
  matrix [NT,M] pop; 
  matrix [NT,P] x;
}
transformed data {
}
parameters {
  matrix[P,M] alpha; 
}
transformed parameters {
  matrix<lower=-10,upper=10>[NT,M] linear_predictor;

  for(m in 1:M) {
    linear_predictor[,m] = x * alpha[,m]+pop[,m];
  }
  
}
model {

  for(m in 1:M) {
    alpha[,m] ~ normal(0.0,1.0);
    y    [,m] ~ poisson_log(linear_predictor[1:N ,m]); 
  }
}
generated quantities {
  matrix[N,M] y_pred;
  for(n in 1:N)
    for(m in 1:M) {
      y_pred[n,m] = poisson_log_rng(linear_predictor[n,m]); 
    }
}
