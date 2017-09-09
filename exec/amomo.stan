// A-MOMO algorithm

data {
  int<lower=1> N ; // observations per time series
  int<lower=1> M ; // number of age- and region groups
  int<lower=1> P ; // covariates including intercept
  int<lower=1> NP; // number of predictions
  int<lower=1> NT; // N+NP
  real<lower=0,upper=1>  penalty[3]; // penalty control terms for estimates
  int y  [N,M];
  matrix [NT,M] pop;
  matrix [NT,P] x;
}
transformed data {
}
parameters {
  matrix[P,M] alpha_eps;
  vector[P] alpha_group; // (dis)similarity between the M groups
  vector[M] alpha_param; // (dis)similarity between the P parameters
  vector<lower=0> [3] shrinkage;
}
transformed parameters {
  matrix<lower=-10,upper=10>[NT,M] linear_predictor;
  matrix[P,M] alpha_real;
  for(m in 1:M) {
    for(p in 1:P) {
      alpha_real[p,m]=penalty[1]*alpha_group[p]+penalty[2]*alpha_param[m]+penalty[3]*alpha_eps[p,m];
    }
    linear_predictor[,m] = x * alpha_real[,m]+pop[,m];
  }
}
model {
  shrinkage   ~ normal(0.0,1.0);
  alpha_param ~ normal(0.0,shrinkage[1]);
  alpha_group ~ normal(0.0,shrinkage[2]);
  for(m in 1:M) {
    alpha_eps[,m] ~ normal(0.0,shrinkage[3]);
    y        [,m] ~ poisson_log(linear_predictor[1:N ,m]);
  }
}
generated quantities {
  matrix[NT,M] y_pred;
  for(n in 1:NT)
    for(m in 1:M) {
      y_pred[n,m] = poisson_log_rng(linear_predictor[n,m]);
    }
}
