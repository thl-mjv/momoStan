data {
  int<lower=1> N; // observations per time series
  int<lower=1> Na; // number of age- and region groups
  int<lower=1> P ; // covariates including intercept
  int<lower=1> Q ; // temperature
  int y  [N,Na];
  matrix [N,Na] pop; 
  matrix [N,P ] x;
  real z[N,Q,Na];
}
transformed data {
  int yy [N,Na]; // copy of the response
  for(i in 1:N) for(j in 1:Na) yy[i,j] = y[i,j];
}
parameters {
  matrix[P,Na] alphay; // baseline characteristics
  matrix[P,Na] betay;  // baseline characteristics
  matrix<lower=0>[Q,Na] alphat; // effect of temp
  vector<lower=0>[Na] sigma; // shrinkage for the temp effects
}

transformed parameters {
  matrix<lower=-10,upper=10>[N,Na] fy ;
  matrix<lower=-10,upper=10>[N,Na] fyx;
  matrix<lower=-10,upper=10>[N,Na] fyz;

  for(i in 1:Na) {
    fyx[,i] = x * alphay[,i] + pop[,i] + to_matrix(z[,,i])*alphat[,i];
    fyz[,i] = x * alphay[,i] + pop[,i];
    fy [,i] = x * betay [,i] + pop[,i];
  }  
}
model {
  for(i in 1:Na) {
    sigma [ i] ~ normal(0.0,0.1);
    alphay[,i] ~ normal(0.0,1.0);
    betay [,i] ~ normal(0.0,1.0);
    alphat[,i] ~ normal(0.0,sigma[i]);
    y  [,i] ~ poisson_log(fyx[,i]);
    yy [,i] ~ poisson_log(fy [,i]);
  }
}

generated quantities {
  matrix[N,Na]  y_pred;
  matrix[N,Na]  yb_pred;
  matrix[N,Na]  yy_pred;
  for(i in 1:N)
    for(j in 1:Na) {
        y_pred[i,j] = poisson_log_rng(fyx[i,j]); 
       yb_pred[i,j] = poisson_log_rng(fyz[i,j]); 
       yy_pred[i,j] = poisson_log_rng(fy [i,j]); 
    }
}
