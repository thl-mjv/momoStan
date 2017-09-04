// AR version?

data {
  int<lower=1> N;
  int<lower=1> N_pred;
  int y[N];
  real x[N];
  real x_pred[N_pred];
}
transformed data {
  int<lower=1> N_tot;
  int<lower=1> k;
  real x_tot[N_pred + N];
  N_tot = N_pred + N;
  k = 1;
  for (n in 1:N) {
    x_tot[k] = x[n];
    k = k + 1;
  }
  for (n in 1:N_pred) {
    x_tot[k] = x_pred[n];
    k = k + 1;
  }
}
parameters {
  real<lower=1,upper=10> alpha0; // intercept
  real alpha1; // AR-1 <lower=-1,upper=1>
  real alpha2; // AR-2 <lower=-1,upper=1>
  real eta[N_tot] ; // innovations
}
transformed parameters {
  real<lower=-19,upper=19>  f[N_tot];
  f[1] = eta[1]/sqrt(1-alpha1*alpha1-alpha2*alpha2)+alpha0;
  f[2] = eta[2]*sqrt(2*(1-alpha1*alpha2))/sqrt(1-alpha1*alpha1-alpha2*alpha2)+f[1];
  for(i in 3:N_tot)
    f[i]=alpha1*(f[i-1]-alpha0)+alpha2*(f[i-2]-alpha0)+eta[i]+alpha0;
  
}
model {
  alpha0 ~ normal(7, .1);
  alpha1 ~ normal(0, .1);
  alpha2 ~ normal(0, .1);
  eta    ~ normal(0, .1);
  y      ~ poisson_log(f[1:N]);
}
generated quantities {
  int<lower=0> y_pred[N_tot];
  for(i in 1:N_tot)
    y_pred[i] = poisson_log_rng(f[i]);
}
