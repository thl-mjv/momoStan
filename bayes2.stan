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
  real<lower=0> length_scale;
  real<lower=0> alpha;
  vector[N_tot] eta;
}
transformed parameters {
  vector[N_tot] f;
  {
    matrix[N_tot, N_tot] L;
    matrix[N_tot, N_tot] K;
    K = cov_exp_quad(x_tot, alpha, length_scale);
    for (n in 1:N_tot)
      K[n, n] = K[n, n] + 1e-12;
    L = cholesky_decompose(K);
    f = L * eta;
  }
}
model {
  length_scale ~ gamma(2, 20);
  alpha ~ normal(0, 1);
  eta ~ normal(0, 1);
  y ~ poisson_log(f[1:N]);
}
generated quantities {
  int<lower=0> y_pred[N_tot];
  for(i in 1:N_tot)
    y_pred[i] = poisson_log_rng(f[i]);
}
