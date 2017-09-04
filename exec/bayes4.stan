// regression version?

data {
  int<lower=1> Nt; // observations per time series
  int<lower=1> Na; // number of age- and region groups
  int<lower=1> P ; // covariates excluding intercept
  int<lower=1> Np; // number of predictions
  int<lower=1> Nx; // number of 2nd level covariates
  int y  [Nt,Na];
  real z  [Nt,Na];
  real pop[Nt,Na]; // cause we model this as normal
  matrix[Nt,P] x;
  matrix[Np,P] xp;
  matrix[Na,Nx] xx;
}
transformed data {
  int<lower=1> N;
  int<lower=1> k;
  matrix[Nt + Np,P] xt;
  N = Nt + Np;
  k = 1;
  for (n in 1:Nt) {
    for(i in 1:P) {
      xt[k,i] = x[n,i];
    }
    k = k + 1;
  }
  for (n in 1:Np) {
    for(i in 1:P) {
      xt[k,i] = xp[n,i];
    }
    k = k + 1;
  }
}
parameters {
  // deaths
  matrix[P,Na] alphay;     // first levels
  real deltay[Na];   // effect of temp
  real<lower=0> sigmay[P+1]; // shrinkage
  matrix [Nx,P+1] betay;    // hierarchical

  // temp
  matrix[P,Na] alphat;     // first levels
  vector<lower=0>[Na] gammat; // temp residual vriance
  real<lower=0> sigmat[P+1]; // shrinkage
  matrix [Nx,P+1] betat;    // hierarchical

  // population
  matrix[P,Na] alphap;     // first levels
  vector<lower=0>[Na] gammap; // pop residual vriance
  real<lower=0> sigmap[P+1]; // shrinkage
  matrix [Nx,P+1] betap;    // hierarchical


}
transformed parameters {
  matrix<lower=-10,upper=10>[N,Na] fy;
  matrix                    [N,Na] ft;
  matrix                    [N,Na] fp;

  for(i in 1:Na) {
    fy[,i] = xt * alphay[,i]; // offset and the effect of the temperature ; //  + ft[,i]*deltay[i] + fp[,i] //  + z[,i]*deltay[i] + pop[,i]
    ft[,i] = xt * alphat[,i];
    fp[,i] = xt * alphap[,i];
  }
  
}
model {

  for(i in 1:P) {
    sigmay[ i] ~ normal(0.0,1.0);
    betay [,i] ~ normal(0.0,1.0);
    alphay[i,] ~ normal(xx*betay[,i], sigmay[i]);  

    sigmat[ i] ~ normal(0.0,1.0);
    betat [,i] ~ normal(0.0,1.0);
    alphat[i,] ~ normal(xx*betat[,i], sigmat[i]);  

    sigmap[ i] ~ normal(0.0,1.0);
    betap [,i] ~ normal(0.0,1.0);
    alphap[i,] ~ normal(xx*betap[,i], sigmap[i]);  
  }
  
  for(i in 1:Na) {
    deltay[i] ~ normal(xx*betay[,P+1],sigmay[P+1]);
    gammat[i] ~ normal(xx*betat[,P+1],sigmat[P+1]);
    gammap[i] ~ normal(xx*betap[,P+1],sigmap[P+1]);
    y  [,i] ~ poisson_log(fy[1:Nt,i]); // +pop[1:Nt,i]); # z[1:Nt,i]*deltay[i]+
    z  [,i] ~ normal     (ft[1:Nt,i],gammat[i]);
    pop[,i] ~ normal     (fp[1:Nt,i],gammap[i]);
  }
}
generated quantities {
  matrix[N,Na] y_pred;
  matrix[N,Na] z_pred;
  matrix[N,Na] p_pred;
  for(i in 1:N)
    for(j in 1:Na) {
      y_pred[i,j] = poisson_log_rng(fy[i,j]); // +z_pred[i,j]*deltay[j]+p_pred[i,j]
      z_pred[i,j] =      normal_rng(ft[i,j],gammat[j]);
      p_pred[i,j] =      normal_rng(fp[i,j],gammap[j]);
    }
}
