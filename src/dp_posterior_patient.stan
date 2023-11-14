// posterior sampling of the parameter where each patient is treated as one observation
functions {
  // when nrows per patient==1
  real wphr1_lpdf(real t, row_vector x, real d, real lambda, 
                 real beta0, vector beta, real eta, real epsilon){
    real y;
    real a;
    y = (log(lambda) + beta0 + dot_product(beta, x) 
    + (lambda - 1) * log(t))*d // hazard 
    - exp(beta0 + dot_product(beta, x)) * t^lambda; // survival
    a = epsilon / (2 * eta);
    if (y <= 0 && y >= -eta) 
      return a * y;
    if (y > 0) 
      return 0;
    else
      return -a * eta;
  }
  
  real wphr_lpdf(vector t, int m, int N, matrix x, vector d, real lambda, 
                 real beta0, vector beta, real eta, real epsilon){
    real y;
    real a;
    y=0;
    for (i in 1:m){
      
      y += (log(lambda) + beta0 + dot_product(beta, x[i]) 
      + (lambda - 1) * log(t[i]))*d[i] // hazard 
      - exp(beta0 + dot_product(beta, x[i])) * t[i]^lambda; // survival
    }
    y = y/m; // divide on the number of transition per patient
    a = epsilon / (2 * eta);
    if (y <= 0 && y >= -eta) 
      return a * y;
    if (y > 0) 
      return 0;
    else
      return -a * eta;
  }
}

data {
  int<lower=0> M; // total number of patients
  int<lower=0> N; // number of rows
  int m[M];       // number of rows for each patient 
  matrix[N, 3] X; // array with vectors of covariates for each observation
  vector[N] t; // vector with trans times for each observation
  vector[N] d; // vector with boolean 0=right-censored 1=observed
  real<lower=0> eta; // 
  real<lower=0> epsilon; // epsilon dp
}

parameters {
  real<lower=0> lambda;  // shape parameter
  real beta0;   // intercept
  vector[3] beta; // age, sex, and treatment group coefficients
}

model {
  // prior
  lambda ~ normal(0, 10);
  beta0 ~ normal(0, 10);
  beta ~ normal(0, 10);      // vectorized
  
  // likelihood
  int count;
  count = 1;
  for (i in 1:M){
// trick to use ragged data structures
      if (m[i]==1){
        target += wphr1_lpdf(t[count] |X[count,1:3], 
          d[count], lambda, beta0, beta, eta, epsilon);
          count += m[i];
      }
      else{
        int start_idx = count;
        int end_idx = count + m[i] - 1;
        target += wphr_lpdf(t[start_idx:end_idx] | m[i], N, X[start_idx:end_idx, 1:3], 
          d[start_idx:end_idx], lambda, beta0, beta, eta, epsilon);
        count += m[i];
      }
  }
}

