functions {
  int neg_binomial_2_log_safe_rng(real eta, real phi) {
    real phi_div_exp_eta;
    real gamma_rate;
    phi_div_exp_eta = phi/exp(eta);
    gamma_rate = gamma_rng(phi, phi_div_exp_eta);
    if (gamma_rate >= exp(20.79))
      return -9;
    return poisson_rng(gamma_rate);
  }
}
data {
  int<lower=1> N;
  int qty_sld[N];
  vector[N] price;
}
parameters {
  real alpha;
  real beta;
  real<lower=0> phi;
}
model {
  beta ~ normal(0, 1);
  alpha ~ normal(0, 10);
  phi ~ normal(0, 1);
  
  qty_sld ~ neg_binomial_2_log(alpha + beta * price,
                               phi);
} 
generated quantities {
  vector[N] pp_y;
  
  for (n in 1:N) 
    pp_y[n] = neg_binomial_2_log_safe_rng(alpha 
                                          + beta * price[n],
                                          phi);
}
