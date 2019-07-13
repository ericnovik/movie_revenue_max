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
  int<lower=1> J;
  int<lower=1, upper=J> ttl_idx[N];
}
parameters {
  real alpha;
  real<lower=0> sigma_alpha;
  vector[J] alphas;
  real beta;
  real<lower=0> phi;
}
model {
  beta ~ normal(0, 1);
  alphas ~ normal(alpha, sigma_alpha);
  sigma_alpha ~ normal(0, 1);
  alpha ~ normal(0, 1);
  phi ~ normal(0, 1);
  
  qty_sld ~ neg_binomial_2_log(alphas[ttl_idx] + beta * price,
                               phi);
} 
generated quantities {
  vector[N] pp_y;
  
  for (n in 1:N) 
    pp_y[n] = neg_binomial_2_log_safe_rng(alphas[ttl_idx[n]] + beta * price[n],
                                          phi);
}
