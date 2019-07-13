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
  int<lower=1> M;
  int<lower=1, upper=J> ttl_idx[N];
  int<lower=1, upper=M> mo_idx[N];
}
parameters {
  real alpha;
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_beta;
  vector[J] alphas_raw;
  real beta;
  vector[J] betas_raw;
  real<lower=0> sigma_mo;
  vector[M] mos;
  real<lower=0> phi;
}
transformed parameters {
  vector[J] alphas;
  vector[J] betas;
  
  alphas = alpha + sigma_alpha * alphas_raw;
  betas = beta + sigma_beta * betas_raw;
}
model {
  beta ~ normal(0, 1);
  alphas_raw ~ normal(0, 1);
  betas_raw ~ normal(0, 1);
  sigma_alpha ~ normal(0, 1);
  sigma_beta ~ normal(0, 1);
  sigma_mo ~ normal(0, 1);
  alpha ~ normal(0, 1);
  mos ~ normal(0, sigma_mo);
  phi ~ normal(0, 1);
  
  qty_sld ~ neg_binomial_2_log(alphas[ttl_idx] + betas[ttl_idx] .* price + mos[mo_idx],
                               phi);
} 
generated quantities {
  vector[N] pp_y;
  
  for (n in 1:N) 
    pp_y[n] = neg_binomial_2_log_safe_rng(alphas[ttl_idx[n]]
                                          + betas[ttl_idx[n]] * price[n]
                                          + mos[mo_idx[n]], 
                                          phi);
  
}
