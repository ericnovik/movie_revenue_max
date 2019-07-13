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
  int<lower=1> W;
  int<lower=1, upper=J> ttl_idx[N];
  int<lower=1, upper=M> mo_idx[N];
  int<lower=1,upper=W> wk_idx[N];
  vector[J] rel_year;
  vector[N] xmas_ind;
  matrix[J,2] meta_data;
}
transformed data {
  vector[20] hypo_prices;
  matrix[J,2] scaled_meta;
  int max_mo;
  
  for (i in 1:cols(meta_data)) 
    scaled_meta[,i] = (meta_data[,i] - mean(meta_data[,i])) / sd(meta_data[,i]);
  
  for (i in 1:20)
    hypo_prices[i] = i;
  max_mo = max(mo_idx);
}
parameters {
  real alpha;
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_beta;
  vector[J] alphas_raw;
  real beta;
  vector[2] gamma_meta;
  vector[J] betas_raw;
  real<lower=0> sigma_mo;
  real<lower=0> sigma_wk;
  vector[W] wk_raw;
  vector[M] mos_raw;
  real gamma_rel;
  real<lower=0> phi;
}
transformed parameters {
  vector[J] alphas;
  vector[J] betas;
  vector[W] wk;
  vector[M] mos;
  
  mos = sigma_mo * mos_raw;
  wk = sigma_wk * wk_raw;
  betas = beta + scaled_meta * gamma_meta + sigma_beta * betas_raw;
  alphas = alpha + gamma_rel * rel_year + sigma_alpha * alphas_raw;
}
model {
  beta ~ normal(0, 1);
  alphas_raw ~ normal(0, 1);
  betas_raw ~ normal(0, 1);
  sigma_alpha ~ normal(0, 1);
  gamma_meta ~ normal(0, 1);
  sigma_beta ~ normal(0, 1);
  sigma_mo ~ normal(0, 1);
  sigma_wk ~ normal(0, 1);
  gamma_rel ~ normal(0, 1);
  alpha ~ normal(0, 1);
//  mos ~ normal(0, sigma_mo);
  mos_raw ~ normal(0, 1);
  wk_raw ~ normal(0, 1);
  phi ~ normal(0, 1);
  
  qty_sld ~ neg_binomial_2_log(alphas[ttl_idx] + betas[ttl_idx] .* price
                               + mos[mo_idx] + wk[wk_idx],
                               phi);
} 
generated quantities {
  vector[N] pp_y;
  real new_wk;
  matrix[J, rows(hypo_prices)] hypo_rev;
  matrix[J, rows(hypo_prices)] hypo_qty_sld;
  new_wk = normal_rng(0, sigma_wk);
  
  for (n in 1:N) 
    pp_y[n] = neg_binomial_2_log_safe_rng(alphas[ttl_idx[n]]
                                          + betas[ttl_idx[n]] * price[n]
                                          + mos[mo_idx[n]]
                                          + wk[wk_idx[n]],
                                          phi);
  for (j in 1:J) {
    for (p in 1:rows(hypo_prices)) {
      hypo_qty_sld[j,p] = neg_binomial_2_log_safe_rng(alphas[j]
                                                      + betas[j] * hypo_prices[p]
                                                      + mos[max_mo]
                                                      + new_wk,
                                                      phi); 
      hypo_rev[j,p] =  hypo_qty_sld[j,p] * hypo_prices[p];
    }
  }
}
