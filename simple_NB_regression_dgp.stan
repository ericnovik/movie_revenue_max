data {
  int<lower=1> N;
  int qty_sld[N];
  vector[N] price;
  int<lower=1> J;
  int<lower=1> t;
  int<lower=1, upper=J> ttl_idx[N];
  real alpha;
  real beta;
  real phi;
}
model {
} 
generated quantities {
  vector[N] y_gen;
  
  for (n in 1:N) 
    y_gen[n] = neg_binomial_2_log_rng(alpha + beta * price[n], phi);
}
