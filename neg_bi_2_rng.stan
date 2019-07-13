functions {
  int neg_bi_rng(real log_mu, real sig) {
    return neg_binomial_2_log_rng(log_mu, sig);
  }
}
data { }
model { }
