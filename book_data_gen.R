library(lubridate)
library(rstan)
library(dplyr)
library(zoo)

cov_exp_quad <- function(x, alpha, len_scale, nug) {
  dim_x = length(x)
  cov_ret <- matrix(NA_real_, dim_x, dim_x)
  for (i in 1:(dim_x - 1)) {
    cov_ret[i, i] = alpha ^ 2 + nug
    for (j in (i + 1):dim_x) {
      r_x = (x[i] - x[j])^2
      cov_ret[i, j] = alpha^2 * exp(-1 / (2 * len_scale^2) * r_x)
      cov_ret[j,i] = cov_ret[i, j]
    }
  }
  cov_ret[dim_x, dim_x] = alpha ^ 2 + nug
  return(cov_ret)
}

tpois <- function(mu, lower, upper) {
  samp <- rpois(1, mu)
  while (samp < lower | samp > upper)
    samp <- rpois(1, mu)
  return(samp)
}

rtpois <- function(N, mu, lower, upper) {
  samps <- rep(NA_real_)
  for (n in 1:N)
    samps[n] <- tpois(mu, lower, upper)
  return(samps)
}

gen_data <- function(N_books, N_wks) {
  set.seed(123)
  book_inds <- as.vector(sapply(1:N_books, rep, N_wks))
  wk_inds <- rep(1:N_wks, N_books)
  dates <- as_date('2015-01-01') + (1:N_wks)*7
  n_price_changes <- rtpois(N_books,3,1,5)
  starting_prices <- rtpois(N_books,9,6,15)
  prices <- sapply(n_price_changes,function(x) sort(sample(2:N_wks,x, replace=F)))
  book_prices <- c()
  len_books <- 5
  sig_books <- 0.5
  noise_books <- 0.15
  book_eta <- matrix(rt(N_wks * N_books, df = 4),N_wks,N_books)
  K <- cov_exp_quad(x = 1:N_wks, alpha = sig_books, len_scale = len_books, nug = noise_books)
  L_K <- t(chol(K))
  book_series <- L_K %*% book_eta
  book_series <- as.vector(book_series)
  for (n in 1:N_books) {
    changes <- prices[[n]]
    starting_price <- starting_prices[n]
    level_changes <- 2*rbinom(length(changes), 1, 0.5) - 1
    levels <- cumsum(c(starting_price, level_changes))
    df_levels <- data.frame(price = levels, ind = c(1,changes))
    p_series <- rep(NA_real_,N_wks)
    p_series[df_levels$ind] <- df_levels$price
    p_series <- zoo::na.locf(p_series)
    book_prices <- c(book_prices, p_series)
  }
  overall_mu <- log(30)
  sd_alpha <- 1/10
  overall_beta <- -0.2
  sd_beta <- 0.05
  alphas <- overall_mu + sd_alpha * rnorm(N_books)
  betas <- overall_beta + sd_beta * rnorm(N_books)
  wk_noise <- rnorm(N_wks) * 0.25
  mus <- alphas[book_inds] + betas[book_inds] * book_prices + wk_noise[wk_inds] + book_series
  data <- rpois(length(mus), exp(mus))
  df_ret <- data.frame(y = data, prices = book_prices, wk_ind = wk_inds, dates = dates[wk_inds],
                       book = book_inds)
  return(df_ret)
}

df_test <- gen_data(100, 50)
with(filter(df_test, book == 5), plot(dates, y, type='l'))
