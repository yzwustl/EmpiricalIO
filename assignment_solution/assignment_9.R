library(doParallel)
library(dplyr)
library(reshape2)
library(Matrix)
library(ggplot2)
options(scipen = 200)
########## Part I: Simulate the data ##########
# set seed
set.seed(1)
# number of auctions
T <- 100
# parameter of value distribution
alpha <- 2
beta <- 2
# parameters of number of potential bidders
lambda <- 10
# number of bidders
N <- rpois(T, lambda)
N <- ifelse(N == 1, 2, N)
# draw valuations
valuation <- foreach(tt = 1:T, .combine = "rbind") %do% {
    n_t <- N[tt]
    header <- expand.grid(t = tt,i = 1:n_t) 
    return(header)
  }
valuation <- valuation %>% dplyr::mutate(x = rbeta(length(i),alpha,beta))
# draw reserve prices
reserve <- 0.2
reserve <- data.frame(t = 1:T, r = reserve)
# winning bids for second price auction
compute_winning_bids_second<-function(valuation, reserve){
  winning_bids <- valuation %>% 
    left_join(reserve, by = "t") %>% 
    group_by(t) %>%
    arrange(desc(x)) %>%
    summarise(n = max(i),
              m = length(x[x>=r]),
              maximum_bids = max(x),
              r = mean(r),
              w = nth(sort(x, decreasing = TRUE), 2)) %>%
    arrange(t) %>%
    mutate(w = ifelse(maximum_bids < r, "NA", w)) %>%
    select(-maximum_bids)
  return(winning_bids)
}
df_second_w <- compute_winning_bids_second(valuation = valuation, reserve = reserve)
 
# bids for first price auction
bid_first <- function(x, r, alpha, beta, n) {
  F_x <- function(y) pbeta(y, alpha, beta)^(n-1)
  if (x < r){
    return(0)
  }
  else{
    numerator <- integrate(F_x, r, x)$value
    denominator <- F_x(x)
    bid <- x - (numerator/denominator)
    return(bid)
  }
}
compute_bids_first<-function(valuation, reserve, alpha, beta){
  first_bids <- valuation %>% 
    left_join(reserve, by = "t") %>%
    group_by(t) %>%
    mutate(n = max(i),
           m = length(x[x>=r]),
           r = mean(r)) %>%
    ungroup() %>%
    rowwise() %>%
    mutate(b = bid_first(x, r, alpha, beta,n))
  return(first_bids)
}
df_first <- compute_bids_first(valuation, reserve, alpha, beta)

# winning bids for first price auction
compute_winning_bids_first<-function(valuation, reserve, alpha, beta){
  df_first <- compute_bids_first(valuation, reserve, alpha, beta)
  winning_bids <- df_first %>%
    group_by(t) %>%
    summarise(n = max(i),
              m = length(x[x>=r]),
              r = mean(r),
              w = max(b)) %>%
    arrange(t) %>%
    mutate(winning_bids = ifelse(w < r, "NA", w))
  return(winning_bids)
}
df_first_w <-compute_winning_bids_first(valuation, reserve, alpha, beta)
 
########## Part II: Parameter Estimation ##########
# compute probability density for winning bids from a second-price auction
w <- df_second_w[1, ]$w
r <- df_second_w[1, ]$r
m <- df_second_w[1, ]$m
n <- df_second_w[1, ]$n
compute_p_second_w<-function(w, r, m, n, alpha, beta){
  F_x <- function(y) pbeta(y, alpha, beta)
  f_x <- function(y) dbeta(y, alpha, beta)
  if (m==1) {
    p <- n * F_x(r)^(n-1) * (1 - F_x(r)) 
  }
  else {
    p <- n * (n-1) * F_x(w)^(n-2) * (1 - F_x(w)) * f_x(w)
  }
  return(p)
}
# compute non-participation probability
compute_m0 <- function(r, n, alpha, beta){
  F_x <- function(y) pbeta(y, alpha, beta)
  p <- F_x(r)^(n)
  return(p)
}
# compute log-likelihood for winning bids from second-price auctions
compute_loglikelihood_second_price_w<-function(theta, df_second_w){
  alpha <- theta[1]
  beta <- theta[2]
  output <- df_second_w %>%
    rowwise() %>%
    mutate(p = compute_p_second_w(w, r, m, n, alpha, beta),
           m0 = compute_m0(r, n, alpha, beta)) %>%
    mutate(loglikelihood = log(p/(1-m0)))
  likelihood <- mean(output$loglikelihood)
  return(likelihood)
}
theta <- c(alpha, beta)
result_second_parametric <-
  optim(
    par = theta,
    fn = compute_loglikelihood_second_price_w,
    df_second_w = df_second_w,
    method = "L-BFGS-B",
    control = list(fnscale = -1)
  )
comparison <- data.frame(true = theta,estimate = result_second_parametric$par)

# inverse bid function for first-price auction
r <- df_first_w[1, "r"] %>%
  as.numeric()
n <- df_first_w[1, "n"] %>%
  as.integer()
b <- 0.5 * r + 0.5 
x <- 0.5

# compute inverse bid equation
inverse_bid_equation<-function(x, b, r, alpha, beta, n){
  beta_x <- bid_first(x,r,alpha,beta,n)
  dif    <- beta_x - b
  return(dif)
}

# compute inverse bid function
inverse_bid_first<-function(b, r, alpha, beta, n){
  x <- uniroot(f = inverse_bid_equation, 
               lower = r,
               upper = 1,
               alpha = alpha,
               beta = beta, 
               r = r,
               n = n, 
               b = b)
  x <- x$root
  return(x)
}
 
# compute probability density for a winning bid from a first-price auction
compute_p_first_w<-function(w, r, alpha, beta, n){
  if (w > bid_first(1,r,alpha,beta,n)){
    return(1e-6)
  }
  F_x <- function(y) pbeta(y, alpha, beta) 
  eta <- function(y) inverse_bid_first(y,r,alpha,beta,n)
  numerator <- n * F_x(eta(w))^n
  denominator <- (n-1) * (eta(w) - w)
  h <- numerator/denominator
  return(h)
}
w <- 0.5
compute_p_first_w(
  w = w, 
  r = r, 
  alpha = alpha, 
  beta = beta, 
  n = n
)
upper <- 
  bid_first(
    x = 1, 
    r = r, 
    alpha = alpha, 
    beta = beta, 
    n = n
  )
# compute log-likelihood for winning bids for first-price auctions
compute_loglikelihood_first_price_w<-function(theta, df_first_w){
  alpha <- theta[1]
  beta <- theta[2]
  F_x <- function(y) pbeta(y, alpha, beta)
  output <- df_first_w %>%
    rowwise() %>%
    mutate(h = compute_p_first_w(w, r, alpha, beta, n),
           d = 1 - F_x(r)^n) %>%
    mutate(loglikelihood = log(h/d))
  likelihood <- mean(output$loglikelihood)
  return(likelihood)
}
compute_loglikelihood_first_price_w(
  theta = theta, 
  df_first_w = df_first_w
) 
result_first_parametric <-
  optim(
    par = theta,
    fn = compute_loglikelihood_first_price_w,
    df_first_w = df_first_w,
    method = "Nelder-Mead",
    control = list(fnscale = -1)
  )
comparison <-
  data.frame(
    true = theta,
    estimate = result_first_parametric$par
  )
comparison

# Non-parametric estimation
F_b <- ecdf(df_first$b)
f_b <- approxfun(density(df_first$b))

H_b <- function(b, n, F_b){
  return(F_b(b)^(n-1)) 
}

h_b <- function(b, n, F_b, f_b){
 return((n-1) * F_b(b)^(n-2) * f_b(b))
}

compute_implied_valuation<-function(b, n, r, F_b, f_b){
  H_b <- H_b(b, n, F_b)
  h_b <- h_b(b, n, F_b, f_b)
  if (b < r) {
    x <- 0
  } else {
    x <- b + H_b / h_b
  }
  return(x)
}
r <- df_first[1, "r"]
n <- df_first[1, "n"]
compute_implied_valuation(b = 0.4,n = n, r = r, F_b = F_b,f_b = f_b)

valuation_implied <- df_first %>% 
  rowwise() %>%
  mutate(x = compute_implied_valuation(b, n, r, F_b, f_b)) %>%
  ungroup() %>%
  select(x) %>%
  mutate(type = "estimate")

valuation_true <- valuation %>% dplyr::select(x) %>% mutate(type = "true")

valuation_plot <- rbind(valuation_true,valuation_implied)
ggplot(
  valuation_plot,
  aes(
    x = x, 
    color = type
  )
) + 
  stat_ecdf() + 
  theme_classic()



