library(dplyr)
library(evd)
library(doParallel)
options(scipen = 200)
#### Generate Data ####
set.seed(1)
J <- 10  # number of products
K <- 3   # dimension of product characteristics including the intercept
T <- 100 # number of markets
N <- 500 # number of consumers per market
L <- 500 # number of Monte Carlo

# set parameters of interests
beta <- rnorm(K); 
beta[1] <- 4

sigma <- abs(rnorm(K))

mu <- 0.5
omega <- 1

price_xi <- 1
prop_jt <- 0.6
sd_x <- 0.5
sd_c <- 0.05
sd_p <- 0.05

# Generate product-market characteristics
j  <- 1:10
x_1 <- rep(1, J)
x_2 <- rnorm(J,0,sd_x)
x_3 <- rnorm(J,0,sd_x)
x  <- data.frame(j,x1,x2,x3)
x <- rbind(c(0,0,0,0),x)

# Generate market-product data
M <- expand.grid(j=1:J,t=1:T)
M$xi <- 0
M$c <- exp(rnorm(J*T,0,sd_c))
M$p <- M$c + exp(rnorm(J*T,price_xi * M$xi,sd_p))
M <- M %>% group_by(t) %>% sample_frac(prop_jt)

# add rows for each  with other columns to be 0
outside <- data.frame(
    j = 0, 
    t = 1:T, 
    xi = 0, 
    c = 0, 
    p = 0
)
M <- rbind(M,outside) %>% dplyr::arrange(t, j)

# generate consumer-market data, idiocyncratic shocks
V <- expand.grid(i=1:N,t=1:T)
V$v_x_1 <- rnorm(N*T,0,1) # beta_1_variation
V$v_x_2 <- rnorm(N*T,0,1) # beta_2_variation
V$v_x_3 <- rnorm(N*T,0,1) # beta_3_variation
V$v_p <- rnorm(N*T,0,1)   # alpha variation

# join M,V,x
df <- expand.grid(t = 1:T, i = 1:N, j = 0:J)
df <- left_join(df, M, by = c("j", "t"))
df <- left_join(df, V, by = c("i", "t"))
df <- left_join(df, x, by = "j")

df <- df %>% dplyr::filter(!is.na(p)) %>% dplyr::arrange(t,i,j)

# # draw idiosyncratic shocks
e <- evd::rgev(dim(df)[1])
df$e <- e

#### Utility Function ####
compute_indirect_utility <- function(df, beta, sigma, mu, omega){
  u = (beta[1]+sigma[1]*df$v_x_1) * df$x1 + 
      (beta[2]+sigma[2]*df$v_x_2) * df$x2 + 
      (beta[3]+sigma[3]*df$v_x_3) * df$x3 - 
      exp(mu+omega*df$v_p) * df$p  + df$xi + df$e
  return(u)
}
u <- 
  compute_indirect_utility(
    df = df, 
    beta = beta, 
    sigma = sigma, 
    mu = mu, 
    omega = omega
  )
head(u)
#### Choice Function ####
compute_choice <- function(X, M, V, e, beta, sigma, mu, omega){
  # Part I: Construct df
  df <- expand.grid(t = 1:T, i = 1:N, j = 0:J)
  df <- left_join(df, M, by = c("j", "t"))
  df <- left_join(df, V, by = c("i", "t"))
  df <- left_join(df, x, by = "j")
  df <- df %>% dplyr::filter(!is.na(p)) %>% dplyr::arrange(t,i,j)
  df$e <- e
  
  # Part II: Compute indirect utility
  df$u <- compute_indirect_utility(df, beta, sigma, mu, omega)
  
  # Part III: Derive Choice, not probability
  df <- df %>% group_by(i,t) %>% mutate(q = (u==max(u))) %>% ungroup()
  df$q <- as.numeric(df$q)
  df_choice <- df
  return(df_choice)
}
df_choice <- compute_choice(df, M, V, e, beta, sigma, mu, omega)
summary(df_choice)

#### Share Function ####
compute_share <- function(X, M, V, e, beta, sigma, mu, omega){
  # Part I: Construct df
  df <- expand.grid(t = 1:T, i = 1:N, j = 0:J)
  df <- left_join(df, M, by = c("j", "t"))
  df <- left_join(df, V, by = c("i", "t"))
  df <- left_join(df, x, by = "j")
  df <- df %>% dplyr::filter(!is.na(p)) %>% dplyr::arrange(t,i,j)
  df$e <- e
  
  # Part II: Compute Choice
  df_choice <- compute_choice(X, M, V, e, beta, sigma, mu, omega)
  
  # Part III: Compute Share and log difference with j=0
  df_share <- df_choice %>% group_by(t, j) %>% summarise(q = sum(q), share = sum(q)/N, .groups = "drop")
  
  # Part IV: Compute log difference
  df_share <- df_share %>% 
    group_by(t) %>% 
    mutate(y = log(share / share[j == 0])) %>% 
    ungroup()
  
  # Merge with M
  df_share <- left_join(df_share, M,by = c("j", "t"))
  
  # Merge with x
  df_share <- left_join(df_share, x, by = "j")
  
  return(df_share)
}
df_share <- compute_share(df, M, V, e, beta, sigma, mu, omega)

#### Parameter Estimation ####

# Model 1: OLS
ols <- lm(y ~ x1 + x2 + x3 + p, data = df_share)
summary(ols)

# Simulate Based Estimator
V_mcmc <- expand.grid(i=1:N,t=1:T)
V_mcmc$v_x_1 <- rnorm(N*T,0,1) # beta_1_variation
V_mcmc$v_x_2 <- rnorm(N*T,0,1) # beta_2_variation
V_mcmc$v_x_3 <- rnorm(N*T,0,1) # beta_3_variation
V_mcmc$v_p <- rnorm(N*T,0,1)   # alpha variation

e_mcmc <- evd::rgev(dim(df)[1])
 
df_mcmc <- expand.grid(t = 1:T, i = 1:L, j = 0:J) 
df_mcmc <- left_join(df_mcmc, V_mcmc, by = c("i", "t"))
df_mcmc <- left_join(df_mcmc, M, by = c("j", "t"))
df_mcmc <- left_join(df_mcmc, x, by = "j")
df_mcmc <- df_mcmc %>% dplyr::filter(!is.na(p)) %>% dplyr::arrange(t,i,j)
df_mcmc$e <- e_mcmc

# vectorize parameters
theta <- c(beta, sigma, mu, omega)

# objective_function
compute_nlls_objective_a3 <- function(theta, df_share, X, M, V_mcmc, e_mcmc){
  beta <- theta[1:3]
  sigma <- theta[4:6]
  mu <- theta[7]
  omega <- theta[8]
  
  df_share_hat <- compute_share(X, M, V_mcmc, e_mcmc, beta, sigma, mu, omega)
  
  nlls <- mean((df_share_hat$share - df_share$share)^2)
  return(nlls)
}

registerDoParallel()

# Optimiazation
result_NLLS <- optim(par = theta, fn = compute_nlls_objective_a3,
                     method = "Nelder-Mead",
                     df_share = df_share, 
                     X = X, 
                     M = M, 
                     V_mcmc = V_mcmc, 
                     e_mcmc = e_mcmc)
result <- data.frame(true = theta, estimates = result_NLLS$par)
result
