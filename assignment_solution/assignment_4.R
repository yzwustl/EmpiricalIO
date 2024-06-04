library(dplyr)
library(evd)
library(purrr)
library(doParallel)
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
sd_x <- 2
sd_xi <- 0.1
sd_c <- 0.5
sd_p <- 0.01

# Generate product-market characteristics
j  <- 1:10
x_1 <- rep(1, J)
x_2 <- rnorm(J,0,sd_x)
x_3 <- rnorm(J,0,sd_x)
x  <- data.frame(j,x_1,x_2,x_3)
x <- rbind(c(0,0,0,0),x)

# Generate market-product data
M <- expand.grid(j=1:J,t=1:T)
M$xi <- rnorm(J*T,0,sd_xi)
M$c <- exp(rnorm(J*T,0,sd_c))
M$p <- M$c + exp(rnorm(J*T,price_xi * M$xi,sd_p))
M <- M %>% group_by(t) %>% sample_frac(size = purrr::rdunif(1, J) / J) %>% dplyr::ungroup()

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

e <- evd::rgev(dim(df)[1])
df$e <- e

# Utility Function 
compute_indirect_utility <- function(df, beta, sigma, mu, omega){
  u = (beta[1]+sigma[1]*df$v_x_1) * df$x_1 + 
    (beta[2]+sigma[2]*df$v_x_2) * df$x_2 + 
    (beta[3]+sigma[3]*df$v_x_3) * df$x_3 - 
    (exp(mu+omega*df$v_p) * df$p) + df$xi
  return(u)
}
u <- compute_indirect_utility(df, beta, sigma, mu, omega)

# Smooth Choice Function
compute_choice_smooth <- function(df, beta, sigma, mu, omega){
  # Part I: Compute indirect utility
  df$u  <- compute_indirect_utility(df, beta, sigma, mu, omega)
  df$eu <- exp(df$u)
  
  # Part II: Compute share, group by market and consumer
  df <- df %>% group_by(t,i) %>% mutate(q = eu / sum(eu))
  
  return(df)
}

# Smooth Share Function 
compute_share_smooth <- function(df, beta, sigma, mu, omega){
  # Part I: Compute Choice
  df <- compute_choice_smooth(df, beta, sigma, mu, omega)
  
  # Part II: Compute Share and log difference with j=0
  df_share <- df %>% group_by(t, j) %>% summarise(q = sum(q), s = sum(q)/N, .groups = "drop")
  
  # Part III: Compute log difference
  df_share <- df_share %>% 
    group_by(t) %>% 
    mutate(y = log(s/s[j == 0])) %>% 
    ungroup()
  
  # Merge with M
  df_share <- left_join(df_share, M,by = c("j", "t"))
  
  # Merge with x
  df_share <- left_join(df_share, x, by = "j")
  
  return(df_share)
}
df_share_smooth <- compute_share_smooth(
    df,
    beta = beta, 
    sigma = sigma, 
    mu = mu, 
    omega = omega
  )
summary(df_share_smooth)
#### Estimation of Parameters ####
V_mcmc <- expand.grid(i=1:L,t=1:T)
V_mcmc$v_x_1 <- rnorm(L*T,0,1) # beta_1_variation
V_mcmc$v_x_2 <- rnorm(L*T,0,1) # beta_2_variation
V_mcmc$v_x_3 <- rnorm(L*T,0,1) # beta_3_variation
V_mcmc$v_p <- rnorm(L*T,0,1)   # alpha variation
head(V_mcmc)

df_mcmc <- expand.grid(t = 1:T, i = 1:L, j = 0:J)
df_mcmc <- left_join(df_mcmc, M, by = c("j", "t"))
df_mcmc <- left_join(df_mcmc, V_mcmc, by = c("i", "t"))
df_mcmc <- left_join(df_mcmc, x, by = "j")
df_mcmc <- df_mcmc %>% dplyr::filter(!is.na(p)) %>% dplyr::arrange(t,i,j)

e_mcmc <- evd::rgev(dim(df_mcmc)[1])
head(e_mcmc)
df_mcmc$e <- e_mcmc
# set parameters
theta <- c(beta, sigma, mu, omega)

# Compute Delta
XX <- as.matrix(dplyr::select(df_share_smooth, dplyr::starts_with("x_")))
pp <- as.matrix(dplyr::select(df_share_smooth, p)) 
xi <- as.matrix(dplyr::select(df_share_smooth, xi))
alpha <- - exp(mu + omega^2/2)
delta <- XX %*% as.matrix(beta) + pp * alpha + xi
delta <- dplyr::select(df_share_smooth, t, j) %>%
  dplyr::mutate(delta = as.numeric(delta))

# Compute indirect utility
compute_indirect_utility_delta <- function(df, delta, sigma, mu, omega){
  df_all  <- left_join(df, delta, by = c("t", "j"))
  u_delta <- df_all$delta + 
             sigma[1] * df_all$v_x_1 * df_all$x_1 + 
             sigma[2] * df_all$v_x_2 * df_all$x_2 + 
             sigma[3] * df_all$v_x_3 * df_all$x_3 - 
             (exp(mu+omega*df_all$v_p) - exp(mu+omega^2/2)) * df_all$p
  return(u_delta)
}

u_delta <-
  compute_indirect_utility_delta(
    df, 
    delta, 
    sigma,
    mu, 
    omega
  )

# Smooth Choice Function
compute_choice_smooth_delta <- function(df, delta, sigma, mu, omega){
  # Part I: Compute indirect utility
  df$u  <- compute_indirect_utility_delta(df, delta, sigma, mu, omega)
  df$eu <- exp(df$u)
  
  # Part II: Compute share, group by market and consumer
  df <- df %>% group_by(t,i) %>% mutate(q = eu / sum(eu))
  
  return(df)
}
 
# share function
compute_share_smooth_delta <- function(df, delta, sigma, mu, omega){
  # Part I: Compute Choice
  df <- compute_choice_smooth_delta(df, delta, sigma, mu, omega)
  
  # Part II: Compute Share and log difference with j=0
  df_share <- df %>% group_by(t, j) %>% summarise(q = sum(q), s = sum(q)/N, .groups = "drop")
  
  # Part III: Compute log difference
  df_share <- df_share %>% 
    group_by(t) %>% 
    mutate(y = log(s/s[j == 0])) %>% 
    ungroup()
  
  # Merge with M
  df_share <- left_join(df_share, M,by = c("j", "t"))
  
  # Merge with x
  df_share <- left_join(df_share, x, by = "j")
  
  return(df_share)
}
df_share_smooth_delta <-
  compute_share_smooth_delta(
    df,
    delta, 
    sigma, 
    mu, 
    omega
  ) 
df_share_smooth_delta

# Find Delta:
solve_delta <- function(df_share_smooth, df, delta_init, sigma, mu, omega, kappa, lambda){
  diff <- 100  # Initialize diff to be larger than lambda
  delta_new <- delta_init  # Initialize delta to the initial values
  
  while (diff > lambda) {
    df_predicted <- compute_share_smooth_delta(df, delta_new, sigma, mu, omega)
    T_delta <- delta_new$delta + kappa * log(df_share_smooth$s / df_predicted$s)
    diff    <- max(abs(T_delta - delta_new$delta))
    delta_new$delta <- T_delta  # Update delta for the next iteration
  }
  return(delta_new)
}

kappa <- 1
lambda <- 1e-2

delta_new <- solve_delta(df_share_smooth, df, delta, sigma, mu, omega, kappa, lambda)
summary(delta_new$delta - delta$delta)

# Use MCMC shocks to estimate delta
delta_new2 <-
  solve_delta(
    df_share_smooth, 
    df_mcmc,
    delta, 
    sigma, 
    mu, 
    omega, 
    kappa, 
    lambda
  )
summary(delta_new2$delta - delta$delta)

# Linear Parameter Estimation
compute_theta_linear <- function(df_share_smooth, delta, mu, omega, Psi) {
  # Extract relevant variables from the data
  df_share_smooth <- df_share_smooth %>% 
    dplyr::filter(j != 0) %>% 
    dplyr::arrange(t, j)
  
  delta <- delta %>% 
    dplyr::filter(j != 0) %>% 
    dplyr::arrange(t, j)
  
  X <- as.matrix(df_share_smooth[, grep("^x_", names(df_share_smooth))])
  c <- as.matrix(df_share_smooth$c)
  p <- as.matrix(df_share_smooth$p)
  delta_v <- as.matrix(delta$delta)
  W = cbind(X,c)
  
  # Compute alpha0
  alpha0 <- -exp(mu + omega^2 / 2)
  
  # Compute beta0
  beta_0 <- solve(t(X) %*% W %*% solve(Psi) %*% t(W) %*% X) %*% t(X) %*% W %*% solve(Psi) %*% t(W) %*% (delta_v - alpha0*p)
  
  return(beta_0)
}
Psi <- diag(length(beta) + 1)
theta_linear <-
  compute_theta_linear(
    df_share_smooth, 
    delta, 
    mu, 
    omega, 
    Psi
  ) 
cbind(
  theta_linear, 
  beta
)

# Compute the fixed effects
solve_xi <- function(df_share_smooth, delta, beta, mu, omega){
  delta$eps = delta$delta - (as.matrix(df_share_smooth[, grep("^x_", names(df_share_smooth))]) %*% beta - exp(mu + omega^2 / 2) * as.matrix(df_share_smooth$p))
  delta = delta %>% dplyr::filter(j != 0) %>% dplyr::arrange(t, j)
  xi  = delta$eps
  return(xi)
}
xi_new <- 
  solve_xi(
    df_share_smooth, 
    delta, 
    beta, 
    mu, 
    omega
  )
head(xi_new)

# GMM Objective: BLP Algorithm
theta_nonlinear <- c(sigma,mu,omega)

compute_gmm_objective_a4 <- function(theta_nonlinear, delta, df_share_smooth, Psi, df_mcmc,kappa,lambda){
  #Compute delta
  delta_new_inf <- solve_delta(df_share_smooth, df_mcmc, delta, theta_nonlinear[1:3],theta_nonlinear[4],theta_nonlinear[5], kappa, lambda)
  
  #Compute theta_1
  theta_1 <- compute_theta_linear(df_share_smooth, delta_new_inf, theta_nonlinear[4],theta_nonlinear[5], Psi)
  
  #Compute xi
  xi <- solve_xi(df_share_smooth, delta_new_inf, theta_1, theta_nonlinear[4], theta_nonlinear[5]) 
  xi <- as.matrix(xi)
  
  # Get W
  df_share_smooth <- df_share_smooth %>% dplyr::filter(j != 0) %>% dplyr::arrange(t, j)
  
  X <- as.matrix(df_share_smooth[, grep("^x_", names(df_share_smooth))])
  c <- as.matrix(df_share_smooth$c)
  W = cbind(X,c)
  
  # Compute Obj
  obj = t(xi) %*% W %*% solve(Psi) %*% t(W) %*% xi
  
  return(obj)
}

objective <-
  compute_gmm_objective_a4(
    theta_nonlinear, 
    delta, 
    df_share_smooth, 
    Psi, 
    df_mcmc,
    kappa,
    lambda
  ) 

# Final Estimation
library(dplyr)
library(evd)
library(purrr)
library(doParallel)
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
sd_x <- 2
sd_xi <- 0.1
sd_c <- 0.5
sd_p <- 0.01

# Generate product-market characteristics
j  <- 1:10
x_1 <- rep(1, J)
x_2 <- rnorm(J,0,sd_x)
x_3 <- rnorm(J,0,sd_x)
x  <- data.frame(j,x_1,x_2,x_3)
x <- rbind(c(0,0,0,0),x)

# Generate market-product data
M <- expand.grid(j=1:J,t=1:T)
M$xi <- rnorm(J*T,0,sd_xi)
M$c <- exp(rnorm(J*T,0,sd_c))
M$p <- M$c + exp(rnorm(J*T,price_xi * M$xi,sd_p))
M <- M %>% group_by(t) %>% sample_frac(size = purrr::rdunif(1, J) / J) %>% dplyr::ungroup()

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

e <- evd::rgev(dim(df)[1])
df$e <- e

# Utility Function 
compute_indirect_utility <- function(df, beta, sigma, mu, omega){
  u = (beta[1]+sigma[1]*df$v_x_1) * df$x_1 + 
    (beta[2]+sigma[2]*df$v_x_2) * df$x_2 + 
    (beta[3]+sigma[3]*df$v_x_3) * df$x_3 - 
    (exp(mu+omega*df$v_p) * df$p) + df$xi
  return(u)
}
u <- compute_indirect_utility(df, beta, sigma, mu, omega)
# Smooth Choice Function
compute_choice_smooth <- function(df, beta, sigma, mu, omega){
  # Part I: Compute indirect utility
  df$u  <- compute_indirect_utility(df, beta, sigma, mu, omega)
  df$eu <- exp(df$u)
  
  # Part II: Compute share, group by market and consumer
  df <- df %>% group_by(t,i) %>% mutate(q = eu / sum(eu))
  
  return(df)
}

# share function 
compute_share_smooth <- function(df, beta, sigma, mu, omega){
  # Part I: Compute Choice
  df <- compute_choice_smooth(df, beta, sigma, mu, omega)
  
  # Part II: Compute Share and log difference with j=0
  df_share <- df %>% group_by(t, j) %>% summarise(q = sum(q), s = sum(q)/N, .groups = "drop")
  
  # Part III: Compute log difference
  df_share <- df_share %>% 
    group_by(t) %>% 
    mutate(y = log(s/s[j == 0])) %>% 
    ungroup()
  
  # Merge with M
  df_share <- left_join(df_share, M,by = c("j", "t"))
  
  # Merge with x
  df_share <- left_join(df_share, x, by = "j")
  
  return(df_share)
}
df_share_smooth <- compute_share_smooth(
  df,
  beta = beta, 
  sigma = sigma, 
  mu = mu, 
  omega = omega
)
summary(df_share_smooth)
#### Estimation of Parameters ####
V_mcmc <- expand.grid(i=1:L,t=1:T)
V_mcmc$v_x_1 <- rnorm(L*T,0,1) # beta_1_variation
V_mcmc$v_x_2 <- rnorm(L*T,0,1) # beta_2_variation
V_mcmc$v_x_3 <- rnorm(L*T,0,1) # beta_3_variation
V_mcmc$v_p <- rnorm(L*T,0,1)   # alpha variation
head(V_mcmc)

df_mcmc <- expand.grid(t = 1:T, i = 1:L, j = 0:J)
df_mcmc <- left_join(df_mcmc, M, by = c("j", "t"))
df_mcmc <- left_join(df_mcmc, V_mcmc, by = c("i", "t"))
df_mcmc <- left_join(df_mcmc, x, by = "j")
df_mcmc <- df_mcmc %>% dplyr::filter(!is.na(p)) %>% dplyr::arrange(t,i,j)

e_mcmc <- evd::rgev(dim(df_mcmc)[1])
head(e_mcmc)
df_mcmc$e <- e_mcmc
# set parameters
theta <- c(beta, sigma, mu, omega)

# Compute Delta
XX <- as.matrix(dplyr::select(df_share_smooth, dplyr::starts_with("x_")))
pp <- as.matrix(dplyr::select(df_share_smooth, p)) 
xi <- as.matrix(dplyr::select(df_share_smooth, xi))
alpha <- - exp(mu + omega^2/2)
delta <- XX %*% as.matrix(beta) + pp * alpha + xi
delta <- dplyr::select(df_share_smooth, t, j) %>%
  dplyr::mutate(delta = as.numeric(delta))

# Compute indirect utility
compute_indirect_utility_delta <- function(df, delta, sigma, mu, omega){
  df_all  <- left_join(df, delta, by = c("t", "j"))
  u_delta <- df_all$delta + 
    sigma[1] * df_all$v_x_1 * df_all$x_1 + 
    sigma[2] * df_all$v_x_2 * df_all$x_2 + 
    sigma[3] * df_all$v_x_3 * df_all$x_3 - 
    (exp(mu+omega*df_all$v_p) - exp(mu+omega^2/2)) * df_all$p
  return(u_delta)
}

u_delta <-
  compute_indirect_utility_delta(
    df, 
    delta, 
    sigma,
    mu, 
    omega
  )

# Smooth Choice Function
compute_choice_smooth_delta <- function(df, delta, sigma, mu, omega){
  # Part I: Compute indirect utility
  df$u  <- compute_indirect_utility_delta(df, delta, sigma, mu, omega)
  df$eu <- exp(df$u)
  
  # Part II: Compute share, group by market and consumer
  df <- df %>% group_by(t,i) %>% mutate(q = eu / sum(eu))
  
  return(df)
}

# share function
compute_share_smooth_delta <- function(df, delta, sigma, mu, omega){
  # Part I: Compute Choice
  df <- compute_choice_smooth_delta(df, delta, sigma, mu, omega)
  
  # Part II: Compute Share and log difference with j=0
  df_share <- df %>% group_by(t, j) %>% summarise(q = sum(q), s = sum(q)/N, .groups = "drop")
  
  # Part III: Compute log difference
  df_share <- df_share %>% 
    group_by(t) %>% 
    mutate(y = log(s/s[j == 0])) %>% 
    ungroup()
  
  # Merge with M
  df_share <- left_join(df_share, M,by = c("j", "t"))
  
  # Merge with x
  df_share <- left_join(df_share, x, by = "j")
  
  return(df_share)
}
df_share_smooth_delta <-
  compute_share_smooth_delta(
    df,
    delta, 
    sigma, 
    mu, 
    omega
  ) 
df_share_smooth_delta

# Find Delta:
solve_delta <- function(df_share_smooth, df, delta, sigma, mu, omega, kappa, lambda){
  diff <- 100  # Initialize diff to be larger than lambda
  delta_new <- delta  # Initialize delta to the initial values
  
  while (diff > lambda) {
    df_predicted <- compute_share_smooth_delta(df, delta_new, sigma, mu, omega)
    T_delta <- delta_new$delta + kappa * log(df_share_smooth$s / df_predicted$s)
    diff    <- max(abs(T_delta - delta_new$delta))
    print(diff)
    delta_new$delta <- T_delta  # Update delta for the next iteration
  }
  return(delta_new)
}

kappa <- 1
lambda <- 1e-2

delta_new <- solve_delta(df_share_smooth, df, delta, sigma, mu, omega, kappa, lambda)
summary(delta_new$delta - delta$delta)

# Use MCMC shocks to estimate delta
delta_new2 <-
  solve_delta(
    df_share_smooth, 
    df_mcmc,
    delta, 
    sigma, 
    mu, 
    omega, 
    kappa, 
    lambda
  )
summary(delta_new2$delta - delta$delta)

# Linear Parameter Estimation
compute_theta_linear <- function(df_share_smooth, delta, mu, omega, Psi) {
  # Extract relevant variables from the data
  df_share_smooth <- df_share_smooth %>% 
    dplyr::filter(j != 0) %>% 
    dplyr::arrange(t, j)
  
  delta <- delta %>% 
    dplyr::filter(j != 0) %>% 
    dplyr::arrange(t, j)
  
  X <- as.matrix(df_share_smooth[, grep("^x_", names(df_share_smooth))])
  c <- as.matrix(df_share_smooth$c)
  p <- as.matrix(df_share_smooth$p)
  delta_v <- as.matrix(delta$delta)
  W = cbind(X,c)
  
  # Compute alpha0
  alpha0 <- -exp(mu + omega^2 / 2)
  
  # Compute beta0
  beta_0 <- solve(t(X) %*% W %*% solve(Psi) %*% t(W) %*% X) %*% t(X) %*% W %*% solve(Psi) %*% t(W) %*% (delta_v - alpha0*p)
  
  return(beta_0)
}
Psi <- diag(length(beta) + 1)
theta_linear <-
  compute_theta_linear(
    df_share_smooth, 
    delta, 
    mu, 
    omega, 
    Psi
  ) 
cbind(
  theta_linear, 
  beta
)

# Compute the fixed effects
solve_xi <- function(df_share_smooth, delta, beta, mu, omega){
  delta$eps = delta$delta - (as.matrix(df_share_smooth[, grep("^x_", names(df_share_smooth))]) %*% beta - exp(mu + omega^2 / 2) * as.matrix(df_share_smooth$p))
  delta = delta %>% dplyr::filter(j != 0) %>% dplyr::arrange(t, j)
  xi  = delta$eps
  return(xi)
}
xi_new <- 
  solve_xi(
    df_share_smooth, 
    delta, 
    beta, 
    mu, 
    omega
  )
head(xi_new)

# GMM Objective: BLP Algorithm
theta_nonlinear <- c(sigma,mu,omega)

compute_gmm_objective_a4 <- function(theta_nonlinear, delta, df_share_smooth, Psi, df_mcmc,kappa,lambda){
  #Compute delta
  delta_new_inf <- solve_delta(df_share_smooth, df_mcmc, delta, theta_nonlinear[1:3],theta_nonlinear[4],theta_nonlinear[5], kappa, lambda)
  
  #Compute theta_1
  theta_1 <- compute_theta_linear(df_share_smooth, delta_new_inf, theta_nonlinear[4],theta_nonlinear[5], Psi)
  
  #Compute xi
  xi <- solve_xi(df_share_smooth, delta_new_inf, theta_1, theta_nonlinear[4], theta_nonlinear[5]) 
  xi <- as.matrix(xi)
  
  # Get W
  df_share_smooth <- df_share_smooth %>% dplyr::filter(j != 0) %>% dplyr::arrange(t, j)
  
  X <- as.matrix(df_share_smooth[, grep("^x_", names(df_share_smooth))])
  c <- as.matrix(df_share_smooth$c)
  W = cbind(X,c)
  
  # Compute Obj
  obj = t(xi) %*% W %*% solve(Psi) %*% t(W) %*% xi
  
  return(obj)
}

objective <-
  compute_gmm_objective_a4(
    theta_nonlinear, 
    delta, 
    df_share_smooth, 
    Psi, 
    df_mcmc,
    kappa,
    lambda
  ) 

# Final Estimation
result <-
  optim(par = theta_nonlinear,
        fn = compute_gmm_objective_a4,
        method = "Nelder-Mead",
        delta = delta, 
        df_share_smooth = df_share_smooth, 
        Psi = Psi, 
        df_mcmc = df_mcmc,
        kappa = kappa,
        lambda = lambda)
comparison <- cbind(theta_nonlinear, abs(result$par))
colnames(comparison) <- c("true", "estimate")
comparison

