options(scipen = 20)
library(doParallel)
library(doRNG)
library(dplyr)

######## Part I: Generate Data ########
registerDoParallel()
registerDoRNG() # register doRNG backend
set.seed(1) # set the seed
M <- 100 # number of markets
N <- 10 # the upper bound of the number of potential entrants
K <- 2 # the dimension of market attributes
L <- 2 # the dimension of potential entrant attributes
R <- 100 # the number of Monte Carlo simulations

# main parameters
beta <- abs(rnorm(K))
alpha <- abs(rnorm(L))
delta <- 1
rho <- abs(rnorm(1))

# auxiliary parameters
x_mu <- 1
x_sd <- 3
z_mu <- 0
z_sd <- 4

# number of potential entrants
E <- purrr::rdunif(M,1,N)
X <- matrix(rnorm(M * K, x_mu, x_sd),nrow = M)
colnames(X) <- paste("x", 1:K, sep = "_")

# entrant attributes
Z <- foreach(m = 1:M) %dorng% {
    Z_m <- matrix(rnorm(E[m] * L, z_mu, z_sd),nrow = E[m])
    colnames(Z_m) <- paste("z", 1:L, sep = "_")
    return(Z_m)
  }

# unobserved market attributes: epsilon
EP <- matrix(rnorm(M),nrow = M)
 
# unobserved entrant attributes: \nu
NU <-foreach(m = 1:M) %dorng% {
    NU_m <- matrix(rnorm(E[m]),nrow = E[m])
    return(NU_m)
  }
m <- 1
N_m <- dim(Z[[m]])[1]
y_m <- as.matrix(rep(1, N_m))
y_m[length(y_m)] <- 0
X_m <- X[m, , drop = FALSE]
Z_m <- Z[[m]]
EP_m <- EP[m, , drop = FALSE]
NU_m <- NU[[m]]

# payoff function
compute_payoff<-function(y_m, X_m, Z_m, EP_m, NU_m, beta, alpha, delta, rho){
  constant <- sum(beta * X_m) - delta * log(sum(y_m)) + rho*EP_m
  constant <- as.vector(constant)
  payoff <- sqrt(1-rho^2) * NU_m + constant + Z_m %*% alpha
  index_0 <- which(y_m == 0)
  payoff[index_0] <- 0
  return(payoff)
}

m <- 1
N_m <- dim(Z[[m]])[1]
y_m <- as.matrix(rep(1, N_m))
y_m[length(y_m)] <- 0
X_m <- X[m, , drop = FALSE]
Z_m <- Z[[m]]
EP_m <- EP[m, , drop = FALSE]
NU_m <- NU[[m]]


compute_payoff(
  y_m = y_m, 
  X_m = X_m, 
  Z_m = Z_m, 
  EP_m = EP_m, 
  NU_m = NU_m, 
  beta = beta, 
  alpha = alpha, 
  delta = delta, 
  rho = rho
)

# sequential entry order
compute_order_sequential_entry<-function(X_m, Z_m, EP_m, NU_m, beta, alpha, rho){
  constant <- sum(beta * X_m) + rho*EP_m
  constant <- as.vector(constant)
  payoff <- sqrt(1-rho^2) * NU_m + constant + Z_m %*% alpha
  return(order(payoff, decreasing = TRUE))
}
compute_order_sequential_entry(
  X_m = X_m,
  Z_m = Z_m, 
  EP_m = EP_m, 
  NU_m = NU_m, 
  beta = beta, 
  alpha = alpha, 
  rho = rho
)
# sequential entry outcome
compute_sequential_entry<-function(X_m, Z_m, EP_m, NU_m, beta, alpha, delta, rho){
  order    <- compute_order_sequential_entry(X_m, Z_m, EP_m, NU_m, beta, alpha, rho)
  constant <- sum(beta * X_m) + rho*EP_m
  constant <- as.vector(constant)
  payoff <- sqrt(1-rho^2) * NU_m + constant + Z_m %*% alpha
  entry_choice <- rep(0, length(payoff))
  for (o in order) {
    index_o <- which(order == o)
    if (payoff[o] - delta * log(index_o) > 0){
      entry_choice[o] <- 1
    }
  }
  return(entry_choice)
}
compute_sequential_entry(
  X_m = X_m, 
  Z_m = Z_m, 
  EP_m = EP_m, 
  NU_m = NU_m, 
  beta = beta, 
  alpha = alpha, 
  delta = delta, 
  rho = rho
)
# best response of simultaneous entry: rho = 0
compute_best_response_simultaneous_entry<-function(y_m,X_m, Z_m, EP_m, NU_m, beta, alpha, delta){
  number_entrant <- length(y_m)
  for (i in 1:number_entrant){
    y_m0 <- y_m
    y_m0[i] <- 1 - y_m0[i]
    payoff_entry <- compute_payoff(y_m, X_m, Z_m, EP_m, NU_m, beta, alpha, delta, 0)
    payoff_noentry <- compute_payoff(y_m0, X_m, Z_m, EP_m, NU_m, beta, alpha, delta, 0)
    payoff_entry_i <- payoff_entry[i]
    payoff_noentry_i <- payoff_noentry[i]
    if (payoff_noentry_i > payoff_entry_i) {
      y_m <- y_m0
    }
  }
  y_m <- as.matrix(y_m)
  return(y_m)
}
compute_best_response_simultaneous_entry(
  y_m = y_m,
  X_m = X_m, 
  Z_m = Z_m, 
  EP_m = EP_m, 
  NU_m = NU_m, 
  beta = beta, 
  alpha = alpha, 
  delta = delta
) 
# equilibrium of simultaneous entry based on best response: an iteration
compute_simultaneous_entry <- function(X_m, Z_m, EP_m, NU_m, beta, alpha, delta){
  # initialization
  number_entrant <- dim(Z_m)[1]
  initial_y_m <- rep(1, number_entrant)
  old_entry <- initial_y_m
  new_entry <- compute_best_response_simultaneous_entry(initial_y_m, X_m, Z_m, EP_m, NU_m, beta, alpha, delta)
  # loop until convergence
  while (!all(new_entry == old_entry)) {
    old_entry <- new_entry
    new_entry <- compute_best_response_simultaneous_entry(old_entry, X_m, Z_m, EP_m, NU_m, beta, alpha, delta)
  }
  return(new_entry)
}
compute_simultaneous_entry(
  X_m = X_m, 
  Z_m = Z_m, 
  EP_m = EP_m, 
  NU_m = NU_m, 
  beta = beta, 
  alpha = alpha,
  delta = delta
) 
# compute for the whole market: sequential entry
compute_sequential_entry_across_markets <- function(X, Z, EP, NU, beta, alpha, delta, rho) {
  entry_sequential <- list()
  for (m in 1:M) {
    X_m <- X[m, , drop = FALSE]
    Z_m <- Z[[m]]
    EP_m <- EP[m, , drop = FALSE]
    NU_m <- NU[[m]]
    entry_sequential[[m]] <- compute_sequential_entry(X_m, Z_m, EP_m, NU_m, beta, alpha, delta, rho)
  }
  return(entry_sequential)
}
Y_sequential <-
  compute_sequential_entry_across_markets(
    X = X, 
    Z = Z, 
    EP = EP, 
    NU = NU, 
    beta = beta, 
    alpha = alpha, 
    delta = delta, 
    rho = rho
  )

# compute for the whole market: payoffs for sequential entry
compute_payoff_across_markets <- function(Y, X, Z, EP, NU, beta, alpha, delta, rho) {
  payoff_sequential <- list()
  for (m in 1:M) {
    X_m <- X[m, , drop = FALSE]
    Z_m <- Z[[m]]
    EP_m <- EP[m, , drop = FALSE]
    NU_m <- NU[[m]]
    y_m <- Y[[m]]
    payoff_sequential[[m]] <- compute_payoff(y_m, X_m, Z_m, EP_m, NU_m, beta, alpha, delta, rho)
  }
  return(payoff_sequential)
}

# compute for the whole market: simultaneous entry
compute_simultaneous_entry_across_markets<-function(X, Z, EP, NU, beta, alpha, delta, rho = 0){
  entry_simultaneous <- list()
  for (m in 1:M) {
    X_m <- X[m, , drop = FALSE]
    Z_m <- Z[[m]]
    EP_m <- EP[m, , drop = FALSE]
    NU_m <- NU[[m]]
    entry_simultaneous[[m]] <- compute_simultaneous_entry(X_m, Z_m, EP_m, NU_m, beta, alpha, delta)
  }
  return(entry_simultaneous)
}
Y_simultaneous <-
  compute_simultaneous_entry_across_markets(
    X = X, 
    Z = Z, 
    EP = EP, 
    NU = NU, 
    beta = beta, 
    alpha = alpha, 
    delta = delta)
Y_simultaneous[[1]]
Y_simultaneous[[M]]
# compute for the whole market: payoffs for simultaneous entry
compute_payoff_simultaneous_entry_across_markets<-function(Y, X, Z, EP, NU, beta, alpha, delta, rho=0){
  payoff_simultaneous <- list()
  for (m in 1:M) {
    X_m <- X[m, , drop = FALSE]
    Z_m <- Z[[m]]
    EP_m <- EP[m, , drop = FALSE]
    NU_m <- NU[[m]]
    y_m <- Y[[m]]
    payoff_simultaneous[[m]] <- compute_payoff(y_m, X_m, Z_m, EP_m, NU_m, beta, alpha, delta, rho)
  }
  return(payoff_simultaneous)
}
payoff_simultaneous <-
  compute_payoff_across_markets(
    Y = Y_simultaneous, 
    X = X, 
    Z = Z, 
    EP = EP, 
    NU = NU, 
    beta = beta, 
    alpha = alpha, 
    delta = delta, 
    rho = 0
  )
min(unlist(payoff_simultaneous))
######## Part II: Estimating Parameters ########
# set seed
set.seed(1)

# unobserved market attributes
EP_mc <- foreach(r = 1:R) %dorng% {
  EP <- matrix(rnorm(M),nrow = M)
  return(EP)
}
head(EP_mc[[1]])

# unobserved entrant attributes: R times NU
NU_mc <- foreach(r = 1:R) %dorng% {
  NU <- foreach(m = 1:M) %do% {
    NU_m <- matrix(rnorm(E[m]),nrow = E[m])
    return(NU_m)
  }
  return(NU)
}
NU_mc[[1]][[1]]

# sequential entry
theta_sequential <- c(beta, alpha, delta,rho)
Y <- Y_sequential 

# Monte Carlo Simulation of sequential entry
compute_monte_carlo_sequential_entry<-function(X, Z, EP_mc, NU_mc, beta, alpha, delta, rho){
 Y_mc <- list()
 for (r in 1:R) {
    EP <- EP_mc[[r]]
    NU <- NU_mc[[r]]
    Y  <- compute_sequential_entry_across_markets(X, Z, EP, NU, beta, alpha, delta, rho)
    Y_mc[[r]] <- Y
  }
  return(Y_mc)
}
 
# compute objective function for sequential entry 
compute_objective_sequential_entry<-function(Y, X, Z, EP_mc, NU_mc, theta){
  Y_mc <- compute_monte_carlo_sequential_entry(X, Z, EP_mc, NU_mc, theta[1:2], theta[3:4], theta[5], theta[6])
  obj <- 0
  for (r in 1:R){
    for (m in 1:M) {
      obj <- obj +  sum(abs(Y_mc[[r]][[m]] - Y[[m]]))^2
    }
  }
  obj <- obj / (R * M)
  return(obj)
}
compute_objective_sequential_entry(Y, X, Z, EP_mc, NU_mc, theta_sequential)

# simultaneous entry
theta_simultaneous <- c(beta, alpha, delta)
Y <- Y_simultaneous

# Compute Monte Carlo Simulation of simultaneous entry
compute_monte_carlo_simultaneous_entry<-function(X, Z, EP_mc, NU_mc, beta, alpha, delta){
  Y_mc <- list()
  for (r in 1:R){
    EP_r <- EP_mc[[r]]
    NU_r <- NU_mc[[r]]
    Y_mc_r  <- compute_simultaneous_entry_across_markets(X, Z, EP_r, NU_r, beta, alpha, delta)
    Y_mc[[r]] <- Y_mc_r
  }
  return(Y_mc)
}
Y_mc <- 
  compute_monte_carlo_simultaneous_entry(
    X = X, 
    Z = Z, 
    EP_mc = EP_mc, 
    NU_mc = NU_mc, 
    beta = beta, 
    alpha = alpha, 
    delta = delta
  )
Y_mc[[1]][[1]]

# obj function for simultaneous entry 
compute_objective_simultaneous_entry<-function(Y, X, Z, EP_mc, NU_mc, theta){
  Y_mc <- compute_monte_carlo_simultaneous_entry(X, Z, EP_mc, NU_mc, theta[1:2], theta[3:4],theta[5])
  obj <- 0
  for (r in 1:R){
    for (m in 1:M) {
      obj <- obj +  sum(Y_mc[[r]][[m]] - Y[[m]])^2
    }
  }
  obj <- obj / (R * M)
  return(obj)
}
 
compute_objective_simultaneous_entry(
  Y = Y_simultaneous, 
  X = X, 
  Z = Z, 
  EP_mc = EP_mc, 
  NU_mc = NU_mc, 
  theta = theta_simultaneous
)

# estimate parameters: sequential
result_sequential <-
  optim(
    par = theta_sequential,
    fn = compute_objective_sequential_entry,
    method = "Nelder-Mead",
    Y = Y_sequential,
    X = X,
    Z = Z,
    EP_mc = EP_mc,
    NU_mc = NU_mc
  )
comparison <-
  data.frame(
    actual = theta_sequential,
    estimate = result_sequential$par
  )
comparison

# estimate parameters: simultaneous
result_simultaneous <-
  optim(
    par = theta_simultaneous,
    fn = compute_objective_simultaneous_entry,
    method = "Nelder-Mead",
    Y = Y_simultaneous,
    X = X,
    Z = Z,
    EP_mc = EP_mc,
    NU_mc = NU_mc)
comparison <-
  data.frame(
    actual = theta_simultaneous,
    estimate = result_simultaneous$par
  )
comparison
