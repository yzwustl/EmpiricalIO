library(doParallel)
library(dplyr)
library(reshape2)
library(ggplot2)
options(scipen = 20)
############# Part I: Generate Data #############
# set seed
set.seed(1)

# set constants 
L <- 5
K <- 1
T <- 100
N <- 1000
lambda <- 1e-10

# set parameters
alpha <- 0.5
beta <- 3
kappa <- 0.1
gamma <- 0.6
delta <- 0.95

# compute profits function (vectorized)
compute_pi <- function(alpha, beta, L, K) {
  s <- 1:L
  a <- 0:K
  pi <- as.vector(t(outer(alpha * log(s), -beta * a, "+")))
  pi <- matrix(pi, nrow = L * (K + 1), ncol = 1)
  return(pi)
}
PI <- compute_pi(alpha = alpha, beta = beta, L = L, K = K)

# compute transition probabilities
compute_G<-function(kappa, gamma, L, K){
  S <- 1:L
  A <- 0:K
  G <- matrix(0,nrow = (K+1)*L, ncol = L)
  for (a in A){
    for (s in S){
      if (1<s & s < L & a==0){
        G[2*s-1,s] <- 1 - kappa
        G[2*s-1,s-1] <- kappa
      }
      if (1<s & s<L & a==1){
        G[2*s,s] <- 1 - gamma - kappa
        G[2*s,s-1] <- kappa
        G[2*s,s+1] <- gamma
      }
      if (s==1 & a==0){
        G[s,s] <- 1
      }
      if (s==1 & a==1){
        G[2*s,s] <- 1 - gamma
        G[2*s,s+1] <- gamma
      }
      if (s==L & a==0){
        G[2*s-1,s] <- 1 - kappa
        G[2*s-1,s-1] <- kappa
      }
      if (s==L & a==1){
        G[2*s,s] <- 1 - kappa
        G[2*s,s-1] <- kappa
      }
    }
  }
  return(G)
}
G <- compute_G(kappa = kappa, gamma = gamma,  L = L, K = K)

# p vector: testing case
p <- matrix(rep(0.5, L * (K + 1)), ncol = 1)

# compute ex_ante value function
compute_exante_value<-function(p, PI, G, L, K, delta){
  # Term 3
  E_p <- -digamma(1) - log(p)
  E_p <- ifelse(is.finite(E_p), E_p, 0)
  term_3 <- PI + E_p
  # Term 2
  sigma_p <- matrix(0,nrow = L, ncol = (K+1)*L)
  for (l in 1:L) {
    sigma_p[l,c((K+1)*l-1,(K+1)*l)] <- p[c((K+1)*l-1,(K+1)*l)]
  }
  # Term 1
  term_1 <- solve(diag(L) - delta * sigma_p %*% G)
  # Ex-ante value function
  return(term_1 %*% sigma_p %*% term_3)
}
V <- compute_exante_value(p = p, PI = PI, G = G, L = L, K = K, delta = delta)
 
# compute choice value: section 7.2.8
compute_choice_value<-function(V, PI, G, delta){
  return(PI + delta * G %*% V)
}
value <- compute_choice_value(V = V, PI = PI, G = G, delta = delta)
 
# compute optimal belief 
compute_ccp<-function(V, PI, G, L, K, delta){
  # compute choice value
  value <- compute_choice_value(V = V, PI = PI, G = G, delta = delta)
  value <- as.data.frame(value)
  names(value)[1] <- "value"
  # value$a  01 01 01
  value$a <- rep(c(0:K), L)
  value$s <- rep(1:L, each = K + 1)
  # compute ccp
  ccp <- value %>%
    group_by(s) %>%
    mutate(p = exp(value) / sum(exp(value)))
  p <- ccp$p
  p <- as.matrix(p,nrow(L*(K+1)),ncol(1))
  return(p)
}
p <- compute_ccp(V = V, PI = PI, G = G, L = L, K = K, delta = delta)

# compute p and V by iteration:
solve_dynamic_decision<-function(PI, G, L, K, delta, lambda){
  # initial guess
  p <- matrix(rep(0.5, L * (K + 1)), ncol = 1)
  V <- compute_exante_value(p = p, PI = PI, G = G, L = L, K = K, delta = delta)
  diff = 100
  # iteration
  while (diff > lambda) {
    p_new <- compute_ccp(V = V, PI = PI, G = G, L = L, K = K, delta = delta)
    V_new <- compute_exante_value(p = p_new, PI = PI, G = G, L = L, K = K, delta = delta)
    diff <- max(abs(p_new - p))
    p <- p_new
    V <- V_new
    print(diff)
  }
  return(list(p = p, V = V))
}
output <- solve_dynamic_decision(PI = PI, G = G, L = L, K = K, delta = delta, lambda = lambda)
p <- output$p
V <- output$V

# simulate dynamic decision for a single firm
s <- 1
seed <- 1
simulate_dynamic_decision<-function(p, s, PI, G, L, K, T, delta, seed){
  set.seed(seed)
  df <- data.frame(t=1:T,s=s,a=0)
  for (t in 1:T) {
    # state
    s_t <- df[t, "s"]
    p_t <- p[c((K+1)*s_t-1,(K+1)*s_t)]
    # generate choice
    a_t <- rmultinom(1, 1, p_t)
    a_t <- which(a_t == 1) - 1
    df[t, "a"] <- a_t
    # simulate next state
    if (t<T){
      g_t <- G[(K + 1) * s_t + a_t - 1,]
      s_t_1 <- rmultinom(1, 1, g_t)
      s_t_1 <- which(s_t_1 == 1)
      df$s[t+1] <- s_t_1 
    }
  }
  df$i <- seed
  return(df)
}
 
# simulate dynamic decision for N firms
simulate_dynamic_decision_across_firms <- function(p, s, PI, G, L, K, T, N, delta) {
  df <- foreach(i = 1:N, .combine = rbind) %dopar% {
    df_i <- simulate_dynamic_decision(p = p, s = s, PI = PI, G = G, L = L, K = K, T = T, delta = delta, seed = i)
    return(df_i)
  }
  return(df)
}
df <- simulate_dynamic_decision_across_firms(
    p = p, 
    s = s, 
    PI = PI,
    G = G, 
    L = L, 
    K = K, 
    T = T, 
    N = N, 
    delta = delta
  )

# estimate ccp
estimate_ccp<-function(df){
  ccp <- df %>%
    group_by(s,a) %>%
    summarise(p = length(a)) %>%
    mutate(p = p / sum(p))
  p <- ccp$p
  p <- as.vector(p)
  return(p)
}
p_est <- estimate_ccp(df = df)
check_ccp <- cbind(p, p_est)
colnames(check_ccp) <- c("true", "estimate")
check_ccp <- check_ccp %>% reshape2::melt()
ggplot(
  data = check_ccp, 
  aes(
    x = Var1, 
    y = value, 
    fill = Var2
  )
) +
  geom_bar(
    stat = "identity", 
    position = "dodge"
  ) +
  labs(fill = "Value") + 
  xlab("action/state") + 
  ylab("probability")  + 
  theme_classic()

# estimate transition matrix
estimate_G<-function(df){
  df$s_1 <- df %>% group_by(i) %>% mutate(s_1 = lead(s)) %>% pull(s_1)
  df <- df[!is.na(df$s_1),]
  G <- df %>%
    group_by(s,a,s_1) %>%
    summarise(n = length(s_1)) %>%
    mutate(n = n / sum(n))
  G <- G[!is.na(G$s_1),]
  K <- max(G$a)
  L <- max(G$s)
  n_comb <- dim(G)[1]
  G_est <- matrix(0,nrow = (K+1)*L, ncol = L)
  for (i in 1:n_comb) {
    s <- G$s[i]
    a <- G$a[i]
    s_1 <- G$s_1[i]
    G_est[(K+1)*s+a-1,s_1] <- G$n[i]
  }
  return(G_est)
}
G_est <- estimate_G(df = df)

#### Part II: Estimate Parameters ####
theta_1 <- c(alpha,beta)
theta_2 <- c(kappa,gamma)
theta <- c(theta_1,theta_2)

# loglikelihood function for NFP
compute_loglikelihood_NFP<-function(theta, df, delta, L, K){
  alpha <- theta[1]
  beta <- theta[2]
  kappa <- theta[3]
  gamma <- theta[4]
  # PI
  PI <- compute_pi(alpha = alpha, beta = beta, L = L, K = K)
  # G
  G <- compute_G(kappa = kappa, gamma = gamma,  L = L, K = K)
  G_extension <- expand.grid(s = 1:L, a = 0:K, s_1 = 1:L)
  G_extension$ccp <- 0
  for (i in 1:dim(G_extension)[1]) {
    s <- G_extension$s[i]
    a <- G_extension$a[i]
    s_1 <- G_extension$s_1[i]
    G_extension$ccp[i] <- G[(K+1)*s+a-1,s_1]
  }
  # solve dynamic decision
  output <- solve_dynamic_decision(PI = PI, G = G, L = L, K = K, delta = delta, lambda = lambda)
  p <- output$p
  p_extension <- expand.grid(s = 1:L, a = 0:K)
  p_extension$p <- 0
  for (i in 1:dim(p_extension)[1]) {
    s <- p_extension$s[i]
    a <- p_extension$a[i]
    p_extension$p[i] <- p[(K+1)*s+a-1]
  }
  # df
  df$s_1 <- df %>% group_by(i) %>% mutate(s_1 = lead(s)) %>% pull(s_1)
  df <- df %>% left_join(G_extension, by = c("s","a","s_1"))
  df <- df %>% left_join(p_extension, by = c("s","a"))
  df$ccp[is.na(df$ccp)] <- 1
  # log likelihood
  df$loglikelihood <- log(df$ccp) + log(df$p)
  loglikelihood <- mean(df$loglikelihood)
  return(loglikelihood)
}
likelihood <- compute_loglikelihood_NFP(theta = theta, df = df, delta = delta, L = L, K = K)

lower <- rep(0, length(theta))
upper <- c(1, 5, 0.2, 0.7)
NFP_result <-
  optim(
    par = theta,
    fn = compute_loglikelihood_NFP,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    control = list(fnscale = -1),
    df = df, 
    delta = delta,
    L = L,
    K = K
  )
compare <-
  data.frame(
    true = theta,
    estimate = NFP_result$par
  ); 
compare

# estimate theta_2
estimate_theta_2<-function(df){
  G_est <- estimate_G(df = df)
  G_extension <- expand.grid(s = 1:L, a = 0:K, s_1 = 1:L)
  G_extension$ccp <- 0
  for (i in 1:dim(G_extension)[1]) {
    s <- G_extension$s[i]
    a <- G_extension$a[i]
    s_1 <- G_extension$s_1[i]
    G_extension$ccp[i] <- G_est[(K+1)*s+a-1,s_1]
  }
  G_kappa <- G_extension %>%
    filter((s == max(s) & s_1 == max(s) - 1) |
             (s > min(s) & s < max(s) & a == 0 & s_1 == s - 1) |
             (s > min(s) & s < max(s) & a == 1 & s_1 == s - 1)) 
  kappa_est <- mean(G_kappa$ccp)
  G_gamma <- G_extension %>%
    filter((s == min(s) & s_1 == min(s) + 1 & a==1) |
            (s > min(s) & s < max(s) & a == 1 & s_1 == s + 1))
  gamma_est <- mean(G_gamma$ccp)
  return(c(kappa_est, gamma_est))
}
theta_2_est <- estimate_theta_2(df = df)

# obj function for CCP
compute_CCP_objective<-function(theta_1, theta_2, p_est, L, K, delta){
  alpha <- theta_1[1]
  beta <- theta_1[2]
  kappa <- theta_2[1]
  gamma <- theta_2[2]
  # PI
  PI <- compute_pi(alpha = alpha, beta = beta, L = L, K = K)
  # G
  G <- compute_G(kappa = kappa, gamma = gamma,  L = L, K = K)
  # compute ccp
  V <- compute_exante_value(p = p_est, PI = PI, G = G, L = L, K = K, delta = delta)
  ccp <- compute_ccp(V = V, PI = PI, G = G, L = L, K = K, delta = delta)
  # objective function
  obj <- mean((ccp - p_est)^2)
  return(obj)
}
lower <- rep(0, length(theta_1))
upper <- c(1, 5)
CCP_result <-
  optim(
    par = theta_1,
    fn = compute_CCP_objective,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    theta_2 = theta_2_est, 
    p_est = p_est, 
    L = L, 
    K = K, 
    delta = delta
  )
compare <-
  data.frame(
    true = theta_1,
    estimate = CCP_result$par
  ); 
compare

