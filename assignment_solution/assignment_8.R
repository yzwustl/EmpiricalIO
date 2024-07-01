library(doParallel)
library(dplyr)
library(reshape2)
library(Matrix)
library(ggplot2)
options(scipen = 200)
########## Part I: Simulate the data ##########
# set seed
set.seed(1)
# set constants 
L <- 5     
K <- 1     
T <- 100
N <- 3      
M <- 1000  
lambda <- 1e-10
# set parameters
alpha <- 1
eta <- 0.3
beta <- 2
kappa <- 0.1
gamma <- 0.6
delta <- 0.95
# action space and state space
compute_action_state_space<-function(K, L, N){
  n_actions <- (K+1)^N
  n_states  <- L^N
  actions <- 0:K
  states <- 1:L
  # Action space
  A <- data.frame(k=rep(1:n_actions,each=N),i = rep(1:N,n_actions))
  expanded_actions <- expand.grid(rep(list(actions), N))
  A$a <- 0
  for (i in 1:n_actions){
    A[A$k==i,]$a <- as.vector(expanded_actions[i,])
  } 
  # State Space
  S <- data.frame(l=rep(1:n_states,each=N),i = rep(1:N,n_states))
  expanded_states <- expand.grid(rep(list(states), N))
  S$s <- 0
  for (i in 1:n_states){
    S[S$l==i,]$s <- as.vector(expanded_states[i,])
  }
  return(list(A=A,S=S))
}
output <- compute_action_state_space(K, L, N)
A <- output$A
A$a <- as.numeric(A$a)
S <- output$S
S$s <- as.numeric(S$s)
m_a <- max(A$k) 
m_s <- max(S$l)

# profits function
compute_PI_game <- function(alpha, beta, eta, A, S) {
  m_a <- max(A$k)
  m_s <- max(S$l)
  N <- max(A$i)
  PI <- vector("list", length = N)
  for (n in 1:N) {
    PI[[n]] <- rep(1, m_a * m_s)
  }
  for (j in 1:m_s) {
    for (i in 1:m_a) {
      A_i <- A[A$k == i, ]
      S_j <- S[S$l == j, ]
      data_ij <- merge(A_i, S_j, by = "i")
      data_ij$pi <- alpha * log(data_ij$s) - eta * log(data_ij$s) * (sum(log(data_ij$s)) - log(data_ij$s)) - beta * data_ij$a
      for (n in 1:N) {
        PI[[n]][(j - 1) * m_a + i] <- data_ij$pi[data_ij$i == n]
      }
    }
  }
  return(PI)
}
PI <- compute_PI_game(alpha, beta, eta, A, S)

# marginal transitional matrix
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
G_marginal <- compute_G(kappa = kappa,gamma = gamma, L = L, K = K)

# joint transitional matrix
compute_G_game <- function(G_marginal, A, S) {
  m_a <- max(A$k)
  m_s <- max(S$l)
  n_firm <- max(A$i)
  output <- matrix(0, nrow = m_a * m_s, ncol = m_s)
  # Loop to get output
  for (k in 1:m_a) {
    for (l in 1:m_s) {
      A_k <- A$a[A$k == k]
      S_l <- S$s[S$l == l]
      transitions <- list()
      for (i in 1:n_firm) {
        sa <- c(A_k[i], S_l[i])
        transitions[[i]] <- G_marginal[(sa[2] - 1) * 2 + sa[1] + 1,]
      }
      # Compute the Cartesian product of transitions and their product
      transition_grid <- expand.grid(transitions)
      transition_grid$t <- apply(transition_grid, 1, prod)
      output[(l - 1) * m_a + k, ] <- transition_grid$t
    }
  }
  return(output)
}
G <- 
  compute_G_game(
    G_marginal = G_marginal, 
    A = A, 
    S = S
  )

# Initial marginal distribution
initialize_p_marginal<-function(A, S){
  p_marginal <- left_join(A, S, by = "i", relationship = "many-to-many")
  p_marginal <- p_marginal %>% select(i,l,a) 
  p_marginal <- unique(p_marginal)
  p_marginal$p <- 0.5
  # order by i,l
  p_marginal <- p_marginal %>% arrange(i, l)
  return(p_marginal)
}
p_marginal <- 
  initialize_p_marginal(
    A = A,
    S = S
  )

# from marginal to joint distribution
compute_p_joint<-function(p_marginal, A, S){
  n_action <- max(A$k)
  n_state <- max(S$l)
  p_joint <- merge(S, A, by = "i")
  p_joint <- p_joint %>% select(i,l,a,k)
  p_joint <- left_join(p_joint,p_marginal,by=c('i','l','a'))
  p_joint <- p_joint %>% group_by(l,k) %>% mutate(p = prod(p))
  p_joint <- p_joint %>% select(l,k,p) %>% unique() %>% arrange(l,k)
  return(p_joint)
}
p_joint <- 
  compute_p_joint(
    p_marginal = p_marginal, 
    A = A, 
    S = S
  )

# from joint to marginal: just to check
compute_p_marginal <- function(p_joint, A, S){
  N <- max(S$i)
  m_a <- max(A$k)
  m_s <- max(S$l)
  p_marginal2 <-
    expand.grid(
      i = 1:N, 
      l = 1:m_s,
      k = 1:m_a
    )
  p_marginal2 <- left_join(p_marginal2, p_joint, by = c("l", "k"))
  p_marginal2 <- left_join(p_marginal2,A,by = c('i','k'))
  p_marginal2 <- p_marginal2 %>% group_by(i, l, a) %>% summarise(p = sum(p))
  return(p_marginal2)
}
 
# compute sigma
compute_sigma <- function(p_marginal, A, S){
  p_joint <- compute_p_joint(p_marginal, A, S)
  m_s <- max(S$l)
  p_joint <- foreach (ll = 1:m_s) %do% {
      p_joint_l <- p_joint %>% filter(ll == l) %>% arrange(k)
      p_joint_l <- t(matrix(p_joint_l$p))
      return(p_joint_l)
  }
  sigma <- bdiag(p_joint)
  return(sigma)
}
sigma <- compute_sigma(p_marginal, A, S)

# compute D
compute_D <- function(p_marginal){
 p_marginal$E <- - digamma(1) - log(p_marginal$p)
 p_marginal$pE <- p_marginal$p * p_marginal$E
 D <- foreach(n = 1:N) %do% {
   p_marginal_n <- p_marginal %>% filter(i == n)
   output <- p_marginal_n %>% select(l,pE) %>% group_by(l) %>% summarise(D = sum(pE))
   D_n <- output$D
   return(D_n)
 }
 return(D)
}
D <- compute_D(p_marginal)
 
# compute ex-ante value
compute_exante_value_game<-function(p_marginal, A, S, PI, G, delta){
  sigma <- compute_sigma(p_marginal, A, S)
  D     <- compute_D(p_marginal)
  # term 1
  term_1 <- diag(dim(sigma)[1]) - delta * sigma %*% G
  V <- foreach (n = 1:length(D)) %do% {
    PI_n <- PI[[n]]
    D_n <- D[[n]]
    term_2_n <- sigma %*% PI_n + D_n
    V_n <- solve(term_1) %*% term_2_n
    return(V_n)
  }
  return(V)
}
V <- 
  compute_exante_value_game(
    p_marginal = p_marginal,
    A = A,
    S = S,
    PI = PI,
    G = G,
    delta = delta
  )

# compute state-action-profile value function
compute_profile_value_game<-function(V, PI, G, delta, S, A){
  m_s <- max(S$l)
  m_a <- max(A$k)
  value <- foreach(n = 1:length(V),.combine = 'rbind') %do% {
    V_n  <-  V[[n]]
    PI_n <- PI[[n]]
    value_n <- PI_n + delta * G %*% V_n
    header <- expand.grid(i = n,l = 1:m_s,k = 1:m_a) %>% dplyr::arrange(i,l,k)
    value_n <- data.frame(header,value = as.numeric(value_n))
    return(value_n)
  }
  return(value)
}
value <- 
  compute_profile_value_game(
    V = V, 
    PI = PI,
    G = G, 
    delta = delta, 
    S = S, 
    A = A
  )
compute_choice_value_game<-function(p_marginal, V, PI, G, delta, A, S){
  value <- compute_profile_value_game(V, PI, G, delta, S, A)
  p_joint <- compute_p_joint(p_marginal, A, S) 
  value <- value %>%
    left_join(p_joint, by = c("l", "k")) %>%
    rename(p_joint = p) %>%
    left_join(A, by = c("i", "k")) %>%
    left_join(p_marginal, by = c("i", "l", "a")) %>%
    mutate(p_others = p_joint / p) %>%
    mutate(value_p = value * p_others) %>%
    group_by(i, l, a) %>% summarise(value = sum(value_p))
 return(value)
}
value <- compute_choice_value_game(p_marginal, V, PI, G, delta, A, S)

# compute conditional choice probability
compute_ccp_game<-function(p_marginal, V, PI, G, delta, A, S){
  value <- compute_choice_value_game(p_marginal, V, PI, G, delta, A, S)
  p_marginal <- value %>% group_by(i,l) %>% mutate(p = exp(value) / sum(exp(value))) %>%
    select(i,l,a,p) %>% arrange(i,l,a) %>% ungroup()
  return(p_marginal)
}
 
# iterate until V converges to solve for p_marginal and V
solve_dynamic_game<-function(PI, G, L, K, delta, lambda, A, S){
  p_marginal <- initialize_p_marginal(A, S)
  V <- compute_exante_value_game(p_marginal, A, S, PI, G, delta)
  # initialize distance
  distance <- 100
  # loop until V converges
  while (distance > lambda){
    V_old <- V
    p_marginal <- compute_ccp_game(p_marginal, V, PI, G, delta, A, S)
    V <- compute_exante_value_game(p_marginal, A, S, PI, G, delta)
    V_check <- do.call(rbind, V)
    V_check_old <- do.call(rbind, V_old)
    distance <- max(abs(unlist(V_check) - unlist(V_check_old)))
    print(distance)
  }
  return(list(V = V, p_marginal = p_marginal))
}
output <-
  solve_dynamic_game(
    PI = PI,
    G = G, 
    L = L,
    K = K,
    delta = delta,
    lambda = lambda,
    A = A,
    S = S
  )
p_marginal <- output$p_marginal 
V <- output$V[[N]] 
p_joint <- 
  compute_p_joint(
    p_marginal = p_marginal, 
    A = A, 
    S = S
  )

# Simulate the dynamic data
seed <- 1
l <- 1
simulate_dynamic_game<-function(p_joint, l, G, N, T, S, A, seed){
  set.seed(seed)
  df <- expand.grid(t = 1:T, i=1:N, l=l, k=1)
  for (tt in 1:max(df$t)) {
    # state
    l_t <- df %>% dplyr::filter(t == tt) 
    l_t <- unique(l_t$l)
    # draw action
    p_t <- p_joint %>% dplyr::filter(l == l_t)
    k_t <- rmultinom(n = 1, size = 1, prob = p_t$p) %>% as.numeric()
    k_t <- which(k_t == 1)
    df[df$t == tt, "k"] <- rep(k_t, N)
    # next state
    g_t <- G[m_a * (l_t - 1) + k_t, ]
    l_t_next <- which(rmultinom(n = 1, size = 1, prob = g_t) == 1)
    df[df$t == tt+1, "l"] <- rep(l_t_next, N)
  }
  df <- 
    df %>% 
    left_join(A, by = c("i", "k")) %>% 
    left_join(S, by = c("i", "l")) %>% 
    select(t, i, l, k, s, a) %>%
    arrange(t,i)
  return(df)
}
 
# decisions
simulate_dynamic_decision_across_markets<-function(p_joint, l, G, N, T, M, S, A){
  df <- foreach (mm=1:M,.combine = 'rbind') %dopar% {
    seed <- mm
    df_m <- simulate_dynamic_game(p_joint = p_joint,l = l,G = G,N = N,T = T,S = S, A = A,seed = seed)
    df_m <- data.frame(m = mm,df_m) 
    return(df_m)
  }
  return(df)
}
df <- 
  simulate_dynamic_decision_across_markets(
    p_joint = p_joint,
    l = l,
    G = G,
    N = N,
    T = T,
    M = M,
    S = S,
    A = A
  )
saveRDS(df, "df.rds")
df <- readRDS("df.rds")
# estimate p_marginal based on df
estimate_ccp_marginal_game<-function(df){
  p_marginal_est <- df %>% select(i,l,a) %>% arrange(i,l,a) %>% 
    group_by(i,l,a) %>% mutate(n = n()) %>% unique() %>%
    group_by(i,l) %>% mutate(p = n / sum(n)) %>% ungroup() %>%
    select(i,l,a,p)
  return(p_marginal)
}
p_marginal_est <- estimate_ccp_marginal_game(df = df)
check_ccp <- p_marginal_est %>% rename(estimate = p) %>%
  left_join(p_marginal, by = c("i","l", "a")) %>%
  rename(true = p) %>%
  filter(a == 1)
ggplot(data = check_ccp, aes(x = true, y = estimate)) +
  geom_point() +
  labs(fill = "Value") +
  xlab("true") + 
  ylab("estimate")  + 
  theme_classic()

# estimate G_marginal based on df
estimate_G_marginal <- function(df){
  n_action <- max(df$a)+1
  n_state <-  max(df$s)
  G_marginal_est <- matrix(0, nrow = n_action * n_state, ncol = n_state)
  G <- df %>% arrange(m,i,t) %>% 
    group_by(m,i) %>% mutate(s_lead = lead(s, 1)) %>% 
    filter(!is.na(s_lead)) %>% dplyr::ungroup() %>%
    group_by(s,a,s_lead) %>% mutate(n = n()) %>% select (s,a,s_lead,n) %>% 
    arrange(s,a) %>% unique() %>% ungroup() %>% group_by(s,a) %>%
    mutate(p = n / sum(n)) %>% ungroup() %>% select(s,a,s_lead,p)
  n_comb <- dim(G)[1]
  for (i in 1:n_comb) {
    s <- G$s[i]
    a <- G$a[i]
    s_1 <- G$s_lead[i]
    G_marginal_est[n_action*s+a-1,s_1] <- G$p[i]
  }
  return(G_marginal_est)
}

G_marginal_est <- estimate_G_marginal(df = df)
check_G <- 
  data.frame(
    type = "true",
    reshape2::melt(G_marginal)
  )
check_G_est <- 
  data.frame(
    type = "estimate",
    reshape2::melt(G_marginal_est)
  )
check_G <- 
  rbind(
    check_G, 
    check_G_est
  )
check_G$variable <- 
  paste(
    check_G$Var1,
    check_G$Var2, 
    sep = "_"
  )
ggplot(
  data = check_G, 
  aes(
    x = variable, 
    y = value,
    fill = type
  )
) +
  geom_bar(
    stat = "identity",
    position = "dodge"
  ) +
  labs(fill = "Value") + 
  xlab("action/state/state") +
  ylab("probability") +
  theme(axis.text.x = element_blank())  + 
  theme_classic()

# Part II: Estimation
theta_1 <- c(alpha, beta, eta)
theta_2 <- c(kappa, gamma)
theta <- c(theta_1, theta_2)
 
estimate_theta_2_game<-function(df){
  G_est <- estimate_G_marginal(df = df)
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
theta_2_est <- estimate_theta_2_game(df) 

# objective function
compute_CCP_objective_game<-function(theta_1, theta_2, p_est, A, S, delta,lambda){
  alpha <- theta_1[1]
  beta <- theta_1[2]
  eta <- theta_1[3]
  kappa <- theta_2[1]
  gamma <- theta_2[2]
  L <- max(S$s)
  K <- max(A$a)
  # PI
  PI <- compute_PI_game(alpha, beta, eta, A,S)
  # G
  G_marginal <- compute_G(kappa,gamma,L,K)
  G <- compute_G_game(G_marginal,A,S)
  # p_marginal
  V <-  compute_exante_value_game(p_marginal_est, A, S, PI, G, delta)
  p_marginal <- compute_ccp_game(p_marginal_est,V,PI,G,delta,A,S)
  # compare p_marginal and p_marginal_est
  distance <- 
    p_marginal %>%
    dplyr::rename(ccp = p) %>%
    dplyr::left_join(
      p_marginal_est,
      by = c("i", "l", "a")
    ) %>% 
    dplyr::mutate(x = (ccp - p)^2) %>%
    dplyr::summarise(
      mean(
        x, 
        na.rm = TRUE
      )
    ) %>%
    as.numeric()
  return(distance)
}
objective <- 
  compute_CCP_objective_game(
    theta_1 = theta_1,
    theta_2 = theta_2, 
    p_est = p_marginal_est,
    A = A,
    S = S, 
    delta = delta, 
    lambda = lambda
  ) 

# Estimation 
lower <-rep(0, length(theta_1))
upper <- c(1, 5, 0.3)
CCP_result <-
  optim(
    par = theta_1,
    fn = compute_CCP_objective_game,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    theta_2 = theta_2_est,
    p_est = p_marginal_est,
    A = A,
    S = S,
    delta = delta,
    lambda = lambda
  )
CCP_result
compare <-
  data.frame(
    true = theta_1,
    estimate = CCP_result$par
  ); 
compare

