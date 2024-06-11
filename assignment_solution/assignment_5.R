new_tempdir <- "C:/Your/Custom/TempDir"
if (!dir.exists(new_tempdir)) {
  dir.create(new_tempdir, recursive = TRUE)
}
tempdir(new_tempdir)
library(dplyr)
library(evd)
library(purrr)
library(doParallel)
library(foreach)
library(ggplot2)
options(scipen = 20)
######## Part I: Generate Data ########
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

# set auxiliary parameters
price_xi <- 1
sd_x <- 2
sd_xi <- 0.5
sd_c <- 0.05
sd_p <- 0.05

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
M$p <- 0
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

# Smooth Derivative Function
compute_derivative_share_smooth <- function(df, beta, sigma, mu, omega) {
  # Part I: Compute Choice
  df <- compute_choice_smooth(df, beta, sigma, mu, omega)
  df <- df %>% filter(j != 0)
  df$alpha_i <- -exp(mu + omega * df$v_p)
  
  # Part II: For each market, compute the derivative of share with respect to price, put results in a list
  T <- length(unique(df$t))
  
  # Initialize parallel backend
  cores <- detectCores()
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  result <- foreach(l = 1:T, .packages = c("dplyr")) %dopar% {
    data_j <- df %>% filter(t == l)
    J_t <- length(unique(data_j$j))
    matrix <- matrix(0, nrow = J_t, ncol = J_t)
    j_list <- unique(data_j$j)
    
    for (m in 1:J_t) {
      for (n in 1:J_t) {
        if (m == n) {
          data_jk <- data_j %>% filter(j == j_list[m])
          data_jk$obj <- data_jk$alpha_i * data_jk$q * (1 - data_jk$q)
          matrix[m, n] <- mean(data_jk$obj)
        } else {
          data_jk <- data_j %>% filter(j == j_list[m] | j == j_list[n])
          data_jk <- data_jk %>%
            arrange(i, j) %>%
            group_by(i) %>%
            mutate(lag_q = lag(q)) %>%
            ungroup() %>%
            filter(!is.na(lag_q))
          data_jk$obj <- data_jk$alpha_i * data_jk$q * data_jk$lag_q
          matrix[m, n] <- -mean(data_jk$obj)
        }
      }
    }
    
    matrix
  }
  
  # Stop parallel backend
  stopCluster(cl)
  
  return(result)
}

# Delta Matrix: Identity Matrix
delta <- 
  foreach (tt = 1:T) %do% {
    J_t <- M %>%
      dplyr::filter(t == tt) %>%
      dplyr::filter(j > 0) 
    J_t <- dim(J_t)[1]
    Delta_t <- diag(rep(1, J_t))
    return(Delta_t)
  }

# Update Price Function
update_price<-function(p, x, M, V, beta, sigma, mu, omega,delta){
  # Part I: construct df based on the p
  M[M$j != 0, "p"] <- p
  df <- expand.grid(t = 1:T, i = 1:N, j = 0:J)
  df <- left_join(df, M, by = c("j", "t"))
  df <- left_join(df, V, by = c("i", "t"))
  df <- left_join(df, x, by = "j")
  df <- df %>% dplyr::filter(!is.na(p)) %>% dplyr::arrange(t,i,j)
  # Part II: compute share and derivative of share
  df_share <- compute_share_smooth(df, beta, sigma, mu, omega)
  df_share <- df_share %>% dplyr::filter(j != 0)
  derivative_share <- compute_derivative_share_smooth(df, beta, sigma, mu, omega)
  # Part III: compute the update of price
  M <- M %>% dplyr::filter(j != 0)
  for (l in 1:T) {
    c_t <- M %>% dplyr::filter(t == l) %>% dplyr::select(c) %>% as.matrix()
    omega_t <- delta[[l]] * derivative_share[[l]]
    s_t <- df_share %>% dplyr::filter(t == l) %>% dplyr::select(s) %>% as.matrix()
    pt1 = c_t + solve(-omega_t) %*% s_t
    M[M$t == l, "p"] <- pt1
  }
  new_p <- M$p
  return(new_p)
}

# set the threshold
lambda <- 1e-6
# set the initial price
p <- M[M$j > 0, "p"]
p_init <- rep(1,dim(p)[1])
p_new <- update_price(
  p = p_init, 
  x = x, 
  M = M, 
  V = V, 
  beta = beta, 
  sigma = sigma, 
  mu = mu, 
  omega = omega,
  delta = delta
)
distance <- 10000
while (distance > lambda) {
  p_old <- p_new
  p_new <- 
    update_price(
      p_old, 
      x, 
      M, 
      V, 
      beta, 
      sigma, 
      mu, 
      omega,
      delta
    )
  distance <- max(abs(p_new - p_old))
  print(distance)
}
p_actual <- p_new
################# Part II: Cost Estimation #################
estimate_marginal_cost <- function(p = p, x = x, M = M, V = V, beta = beta,  sigma = sigma, mu = mu, omega = omega, delta = delta){
  # Part I: construct df based on the p
  M[M$j != 0, "p"] <- p_actual
  df <- expand.grid(t = 1:T, i = 1:N, j = 0:J)
  df <- left_join(df, M, by = c("j", "t"))
  df <- left_join(df, V, by = c("i", "t"))
  df <- left_join(df, x, by = "j")
  df <- df %>% dplyr::filter(!is.na(p)) %>% dplyr::arrange(t,i,j)
  # Part II: compute share and derivative of share
  df_share <- compute_share_smooth(df, beta, sigma, mu, omega)
  df_share <- df_share %>% dplyr::filter(j != 0)
  derivative_share <- compute_derivative_share_smooth(df, beta, sigma, mu, omega)
  # Part III: compute the cost estimation
  M <- M %>% dplyr::filter(j != 0)
  for (l in 1:T) {
    p_t <- M %>% dplyr::filter(t == l) %>% dplyr::select(p) %>% as.matrix()
    omega_t <-  delta[[l]] * derivative_share[[l]]
    s_t <- df_share %>% dplyr::filter(t == l) %>% dplyr::select(s) %>% as.matrix()
    ct = p_t - solve(-omega_t,s_t)
    M[M$t == l, "c"] <- ct
  }
  new_c <- M$c
  return(new_c)
}
marginal_cost_estimate <- 
  estimate_marginal_cost(
    p = p_new, 
    x = x, 
    M = M, 
    V = V, 
    beta = beta, 
    sigma = sigma, 
    mu = mu, 
    omega = omega, 
    delta = delta)
marginal_cost_actual <- M[M$j > 0, ]$c

# plot the estimate vs actual marginal costs
marginal_cost_df <-
  data.frame(
    actual = marginal_cost_actual,
    estimate = marginal_cost_estimate
  )
ggplot(
  marginal_cost_df, 
  aes(
    x = estimate, 
    y = actual
  )
) +
  geom_point() + 
  theme_classic()

################# Part III: Counter factual Analysis #################
# Delta Matrix: Company 1 merges with 2 and 3
delta_counterfactual <- 
  foreach (tt = 1:T) %do% {
    J_t <- M %>% dplyr::filter(t == tt) %>% dplyr::filter(j > 0)
    J_list <- J_t$j
    Jt <- dim(J_t)[1]
    Delta_t <- diag(rep(1, Jt))
    for (j in J_list) {
      if (j == 1 | j == 2 | j == 3) {
        j_index <- which(J_list == j)
        for (l in J_list) {
          if (l == 1 | l == 2 | l == 3) {
            l_index <- which(J_list == l)
            Delta_t[j_index,l_index] <- 1
          }
        }
      }
    }
    return(Delta_t)
  }

# Compute conterfactual price
p_new <- 
  update_price(
    p = p_actual, 
    x = x, 
    M = M, 
    V = V, 
    beta = beta, 
    sigma = sigma, 
    mu = mu, 
    omega = omega, 
    delta = delta_counterfactual
  )
distance <- 10000
while (distance > lambda) {
  p_old <- p_new
  p_new <- 
    update_price(
      p_old,
      x, 
      M, 
      V, 
      beta, 
      sigma, 
      mu, 
      omega, 
      delta_counterfactual
    )
  distance <- max(abs(p_new - p_old))
  print(distance)
}
p_counterfactual <- p_new

# Compute price change
M_counterfactual<-M
M_counterfactual$p_actual<-0
M_counterfactual[M_counterfactual$j>0,]$p_actual<-p_actual

M_counterfactual$p_counterfactual<-0
M_counterfactual[M_counterfactual$j>0,]$p_counterfactual<-p_counterfactual 

p_change<-M_counterfactual[M_counterfactual$j>0,] %>%
  dplyr::mutate(p_change=(p_counterfactual-p_actual)/p_actual)%>%
  dplyr::group_by(j)%>%
  dplyr::summarise(p_change=mean(p_change))%>%
  dplyr::ungroup()
 
# Producer Surplus Function
compute_producer_surplus<-function(p, marginal_cost, x, M, V, beta, sigma, mu, omega){
  # Part I: construct df based on the p
  M[M$j != 0, "p"] <- p
  df <- expand.grid(t = 1:T, i = 1:N, j = 0:J)
  df <- left_join(df, M, by = c("j", "t"))
  df <- left_join(df, V, by = c("i", "t"))
  df <- left_join(df, x, by = "j")
  df <- df %>% dplyr::filter(!is.na(p)) %>% dplyr::arrange(t,i,j)
  # Part II: compute share and derivative of share
  df_share <- compute_share_smooth(df, beta, sigma, mu, omega)
  df_share <- df_share %>% dplyr::filter(j != 0)
  # Part III: compute the producer surplus
  M[M$j != 0, "c"] <- marginal_cost
  M[M$j != 0, "s"] <- df_share$s
  M <- M %>% dplyr::filter(j != 0)
  producer_surplus <- M$s * (M$p - M$c)
  return(producer_surplus)
}
producer_surplus_actual <-
  compute_producer_surplus(
    p = p_actual, 
    marginal_cost = marginal_cost_estimate, 
    x = x, 
    M = M, 
    V = V, 
    beta = beta, 
    sigma = sigma, 
    mu = mu, 
    omega = omega
  )
summary(producer_surplus_actual)
producer_surplus_counterfactual <-
  compute_producer_surplus(
    p = p_counterfactual, 
    marginal_cost = marginal_cost_estimate, 
    x = x, 
    M = M, 
    V = V, 
    beta = beta, 
    sigma = sigma, 
    mu = mu, 
    omega = omega
  )
summary(producer_surplus_counterfactual)

# Compute Consumer Surplus Change
M_counterfactual$producer_surplus_actual<-0
M_counterfactual[M_counterfactual$j>0,]$producer_surplus_actual<-producer_surplus_actual

M_counterfactual$producer_surplus_counterfactual<-0
M_counterfactual[M_counterfactual$j>0,]$producer_surplus_counterfactual<-producer_surplus_counterfactual

change<-M_counterfactual[M_counterfactual$j>0,]%>%
  dplyr::mutate(producer_surplus_change=(producer_surplus_counterfactual-producer_surplus_actual)/producer_surplus_actual)%>%
  dplyr::group_by(j)%>%
  dplyr::summarise(producer_surplus_change=mean(producer_surplus_change))%>%
  dplyr::ungroup()
head(change,n=10) 

# Consumer Surplus Function
compute_consumer_surplus<-function(p, x, M, V, beta,sigma, mu, omega){
  # Part I: construct df based on the p
  M[M$j != 0, "p"] <- p
  df <- expand.grid(t = 1:T, i = 1:N, j = 0:J)
  df <- left_join(df, M, by = c("j", "t"))
  df <- left_join(df, V, by = c("i", "t"))
  df <- left_join(df, x, by = "j")
  df <- df %>% dplyr::filter(!is.na(p)) %>% dplyr::arrange(t,i,j)
  # Part II: compute indirect utility and alpha_i
  df$u <- compute_indirect_utility(df, beta, sigma, mu, omega)
  df$alpha_i <-  - exp(mu + omega * df$v_p)
  # Part III: compute the consumer surplus
  consumer_surplus <- df %>% dplyr::group_by(t, i) %>%
    dplyr::mutate(consumer_surplus = log(sum(exp(u))) / abs(alpha_i)) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(t, i, .keep_all = TRUE)
  return(consumer_surplus$consumer_surplus)
}

consumer_surplus_actual <- 
  compute_consumer_surplus(
    p = p_actual, 
    x = x, 
    M = M, 
    V = V, 
    beta = beta, 
    sigma = sigma, 
    mu = mu, 
    omega = omega
  )
consumer_surplus_counterfactual <- 
  compute_consumer_surplus(
    p = p_counterfactual, 
    x = x, 
    M = M, 
    V = V, 
    beta = beta, 
    sigma = sigma, 
    mu = mu, 
    omega = omega
  )
consumer_surplus_change <- 
  (sum(consumer_surplus_counterfactual) - 
     sum(consumer_surplus_actual)) /
  sum(consumer_surplus_actual)
consumer_surplus_change

