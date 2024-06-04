library(dplyr)
library(evd)
library(purrr)
options(scipen = 20)
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

compute_derivative_share_smooth <- function(df, beta, sigma, mu, omega){
  # Part I: Compute Choice
  df <- compute_choice_smooth(df, beta, sigma, mu, omega)
  df <- df %>% filter(j != 0)
  df$alpha_i <- mu + omega * df$v_p
  
  # Part II: For each market, compute the derivative of share with respect to price, put results in a list
  T <- length(unique(df$t))
  result <- list(length = T)
  for (t in 1:T) {
    print(t)
    data_j <- df %>% filter(t == t)
    J_t <- length(unique(data_j$j))
    matrix <- matrix(0, nrow = J_t, ncol = J_t)
    j_list <- unique(data_j$j)
    for (j in 1:J_t) {
      for (k in 1:J_t) {
        if (j == k) {
          data_jk <- data_j %>% filter(j == j_list[j])
          data_jk$obj <- data_jk$alpha_i * data_jk$q * (1 - data_jk$q)
          matrix[j, k] <- mean(data_jk$obj)
        } else {
          data_jk <- data_j %>% filter(j == j_list[j] | j == j_list[k])
          data_jk <- data_jk %>%
            arrange(i, j) %>%
            group_by(i) %>%
            mutate(lag_q = lag(q)) %>%
            ungroup() %>%
            filter(!is.na(lag_q))
          data_jk$obj <- data_jk$alpha_i * data_jk$q * data_jk$lag_q
          matrix[j, k] <- -mean(data_jk$obj)
        }
      }
    }
    result[[t]] <- matrix
  }
  return(result)
}
result <- compute_derivative_share_smooth(df, beta, sigma, mu, omega)
