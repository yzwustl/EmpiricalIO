options(scipen = 200)
library(dplyr)
library(vtable)
library(np)
library(ggplot2)
###### Setting Parameters ###### 
beta_0 = 1
beta_l = 0.2
beta_k = 0.7
alpha  = 0.7
sigma_eta = 0.2
sigma_v = 0.5
sigma_w = 0.1
delta   = 0.05
j = 1000
t = 10

###### Production and Optimal Labor ######
log_production <- function(l,k,omega,eta,beta_0,beta_l,beta_k){
  log_production <- beta_0 + beta_l * l + beta_k * k + omega + eta
  return(log_production)
}

log_labor_choice <- function(k, wage, omega, beta_0, beta_l, beta_k, sigma_eta){
  log_labor_choice <- (log(wage/beta_l) - beta_0 - beta_k * k - omega - 0.5 * sigma_eta^2) / (beta_l-1)
  return(log_labor_choice)
}

log_labor_choice_error <- function(k, wage, omega, beta_0, beta_l, beta_k, iota, sigma_eta){
  log_labor_choice <- (log(wage/beta_l) - beta_0 - beta_k * k - omega - iota - 0.5 * sigma_eta^2) / (beta_l-1)
  return(log_labor_choice)
}

##### Investment Process #####
gamma = 0.1

investment_choice <- function(k,omega,gamma,delta){
  investment_choice <- (delta + gamma*omega) * exp(k)
  return(investment_choice)
}

##### Simulation for Period 1#####
df <- data.frame(j = 1:1000,t = 1)
df$wage= 0.5

set.seed(1)
df$k <- rnorm(j,1,0.5)
df$omega = rnorm(j, 0,  sigma_v / sqrt(1 - alpha^2))

df$iota = rnorm(j, 0, 0.05)
df$l <- log_labor_choice(df$k, df$wage, df$omega, beta_0, beta_l, beta_k, sigma_eta)
df$l_error <- log_labor_choice_error(df$k, df$wage, df$omega, beta_0, beta_l, beta_k, df$iota, sigma_eta)
df$inv <- investment_choice(df$k,df$omega,gamma,delta)

df$eta = rnorm(j, 0, sigma_eta)
df$y = log_production(df$l,df$k,df$omega,df$eta,beta_0,beta_l,beta_k)
df$y_error = log_production(df$l_error,df$k,df$omega,df$eta,beta_0,beta_l,beta_k) 

##### Simulate for Period 2 to 10 #####
df_all <- df
for (t in 2:10) {
  # change time index
  df$t <- t
  
  # draw wage
  df$wage <- 0.5
  
  # update capital
  df <- 
    df %>%
    dplyr::mutate(
      k = log((1 - delta) * exp(k) + inv)
    )
  
  # update omega
  df <- 
    df %>%
    dplyr::mutate(
      v = rnorm(j, 0, sigma_v),
      omega = alpha * omega + v
    )
  
  # compute labor and investment
  df <- 
    df %>%
    dplyr::mutate(iota = rnorm(j, 0, 0.05)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      l = log_labor_choice(k, wage, omega, beta_0, beta_l, beta_k, sigma_eta),
      l_error = log_labor_choice_error(k, wage, omega, beta_0, beta_l, beta_k, iota, sigma_eta),
      inv = investment_choice(k, omega, gamma, delta)
    ) %>%
    dplyr::ungroup()
  # compute output
  df <- 
    df %>%
    dplyr::mutate(eta = rnorm(j, 0, sigma_eta)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      y = log_production(l, k, omega, eta, beta_0, beta_l, beta_k),
      y_error = log_production(l_error, k, omega, eta, beta_0, beta_l, beta_k) 
    ) %>%
    dplyr::ungroup()
  
  # append
  df_all <- dplyr::bind_rows(df_all, df)
}

####### Summary Statistics Table of All variables
sumtable(data=df_all,digits = 3)

####### Estimation ##########
### OLS
ols <- lm(y_error ~ l_error + k, data = df_all)
summary(ols)

### Olley and Pakes First Stage
# Estimate the optimal bandwidth using npplregbw
 optimal_bw <- npplregbw(
   data = df_all,
   formula = y_error ~ l_error | k + inv
)

# Hard Code Optimal Bandwidth
bandwidth <- matrix(c(0.07355058,0.01435558, 0.03908594, 0.01191551),nrow = 2,ncol = 2,byrow = TRUE)

result_1st <-
  npplreg(
    data = df_all,
    formula = y_error ~ l_error | k + inv,
    bws = optimal_bw
  )

summary(result_1st)

result_1st_plot <-
  data.frame(
    actual = df_all$y_error,
    fitted = fitted(result_1st)
  )

result_1st_plot %>%
  ggplot(
    aes(
      x = fitted, 
      y = actual
    )
  ) +
  geom_point()

#### Second Stage Estimation ####
# choose two columns in df_all
df_all_1st = df_all
df_all_1st$y_error_tilde = df_all$y_error - result_1st$xcoef["l_error"] * df_all$l_error
df_all_1st$phi = fitted(result_1st) - result_1st$xcoef["l_error"] * df_all$l_error

# create lag phi by j using group by 
df_all_1st <- df_all_1st %>%
  group_by(j) %>%
  mutate(phi_t_1 = dplyr::lag(phi))

# Second Stage Moment Function
moment_olleypakes_2nd <- function(alpha,beta_0,beta_k,df_all_1st){
  df_all_1st <- df_all_1st %>%
    group_by(j) %>% mutate(k_t_1 = dplyr::lag(k,1),inv_t_1 = dplyr::lag(inv,1))
  df_all_1st$moment <- df_all_1st$y_error_tilde - beta_0 - beta_k*df_all_1st$k - alpha*(df_all_1st$phi_t_1 - beta_0 - beta_k * df_all_1st$k_t_1)
  df_all_1st$moment1 <- df_all_1st$moment * df_all_1st$k
  df_all_1st$moment2 <- df_all_1st$moment * df_all_1st$k_t_1
  df_all_1st$moment3 <- df_all_1st$moment * df_all_1st$inv_t_1

  moment <- cbind(df_all_1st$moment1,df_all_1st$moment2,df_all_1st$moment3)
  moment <- apply(moment, 2, mean, na.rm = TRUE)
  
  return(moment)
}

objective_olleypakes_2nd <- function(theta,df_all_1st,W){
  alpha  = theta[1]
  beta_0 = theta[2]
  beta_k = theta[3]
  
  moment <- t(moment_olleypakes_2nd(alpha,beta_0,beta_k,df_all_1st)) %*% W %*% moment_olleypakes_2nd(alpha,beta_0,beta_k,df_all_1st)
  
  return(moment)
}

##### Plot: Just varying with alpha #####
alpha_seq = seq(0,1,0.1)

objective_alpha <- c()

for (i in 1:length(alpha_seq)){
  theta = c(alpha[i],beta_0,beta_k)
  objective_alpha[i] = objective_olleypakes_2nd(theta,df_all_1st,diag(3))
}

# plot scatter using ggplot
df_plot <- data.frame(alpha = alpha,objective = objective_alpha)
ggplot(df_plot,aes(x = alpha,y = objective)) + geom_point()

##### GMM Estimation #####
theta <- c(alpha, beta_0, beta_k)
W <- diag(3)
result_2nd <-
  optim(
    par = theta,
    f = objective_olleypakes_2nd,
    method = "L-BFGS-B",
    df_all_1st = df_all_1st,
    W = W
  )
result_2nd
