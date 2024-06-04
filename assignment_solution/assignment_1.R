library(evd)
library(dplyr)
##### Simulate Data
df <- data.frame(i = rep(1:1000, each = 2),k = rep(1:2, 1000),x = rep(0:1, 1000))

set.seed(1)
df$e <- evd::rgumbel(n = nrow(df))

beta <- 0.2
b    <- 1
df$latent <- beta*df$x + df$e

df <- df %>%
  group_by(i) %>%
  mutate(y = as.numeric(latent == max(latent))) %>%
  ungroup()

df <- df %>%
  group_by(i) %>%
  mutate(prob = exp(b * x) / sum(exp(b * x))) 

df <- subset(df,y==1)
loglikehihood <- sum(log(df$prob))

##### Estimate Parameter #####
compute_loglikehihood_a1 <- function(b,df){
  df <- df %>%
    group_by(i) %>%
    mutate(prob = exp(b * x) / sum(exp(b * x))) 
  df <- subset(df,y==1)
  liklihood <- -sum(log(df$prob)) / 1000
  return(liklihood)
}
compute_loglikehihood_a1(1,df)

# Plot for beta from 0 to 1 using ggplot2 scatter
library(ggplot2)
beta <- seq(0,1,0.1)
loglikehihood <- sapply(beta,compute_loglikehihood_a1,df = df)
df1 <- data.frame(beta = beta,loglikehihood = loglikehihood)
ggplot(df1,aes(x = beta,y = loglikehihood)) + geom_point()

# Find the maximum likelihood
result <- optim(par = 0,  # Initial guess for beta
                fn = compute_loglikehihood_a1,
                df = df,
                method = "Brent",
                lower = -1,
                upper = 1)
result$par
result$value
result$convergence
result$message

