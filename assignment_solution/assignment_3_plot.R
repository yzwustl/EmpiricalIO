compute_objective_range <- function(theta_true, param_idx, param_range, df_share, X, M, V, e) {
  objective_values <- numeric(length(param_range))
  
  for (i in seq_along(param_range)) {
    theta <- theta_true
    theta[param_idx] <- param_range[i]
    objective_values[i] <- compute_nlls_objective_a3(theta, df_share, X, M, V, e)
  }
  
  data.frame(param_value = param_range, objective_value = objective_values)
}

# Vectorize true parameters
theta_true <- c(beta, sigma, mu, omega)

# Define parameter names
param_names <- c("beta1", "beta2", "beta3", "sigma1", "sigma2", "sigma3", "mu", "omega")

# Compute objective function for each parameter using actual shocks
obj_actual <- list()
for (i in 1:length(theta_true)) {
  obj_actual[[i]] <- compute_objective_range(theta_true, i, seq(theta_true[i] - 0.5, theta_true[i] + 1.5, by = 0.1), df_share, df, M, V, e)
}

# Plot the objective function for each parameter using actual shocks
par(mfrow = c(2, 4))
for (i in 1:length(theta_true)) {
  plot(obj_actual[[i]]$param_value, obj_actual[[i]]$objective_value, type = "l", xlab = param_names[i], ylab = "Objective")
}

# Compute objective function for each parameter using Monte Carlo shocks
obj_mc <- list()
for (i in 1:length(theta_true)) {
  obj_mc[[i]] <- compute_objective_range(theta_true, i, seq(theta_true[i] - 0.5, theta_true[i] + 1.5, by = 0.1), df_share, df, M, V_mcmc, e_mcmc)
}

# Plot the objective function for each parameter using Monte Carlo shocks
par(mfrow = c(2, 4))
for (i in 1:length(theta_true)) {
  plot(obj_mc[[i]]$param_value, obj_mc[[i]]$objective_value, type = "l", xlab = param_names[i], ylab = "Objective (MC)")
}