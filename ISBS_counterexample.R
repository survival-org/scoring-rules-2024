library(ggplot2)

# Y_true ∼ Gamma(k,λ)
# Y_pred ~ Exp(λ)

# Y_true == Y_pred for k = 1

# ISBS expectation for tau* => Inf
cmp_integral_inf = function(k, lambda) {
  (k - 1.5 + 2^(1 - k)) / lambda
}

# tau => Inf, find counterexample
check_integral_values = function(k_values, lambda_values) {
  for (lambda in lambda_values) {
    I_k1 = cmp_integral_inf(k = 1, lambda)  # baseline at k=1

    for (k in k_values) {
      I_val = cmp_integral_inf(k, lambda)

      if (k != 1 && I_val < I_k1) {
        cat(sprintf("For lambda = %.3f, k = %.3f: I(inf) = %.5f < baseline (k=1) = %.5f\n",
                    lambda, k, I_val, I_k1))
      }
    }
  }
}

k_grid = seq(0.3, 2, by = 0.1)
lambda_grid = c(0.5, 1, 2)
check_integral_values(k_grid, lambda_grid)

# ISBS expectation for tau* finite: numerical approximation
cmp_integral_finite_tau = function(k, lambda, tau) {
  integrand = function(t) {
    # Gamma CDF and survival
    F_Y = pgamma(t, shape = k, rate = lambda, lower.tail = TRUE)
    S_Y = 1 - F_Y

    # Exponential CDF and survival
    F_Yhat = 1 - exp(-lambda * t)
    S_Yhat = 1 - F_Yhat

    S_Yhat^2 * F_Y + F_Yhat^2 * S_Y
  }

  # Use integrate with proper error handling
  res = integrate(integrand, lower = 0, upper = tau,
                  subdivisions = 1000, rel.tol = 1e-10, abs.tol = 1e-12)

  if (res$message != "OK") warning("Integration warning: ", res$message)

  res$value
}

# Check: large finite tau == infinite tau formula
k = 1.17
lambda = 1.3
tau = 1000  # large number to approximate infinity

cmp_integral_finite_tau(k, lambda, tau)
cmp_integral_inf(k, lambda) # same as above, good!

# find counterexample for finite tau*
check_integral_values_finite_tau = function(k_values, lambda_values, tau_values) {
 for (lambda in lambda_values) {
    for (tau in tau_values) {
      I_k1 = cmp_integral_finite_tau(k = 1, lambda, tau)  # baseline

      for (k in k_values) {
        if (k == 1) next
        I_val = cmp_integral_finite_tau(k, lambda, tau)

        if (I_val < I_k1) {
          cat(sprintf("lambda=%.3f, tau=%.3f, k=%.3f: I(tau) = %.6f < baseline(k=1) = %.6f\n",
                      lambda, tau, k, I_val, I_k1))
        }
      }
    }
  }
}

k_grid = seq(0.3, 2, by = 0.1)
lambda_grid = c(0.5, 1, 2)
tau_grid = c(1, 2, 5, 10, 100)
check_integral_values_finite_tau(k_grid, lambda_grid, tau_values = tau_grid)

# plot some gamma survival functions
plot_gamma_survival = function(k_values, lambda=1, t_max=5) {
  t_grid = seq(0, t_max, length.out=500)

  surv_data = do.call(rbind, lapply(k_values, function(k) {
    S_t = 1 - pgamma(t_grid, shape = k, rate = lambda)
    data.frame(t = t_grid, survival = S_t, k = factor(k))
  }))

  ggplot(surv_data, aes(x = t, y = survival, color = k)) +
    geom_line(linewidth = 1) +
    labs(title = paste0("Gamma Survival Functions with rate = ", lambda),
         x = "Time (t)", y = "S(t)",
         color = "Shape k") +
    theme_minimal(base_size = 16)
}

plot_gamma_survival(k_values = c(0.4, 0.6, 1, 1.1), lambda = 0.4, t_max = 15)
