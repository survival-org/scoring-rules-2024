library(mlr3)
library(mlr3proba)

#' Simulate Log-Normal Survival Data with Independent and Administrative Censoring
#'
#' @description
#' Simulates right-censored survival data under a log-normal AFT model with:
#'   - Continuous covariate x1 ~ N(0,1)
#'   - Binary covariate x2 ~ Bernoulli(0.5)
#'   - Optional interaction term x1:x2
#'   - Group-specific log-scale parameter sigma depending on x2
#'   - Independent exponential censoring
#'   - Additional administrative censoring at `admin_time`
#'
#' True model:
#'   log(T_i) = eta_i + sigma_i * Z_i, with Z_i ~ N(0,1)
#'
#' Therefore T_i ~ LogNormal(meanlog = eta_i, sdlog = sigma_i)
#'
#' @param n Integer. Sample size.
#' @param b0,b1,b2,b12 Numeric. Regression coefficients. `b12` is x1:x2 interaction
#' @param sigma Numeric vector of length 2. Scale parameters for x2 = 0 and x2 = 1.
#' @param lambdaC Numeric. Exponential censoring rate.
#' @param times Numeric vector. Time grid for survival matrix output.
#' @param admin_time Numeric. Administrative censoring time.
#'
#' @return List with elements:
#'   - `task`: mlr3proba survival `TaskSurv` object
#'   - `pred`: mlr3proba `PredictionSurv` object containing true survival distribution
#' (n x length(times)) and sign-reverse linear predictor
#'   - `data`: Simulated data.frame
#' @noRd
simulate = function(
    n = 500,
    b0 = 1.0,
    b1 = -0.3,
    b2 = 0.2,
    b12 = 0.5,
    sigma = c(0.5, 0.5),
    lambdaC = 0.08,
    times = seq(0, 10, by = 0.1),
    admin_time = 10
) {
  # ---------------------------
  # 1. Generate covariates
  # ---------------------------
  df = data.frame(
    id = seq_len(n),
    x1 = rnorm(n),
    x2 = rbinom(n, 1, 0.5)
  )

  # Group-specific sigma (depends on x2)
  stopifnot(length(sigma) == 2)
  df$sigma_i = ifelse(df$x2 == 1, sigma[2], sigma[1])

  # Linear predictor (AFT scale: larger eta => larger T => lower risk)
  df$eta = b0 + b1 * df$x1 + b2 * df$x2 + b12 * df$x1 * df$x2

  # -----------------------------------
  # 2. Generate true survival times
  # -----------------------------------
  # log(Y) = eta + sigma * Z
  df$true_y = rlnorm(n, meanlog = df$eta, sdlog = df$sigma_i)
  # df$true_y = exp(df$eta + df$sigma_i * rnorm(n)) # same, more manual

  # -----------------------------
  # 3. Independent censoring
  # -----------------------------
  df$true_c = rexp(n, rate = lambdaC)

  # Administrative censoring
  df$true_c_admin = pmin(df$true_c, admin_time)

  # Observed time and status
  df$time = pmin(df$true_y, df$true_c_admin)
  df$status = as.integer(df$true_y <= df$true_c_admin)

  # -----------------------------------------------------------------
  # 4. Observation-wise and Total Right-Censored True Likelihood
  # -----------------------------------------------------------------
  # For delta = 1 (event):  -log[f(t|x)]
  # For delta = 0 (censored):  -log[S(t|x)]
  # The true RCLL value here depends on the observed times so it's not affected
  # by the choice of the prediction time grid
  log_f = dlnorm(df$time, meanlog = df$eta, sdlog = df$sigma_i, log = TRUE)
  log_S = plnorm(df$time, meanlog = df$eta, sdlog = df$sigma_i, lower.tail = FALSE, log.p = TRUE)
  df$rcll_i = - (df$status * log_f + (1 - df$status) * log_S)
  rcll = mean(df$rcll_i)

  # ---------------------------
  # 4. True survival matrix
  # ---------------------------
  # S(t|x) = 1 - Phi((log(t) - eta) / sigma)
  surv_mat = matrix(NA, nrow = n, ncol = length(times))
  log_times = log(pmax(times, 1e-8))

  z = sweep(
    outer(log_times, df$eta, "-"), # m x n
    2,
    df$sigma_i,
    "/"
  )
  surv_mat = 1 - pnorm(t(z)) # now n x m

  # numerical stability (was getting 1 - 1e-48 ~ 1, but not quite!)
  if (times[1] == 0) {
    surv_mat[, 1] = 1
  }
  # check S(t)
  survdistr::assert_prob(surv_mat, times = times)

  # create mlr3proba survival task
  task = mlr3proba::as_task_surv(
    x = df[, c("time", "status", "x1", "x2")],
    time = "time",
    event = "status",
    id = sprintf("sim_lognormal_%s_%04d", format(Sys.time(), "%Y%m%d%H%M%S"), sample.int(9999, 1))
  )

  # ---------------------------
  # 6. Create PredictionSurv
  # ---------------------------
  pred_list = mlr3proba::surv_return(
    times = times,
    surv = surv_mat,
    lp = -df$eta # negative for risk ordering
  )

  pred = mlr3proba::PredictionSurv$new(
    task = task,
    crank = pred_list$crank,
    distr = pred_list$distr,
    lp = pred_list$lp
  )

  list(
    task = task,
    pred = pred,
    data = df,
    rcll = rcll
  )
}
