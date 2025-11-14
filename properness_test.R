#' Testing properness of scoring rules
#' Execute: `Rscript properness_test.R 20 1000 50 FALSE 2434234`
#' CMD arguments are: `n_sims n_distrs n_samples cens_estimate seed`

library(parallel)
library(tibble)
library(progressr)
source("helper.R")

# SCORES ----
sbs = function(pred_shape, pred_scale, t, delta, tstar, cens_shape, cens_scale, cens_fit = NULL, eps = 1e-5) {
  out = numeric(length(t))

  lhs = delta == 1 & t <= tstar
  rhs = t > tstar

  if (is.null(cens_fit)) {
    # use true censoring distribution
    out[lhs] = pweibull(tstar, pred_shape, pred_scale, FALSE)^2 / pmax(eps, pweibull(t[lhs], cens_shape, cens_scale, FALSE))
    out[rhs] = pweibull(tstar, pred_shape, pred_scale)^2 / pmax(eps, pweibull(tstar, cens_shape, cens_scale, FALSE))
  } else {
    # use estimated censoring distribution (and constant interpolation)
    out[lhs] = pweibull(tstar, pred_shape, pred_scale, FALSE)^2 / pmax(eps, interp_surv(cens_fit, new_times = t[lhs], inter_type = "constant"))
    out[rhs] = pweibull(tstar, pred_shape, pred_scale)^2 / pmax(eps, interp_surv(cens_fit, new_times = tstar, inter_type = "constant"))
  }

  mean(pmax(out, eps))
}

isbs = function(pred_shape, pred_scale, t, delta, cens_shape, cens_scale, cens_fit = NULL, eps = 1e-5) {
  # create 50 times between 5% and 80% quantile of the given (observed) times
  n = 50
  q5 = quantile(t, 0.05)
  q80 = quantile(t, 0.80)
  times = seq(q5, q80, length.out = n) # time points to evaluate SBS(t*)

  # get SBS(t*)
  scores = vapply(times, function(tstar) {
    sbs(pred_shape, pred_scale, t, delta, tstar, cens_shape, cens_scale, cens_fit = cens_fit, eps = eps)
  }, numeric(1))

  # calculate ISBS via trapezoidal integration rule + normalize for time range
  sum(diff(times) * (scores[-n] + scores[-1]) / 2) / (max(times) - min(times))
}

RCLL = function(pred_shape, pred_scale, t, delta, cens_shape, cens_scale, cens_fit = NULL, weighted = FALSE, eps = 1e-5) {
  out = numeric(length(t))

  out[delta] = dweibull(t[delta], pred_shape, pred_scale)
  if (weighted) {
    # divide by survival at outcome time (censoring distr, linear interpolation)
    if (is.null(cens_fit)) {
      out[delta] = out[delta] / pmax(eps, pweibull(t[delta], cens_shape, cens_scale, FALSE))
    } else {
      out[delta] = out[delta] / pmax(eps, interp_surv(cens_fit, new_times = t[delta]))
    }
  }

  out[!delta] = pweibull(t[!delta], pred_shape, pred_scale, FALSE)

  if (weighted) {
    # divide by density at outcome time (censoring distr, linear interpolation)
    if (is.null(cens_fit)) {
      out[!delta] = out[!delta] / pmax(eps, dweibull(t[!delta], cens_shape, cens_scale))
    } else {
      out[!delta] = out[!delta] / pmax(eps, interp_pdf(cens_fit, new_times = t[!delta]))
    }
  }

  mean(-log(pmax(out, eps)))
}

# HELPER EXPERIMENT FUNCTIONS ----
tv_distance_weibull = function(shape1, scale1, shape2, scale2, n_points = 500) {
  # Define a grid over which to compute the PDFs
  x_range = seq(0, 3 * max(c(scale1, scale2)), length.out = n_points)

  # Compute the PDFs of the two distributions
  pdf1 = dweibull(x_range, shape = shape1, scale = scale1)
  pdf2 = dweibull(x_range, shape = shape2, scale = scale2)

  # Filter out NA or Inf values
  valid_indices = which(!is.na(pdf1) & !is.infinite(pdf1) & !is.na(pdf2) & !is.infinite(pdf2))
  x_range = x_range[valid_indices]
  pdf1 = pdf1[valid_indices]
  pdf2 = pdf2[valid_indices]

  # Compute the absolute differences between the PDFs
  diffs = abs(pdf1 - pdf2)

  # Numerically integrate the PDF differences over the range of x (using the
  # trapezoidal rule) and divide by 2 to get the normalized TV distance in [0,1]
  dx = diff(x_range)
  tv_distance = 0.5 * sum((diffs[-length(diffs)] + diffs[-1]) / 2 * dx)

  tv_distance
}

run = function(surv_shape, cens_shape, pred_shape,
               surv_scale, cens_scale, pred_scale,
               num_distrs, num_samples, estimate_cens) {
  # --- Change Number of cores ---
  num_cores = 1

  x = mclapply(seq.int(num_distrs), function(i) {
    # --- Simulate data ---
    # True event & censoring times
    true_y = rweibull(num_samples, surv_shape, surv_scale)
    true_c = rweibull(num_samples, cens_shape, cens_scale)

    # observed outcomes (time, status/delta)
    obs_t = pmin(true_y, true_c)
    obs_d = true_y == obs_t

    # --- Analytic mean survival times ---
    # True and predicted mean survival times for Weibull distributions
    m_true = surv_scale * gamma(1 + 1 / surv_shape)
    m_pred = pred_scale * gamma(1 + 1 / pred_shape)

    # --- Trapezoidal rule mean estimation (finite sample approx.) ---
    t_sorted = sort(obs_t)
    S_true_t = pweibull(t_sorted, shape = surv_shape, scale = surv_scale, lower.tail = FALSE)
    S_pred_t = pweibull(t_sorted, shape = pred_shape, scale = pred_scale, lower.tail = FALSE)
    m_true_est = trap_int(S_true_t, t_sorted)
    m_pred_est = trap_int(S_pred_t, t_sorted)

    # --- Derived diagnostics ---
    prop_cens = mean(!obs_d) # proportion of censoring
    n_events = sum(obs_d)
    max_t = max(obs_t)

    # --- Evaluation time points for SBS ---
    tau_median = median(obs_t)
    tau_10 = unname(quantile(obs_t, 0.1))
    tau_90 = unname(quantile(obs_t, 0.9))

    # --- Diagnostics for SBS ----
    # True and censoring survival at each tau
    S_Y_q10 = pweibull(tau_10, shape = surv_shape, scale = surv_scale, lower.tail = FALSE)
    S_Y_med = pweibull(tau_median, shape = surv_shape, scale = surv_scale, lower.tail = FALSE)
    S_Y_q90 = pweibull(tau_90, shape = surv_shape, scale = surv_scale, lower.tail = FALSE)

    S_C_q10 = pweibull(tau_10, shape = cens_shape, scale = cens_scale, lower.tail = FALSE)
    S_C_med = pweibull(tau_median, shape = cens_shape, scale = cens_scale, lower.tail = FALSE)
    S_C_q90 = pweibull(tau_90, shape = cens_shape, scale = cens_scale, lower.tail = FALSE)

    # tail terms for (Y,C)
    S_Y_tail = pweibull(max_t, shape = surv_shape, scale = surv_scale, lower.tail = FALSE)
    S_C_tail = pweibull(max_t, shape = cens_shape, scale = cens_scale, lower.tail = FALSE)
    epsilon = S_Y_tail * S_C_tail

    # SBS shift
    shift_q10 = shift_fun(S_Y_q10, S_C_q10, epsilon)
    shift_med = shift_fun(S_Y_med, S_C_med, epsilon)
    shift_q90 = shift_fun(S_Y_q90, S_C_q90, epsilon)

    fit = NULL
    if (estimate_cens) {
      # use Kaplan-Meier to estimate G(t) censoring distribution
      fit = survival::survfit(survival::Surv(obs_t, 1 - obs_d) ~ 1)
    }

    c(
      # SBS at median observed time
      sbs(surv_shape, surv_scale, obs_t, obs_d, tau_median, cens_shape, cens_scale, fit), # 1
      sbs(pred_shape, pred_scale, obs_t, obs_d, tau_median, cens_shape, cens_scale, fit), # 2
      # SBS at 10% quantile observed time
      sbs(surv_shape, surv_scale, obs_t, obs_d, tau_10, cens_shape, cens_scale, fit), # 3
      sbs(pred_shape, pred_scale, obs_t, obs_d, tau_10, cens_shape, cens_scale, fit), # 4
      # SBS at 90% quantile observed time
      sbs(surv_shape, surv_scale, obs_t, obs_d, tau_90, cens_shape, cens_scale, fit), # 5
      sbs(pred_shape, pred_scale, obs_t, obs_d, tau_90, cens_shape, cens_scale, fit), # 6
      # RCLL
      RCLL(surv_shape, surv_scale, obs_t, obs_d, cens_shape, cens_scale, fit, FALSE), # 7
      RCLL(pred_shape, pred_scale, obs_t, obs_d, cens_shape, cens_scale, fit, FALSE), # 8
      # wRCLL (weighted RCLL)
      RCLL(surv_shape, surv_scale, obs_t, obs_d, cens_shape, cens_scale, fit, TRUE), # 9
      RCLL(pred_shape, pred_scale, obs_t, obs_d, cens_shape, cens_scale, fit, TRUE), # 10
      # ISBS
      isbs(surv_shape, surv_scale, obs_t, obs_d, cens_shape, cens_scale, fit), # 11
      isbs(pred_shape, pred_scale, obs_t, obs_d, cens_shape, cens_scale, fit), # 12
      # Diagnostics
      prop_cens, n_events, max_t, m_true, m_pred, m_true_est, m_pred_est, # 13-19
      S_Y_q10, S_C_q10, S_Y_med, S_C_med, S_Y_q90, S_C_q90, # 20-25
      S_Y_tail, S_C_tail, epsilon, shift_q10, shift_med, shift_q90 # 26-31
    )
  }, mc.cores = num_cores)

  x = do.call(cbind, x)

  # average over all simulations
  means = rowMeans(x)

  # Differences are always: (True Distribution) - (Predicted Distribution)
  list(
    SBS_median_diff = means[1] - means[2], # same as mean(x[1,] - x[2,])
    SBS_median_sd = sd(x[1,] - x[2,]),
    SBS_q10_diff = means[3] - means[4],
    SBS_q10_sd = sd(x[3,] - x[4,]),
    SBS_q90_diff = means[5] - means[6],
    SBS_q90_sd = sd(x[5,] - x[6,]),
    RCLL_diff = means[7] - means[8],
    RCLL_sd = sd(x[7,] - x[8,]),
    wRCLL_diff = means[9] - means[10],
    wRCLL_sd = sd(x[9,] - x[10,]),
    ISBS_diff = means[11] - means[12],
    ISBS_sd = sd(x[11,] - x[12,]),
    prop_cens = means[13],
    n_events = as.integer(means[14]),
    max_t = means[15],
    m_true = means[16],
    m_pred = means[17],
    m_true_est = means[18],
    m_pred_est = means[19],
    S_Y_q10 = means[20],
    S_C_q10 = means[21],
    S_Y_med = means[22],
    S_C_med = means[23],
    S_Y_q90 = means[24],
    S_C_q90 = means[25],
    S_Y_tail = means[26],
    S_C_tail = means[27],
    epsilon  = means[28],
    shift_q10 = means[29],
    shift_med = means[30],
    shift_q90 = means[31],
    tv_dist = tv_distance_weibull(
      shape1 = surv_shape, scale1 = surv_scale,
      shape2 = pred_shape, scale2 = pred_scale
    )
  )
}

run_experiment = function(num_sims = 20, num_distrs = 1000, num_samples = 1000, estimate_cens = FALSE, seed = NULL) {
  cat(sprintf("#Simulations: %i, #Distributions: %i, #Samples: %i, Estimate Censoring = %s\n\n",
              num_sims, num_distrs, num_samples, estimate_cens))

  # set global seed
  set.seed(seed)

  lower = 0.5
  upper = 5

  simulations = seq.int(num_sims)
  with_progress({
    p = progressor(along = simulations)

    res = lapply(simulations, function(i) {
      p(sprintf("Simulation %i", i))

      # Specify parameters for {Y,C,S}
      surv_shape = runif(1, lower, upper)
      cens_shape = runif(1, lower, upper)
      pred_shape = runif(1, lower, upper)

      cens_scale = runif(1, lower, upper)
      pred_scale = runif(1, lower, upper)
      surv_scale = runif(1, lower, upper)

      result = run(
        surv_shape, cens_shape, pred_shape,
        surv_scale, cens_scale, pred_scale,
        num_distrs, num_samples, estimate_cens
      )

      tibble::as_tibble(
        c(
          list(
            sim = i, # simulation number
            n = num_samples, # number of samples
            surv_shape = surv_shape, # Y
            surv_scale = surv_scale,
            cens_shape = cens_shape, # C
            cens_scale = cens_scale,
            pred_shape = pred_shape, # Y_hat
            pred_scale = pred_scale
          ),
          result # automatically adds all elements of `result`
        )
      )
    })
  })

  dplyr::bind_rows(res)
}

# RUN EXPERIMENT ----
options(progressr.enable = TRUE)
handlers(global = TRUE)
handlers("progress")

# Get command-line arguments
args = commandArgs(trailingOnly = TRUE)

# Parse arguments
num_sims = as.integer(args[1])
num_distrs = as.integer(args[2])
num_samples = as.integer(args[3])
estimate_cens = as.logical(args[4]) # Accepts "TRUE" or "FALSE"
seed = as.integer(args[5])

#' `n_sims` = 10000 # number of independent simulations (different distribution choices for {Y,C,S})
#' `n_distrs` = 1000 # number of random sampled distributions
#' `n_samp` = 50 # how many samples to draw from the distributions
#' `estimate_cens` = FALSE # whether to use an estimated censoring distribution (via Kaplan-Meier) or the true {C} in scores
#' `seed` = 20240402 # seed for reproducibility

res = run_experiment(num_sims, num_distrs, num_samples, estimate_cens, seed = seed)

# SAVE RESULTS ----
saveRDS(res, file = sprintf("res_sims%s_distrs%s_n%s_%i.rds", num_sims, num_distrs, num_samples, estimate_cens))
