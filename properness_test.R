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
    # use estimated censoring distribution
    out[lhs] = pweibull(tstar, pred_shape, pred_scale, FALSE)^2 / pmax(eps, interp_surv(cens_fit, new_times = t[lhs]))
    out[rhs] = pweibull(tstar, pred_shape, pred_scale)^2 / pmax(eps, interp_surv(cens_fit, new_times = tstar))
  }

  mean(pmax(out, eps))
}

isbs = function(pred_shape, pred_scale, t, delta, cens_shape, cens_scale, cens_fit = NULL, eps = 1e-5) {
  # create 50 times between 5% and 95% quantile of the given (observed) times
  n = 50
  q5 = quantile(t, 0.05)
  q95 = quantile(t, 0.95)
  times = seq(q5, q95, length.out = n) # time points to evaluate SBS(t*)

  # get SBS(t*)
  scores = vapply(times, function(tstar) {
    sbs(pred_shape, pred_scale, t, delta, tstar, cens_shape, cens_scale, cens_fit = cens_fit, eps = eps)
  }, numeric(1))

  # calculate ISBS via trapezoidal integration rule + normalize for time range
  sum(diff(times) * (scores[-n] + scores[-1]) / 2) / (max(times) - min(times))
}

RCLL = function(pred_shape, pred_scale, t, delta, cens_shape, cens_scale, cens_fit = NULL, proper = FALSE, eps = 1e-5) {

  out = numeric(length(t))

  out[delta] = dweibull(t[delta], pred_shape, pred_scale)
  if (proper) {
    # divide by survival at outcome time (censoring distr)
    if (is.null(cens_fit)) {
      out[delta] = out[delta] / pmax(eps, pweibull(t[delta], cens_shape, cens_scale, FALSE))
    } else {
      out[delta] = out[delta] / pmax(eps, interp_surv(cens_fit, new_times = t[delta]))
    }
  }

  out[!delta] = pweibull(t[!delta], pred_shape, pred_scale, FALSE)

  if (proper) {
    # divide by density at outcome time (censoring distr)
    if (is.null(cens_fit)) {
      out[!delta] = out[!delta] / pmax(eps, dweibull(t[!delta], cens_shape, cens_scale))
    } else {
      out[!delta] = out[!delta] / pmax(eps, interp_pdf(cens_fit, new_times = t[!delta]))
    }
  }

  mean(-log(pmax(out, eps)))
}

# HELPER EXPERIMENT FUNCTIONS ----
tv_distance = function(P, Q) {
  # Create empirical cumulative distribution functions (ECDFs)
  ecdf_P = ecdf(P)
  ecdf_Q = ecdf(Q)

  # Define a grid over the range of both input vectors
  x_range = seq(min(c(P, Q)), max(c(P, Q)), length.out = 500)

  # Compute the absolute differences between the ECDFs
  differences = abs(ecdf_P(x_range) - ecdf_Q(x_range))

  # Return the Total Variation distance (max difference in CDFs)
  max(differences)
}

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

  #num_cores = detectCores() - 1
  num_cores = 60

  x = mclapply(seq.int(num_distrs), function(i) {
    # True event & censoring times
    true_y = rweibull(num_samples, surv_shape, surv_scale)
    true_c = rweibull(num_samples, cens_shape, cens_scale)

    # observed outcomes (time, status/delta)
    obs_t = pmin(true_y, true_c)
    obs_d = true_y == obs_t

    # calculate proportion of censoring
    prop_cens = sum(!obs_d) / num_samples

    # time to evaluate at (median) for SBS
    tau_median = median(obs_t)
    tau_10 = unname(quantile(obs_t, 0.1))
    tau_90 = unname(quantile(obs_t, 0.9))

    fit = NULL
    if (estimate_cens) {
      # use Kaplan-Meier to estimate G(t) censoring distribution
      fit = survival::survfit(survival::Surv(obs_t, 1 - obs_d) ~ 1)
    }

    c(
      # SBS at median observed time
      sbs(surv_shape, surv_scale, obs_t, obs_d, tau_median, cens_shape, cens_scale, fit),
      sbs(pred_shape, pred_scale, obs_t, obs_d, tau_median, cens_shape, cens_scale, fit),
      # SBS at 10% quantile observed time
      sbs(surv_shape, surv_scale, obs_t, obs_d, tau_10, cens_shape, cens_scale, fit),
      sbs(pred_shape, pred_scale, obs_t, obs_d, tau_10, cens_shape, cens_scale, fit),
      # SBS at 90% quantile observed time
      sbs(surv_shape, surv_scale, obs_t, obs_d, tau_90, cens_shape, cens_scale, fit),
      sbs(pred_shape, pred_scale, obs_t, obs_d, tau_90, cens_shape, cens_scale, fit),
      # RCLL
      RCLL(surv_shape, surv_scale, obs_t, obs_d, cens_shape, cens_scale, fit, FALSE),
      RCLL(pred_shape, pred_scale, obs_t, obs_d, cens_shape, cens_scale, fit, FALSE),
      # rRCLL
      RCLL(surv_shape, surv_scale, obs_t, obs_d, cens_shape, cens_scale, fit, TRUE),
      RCLL(pred_shape, pred_scale, obs_t, obs_d, cens_shape, cens_scale, fit, TRUE),
      # ISBS
      isbs(surv_shape, surv_scale, obs_t, obs_d, cens_shape, cens_scale, fit),
      isbs(pred_shape, pred_scale, obs_t, obs_d, cens_shape, cens_scale, fit),
      # censoring proportion
      prop_cens
    )
  }, mc.cores = num_cores)

  x = do.call(cbind, x)

  # average over all simulations
  means = rowMeans(x)

  list(
    SBS_median = means[1] - means[2],
    SBS_q10 = means[3] - means[4],
    SBS_q90 = means[5] - means[6],
    RCLL = means[7] - means[8],
    rRCLL = means[9] - means[10],
    ISBS = means[11] - means[12],
    prop_cens = means[13],
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

      tibble::tibble(
        # simulation number
        sim = i,
        # number of samples
        n = num_samples,
        # Y
        surv_shape = surv_shape,
        surv_scale = surv_scale,
        # C
        cens_shape = cens_shape,
        cens_scale = cens_scale,
        # S
        pred_shape = pred_shape,
        pred_scale = pred_scale,
        # mean proportion of censoring
        prop_cens = result$prop_cens,
        # total variational distance between Y and Y_hat
        tv_dist = result$tv_dist,
        # score mean differences (score(Y) - score(Y_hat))
        SBS_median_diff = result$SBS_median,
        SBS_q10_diff = result$SBS_q10,
        SBS_q90_diff = result$SBS_q90,
        ISBS_diff = result$ISBS,
        RCLL_diff = result$RCLL,
        rRCLL_diff = result$rRCLL,
        rSBS_median_diff = result$rSBS_median
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
