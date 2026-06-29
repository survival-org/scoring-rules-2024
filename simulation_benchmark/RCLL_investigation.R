#' Investigation: Fixed test set size, varying prediction grid coarseness.
#' Grid defined by proportions (0.02, 0.05, 0.1, 0.2, 0.5, 0.8, 1.0) of unique event times.
#' How much does the density approximation affects the sensitivity of RCLL to model
#' misspecification.
#'
#' Run: `Rscript simulation_benchmark/RCLL_investigation.R`
suppressWarnings({
  library(mlr3) # v1.7.1
  library(mlr3proba) # v0.8.9
  library(mlr3extralearners) # v1.5.1
  library(mlr3pipelines) # v0.10.0
  library(survdistr) # v0.0.3
  library(data.table)
  library(future)
  library(future.apply)
  library(progressr)
})
source("simulation_benchmark/simulate.R")

# Parallelization + progress bars
plan("multisession", workers = 80)
options(progressr.enable = TRUE)
handlers(global = TRUE)
handlers("progress")
set.seed(42)

# Load tasks and learner IDs
tasks = readRDS("simulation_benchmark/tasks.rds")
learner_ids = readRDS("simulation_benchmark/learner_ids.rds")

# Monte Carlo repetitions
n_rsmps = 100
rsmp_id = seq_len(n_rsmps)

# Test set sizes (sampled from DGP for prediction)
n_tests = c(50, 100, 250, 500)

# Proportions of unique event times to use as prediction grid
proportions = c(0.02, 0.05, 0.1, 0.2, 0.5, 0.8, 1)

# Precompute deterministic time grids for each task and proportion
# (ensures same grid across replicates)
# time_grids = list()
# for (task_name in names(tasks)) {
#   train_task = tasks[[task_name]]
#   all_times = train_task$unique_event_times()
#   n_total = length(all_times)
#   for (prop in proportions) {
#     if (prop == 1) {
#       sampled_times = all_times  # full grid
#     } else {
#       n_times = ceiling(prop * n_total)
#       # Deterministic sampling: seed based on task and proportion
#       set.seed(1000 * (which(names(tasks) == task_name)) +
#                100 * (which(proportions == prop)))
#       sampled_times = sort(sample(all_times, size = n_times))
#     }
#     time_grids[[task_name]][[as.character(prop)]] = sampled_times
#   }
# }

# Simulation grid: task, replicate, proportion
sim_grid = data.table::CJ(
  task_id = names(tasks),
  rsmp_id = rsmp_id,
  prop = proportions,
  n_test = n_tests
)
sim_grid[, config_id := .I]

# Measures
# RCLL => right-censored negative log-loss
rcll = msr("surv.rcll")

# measure distance between predictions and true S(t)
# MISE => mean integrated squared error
mise = function(S_true, S_pred, times) {
  dt_vec = diff(times)
  sq_diff = (S_true - S_pred)^2
  mean(rowSums(0.5 * (sq_diff[, -ncol(sq_diff)] + sq_diff[, -1]) * dt_vec))
}

# alternative S(t) distance: MIAE => mean integrated absolute error
miae = function(S_true, S_pred, times) {
  dt_vec = diff(times)
  abs_diff = abs(S_true - S_pred)
  mean(rowSums(0.5 * (abs_diff[, -ncol(abs_diff)] + abs_diff[, -1]) * dt_vec))
}

# Eval function
eval_config = function(row, p) {
  config_id = row$config_id
  cens_id = row$task_id
  rsmp_id = row$rsmp_id
  prop = row$prop
  n_test = row$n_test

  p(sprintf("Config-id: %i, Task: %s, RSMP-id: %i, Prop: %f, n_test: %i",
            config_id, cens_id, rsmp_id, prop, n_test))

  train_task = tasks[[cens_id]]

  # Retrieve precomputed time grid for this task and proportion
  # times = time_grids[[cens_id]][[as.character(prop)]]

  # Seed depends on task, replicate, and proportion
  seed = 100000 * rsmp_id +
         1000 * as.integer(prop * 100) +
         10 * n_test +
         ifelse(cens_id == "low", 1, 2)
  set.seed(seed)

  # Sample prediction grid times (proportion of unique event times)
  uevents = train_task$unique_event_times()
  n_total = length(uevents)
  n_times = ceiling(prop * n_total)
  times = sort(sample(uevents, size = n_times, replace = FALSE))

  # Simulate test data (using the same dense grid for true survival)
  while (TRUE) {
    test_sim = simulate(
      n = n_test,
      b0 = 1.15, b1 = 0.15, b2 = -0.55, b12 = -0.75,
      sigma = c(0.5, 1.5),
      lambdaC = if (cens_id == "low") 0.075 else 0.45,
      times = uevents
    )
    test_task = test_sim$task
    cp = test_task$cens_prop()
    # guard against low sample size (n) and extreme censoring in test set
    # which can cause C-index to be undefined
    if (!is.null(cp) && length(cp) == 1 && !is.na(cp) && cp > 0 && cp < 0.9) break
  }

  true_pred = test_sim$pred
  s_true = true_pred$data$distr # True S(t) on the dense grid (for MISE/MIAE calculation)

  # interpolate true S(t) on the sparse grid for RCLL calculation
  true_pred$data$distr = survdistr::interp(
    x = true_pred$data$distr,
    times = uevents,
    eval_times = times, # times is a subset of uevents, so we just filter columns here
    method = "const_surv",
    output = "surv",
    add_times = TRUE,
    check = FALSE
  )

  # True model results
  true_model_res = suppressWarnings(
    data.table(
      learner_id = "true",
      mise = 0,
      miae = 0,
      rcll = true_pred$score(rcll),
      true_rcll = test_sim$rcll
    )
  )

  # Evaluate learner predictions
  out = list()
  for (learner_id in learner_ids) {
    if (learner_id == "RSF") {
      # load trained RSF model separately (as it is large in size)
      learner = readRDS(sprintf("simulation_benchmark/trained_%s_rsf.rds", cens_id))
    } else {
      learner = readRDS(sprintf("simulation_benchmark/trained_%s.rds", cens_id))[[learner_id]]
    }

    # flexsurv learners → specifically set prediction grid to all unique event times (dense grid)
    if (inherits(learner, "LearnerSurvFlexreg")) {
      learner$param_set$set_values(times = uevents)
    }

    # get predictions
    pred = learner$predict(test_task)

    # Cox / KM / RSF → prediction grid has both unique events + censoring times
    # Interpolate to the event times
    if (!inherits(learner, "LearnerSurvFlexreg")) {
      pred$data$distr = survdistr::interp(
        x = pred$data$distr,
        eval_times = uevents,
        method = "const_surv",
        output = "surv",
        add_times = TRUE,
        check = FALSE
      )
    }

    # calculate MISE and MIAE on the dense events grid (before sparsifying for RCLL evaluation)
    mise = mise(s_true, pred$data$distr, uevents)
    miae = miae(s_true, pred$data$distr, uevents)

    # Sparsify prediction time grid
    pred$data$distr = survdistr::interp(
      x = pred$data$distr, # survival prediction matrix
      times = uevents, # original dense grid (unique event times)
      eval_times = times, # sampled events (sparse grid)
      method = "const_surv",
      output = "surv",
      add_times = TRUE,
      check = FALSE
    )

    out[[learner_id]] = suppressWarnings(
      data.table(
        learner_id = learner_id,
        mise = mise,
        miae = miae,
        rcll = pred$score(rcll)
      )
    )
  }

  combined = rbind(true_model_res, rbindlist(out), fill = TRUE)
  combined[, c("task_id", "rsmp_id", "prop", "n_times", "n_test") :=
           list(cens_id, rsmp_id, prop, length(times), n_test)]
  combined
}

execute_sim = function(bm_grid) {
  row_seq = bm_grid$config_id
  p = progressr::progressor(along = row_seq)
  results = future.apply::future_lapply(
    row_seq,
    function(i) eval_config(bm_grid[i, ], p),
    future.seed = TRUE
  )
  rbindlist(results)
}

with_progress({
  results_dt = execute_sim(sim_grid)
})

saveRDS(results_dt, "results/rcll_sim.rds")
