#' Investigation: Fixed test set size, varying prediction grid coarseness.
#' Grid defined by proportions (0.02, 0.05, 0.1, 0.2, 0.5, 0.8) of unique event times.
#' How much does the density approximation affects the sensitivity of RCLL to model
#' misspecification.
#'
#' Run: `Rscript simulation_benchmark/RCLL_investigation.R`
suppressWarnings({
  library(mlr3) # v1.6.0
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

# Fxied test set size (sampled from DGP for prediction)
n_test = 250

# Proportions of unique event times to use as prediction grid
proportions = c(0.02, 0.05, 0.1, 0.2, 0.5, 0.8)

# Precompute deterministic time grids for each task and proportion
# (ensures same grid across replicates)
time_grids = list()
for (task_name in names(tasks)) {
  train_task = tasks[[task_name]]
  all_times = sort(unique(c(0, train_task$unique_event_times())))  # include 0
  n_total = length(all_times)
  for (prop in proportions) {
    n_times = ceiling(prop * n_total)
    # Deterministic sampling: seed based on task and proportion
    set.seed(1000 * (which(names(tasks) == task_name)) + 100 * (which(proportions == prop)))
    sampled_times = sort(sample(all_times, size = n_times))
    time_grids[[task_name]][[as.character(prop)]] = sampled_times
  }
}

# Simulation grid: task, replicate, proportion
sim_grid = data.table::CJ(task_id = names(tasks), rsmp_id = rsmp_id, prop = proportions)
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

  p(sprintf("Config-id: %i, Task: %s, RSMP-id: %i, Prop: %f",
            config_id, cens_id, rsmp_id, prop))

  train_task = tasks[[cens_id]]

  # Retrieve precomputed time grid for this task and proportion
  times = time_grids[[cens_id]][[as.character(prop)]]

  # Seed depends on task, replicate, and proportion (n_test fixed)
  seed = 100000 * rsmp_id + 1000 * as.integer(prop * 100) + ifelse(cens_id == "low", 1, 2)
  set.seed(seed)

  # Simulate test data (using the same dense grid for true survival)
  while (TRUE) {
    test_sim = simulate(
      n = n_test,
      b0 = 1.15, b1 = 0.15, b2 = -0.55, b12 = -0.75,
      sigma = c(0.5, 1.5),
      lambdaC = if (cens_id == "low") 0.075 else 0.45,
      times = times
    )
    test_task = test_sim$task
    cp = test_task$cens_prop()
    # guard against low sample size (n) and extreme censoring in test set
    # which can cause C-index to be undefined
    if (!is.null(cp) && length(cp) == 1 && !is.na(cp) && cp > 0 && cp < 0.9) break
  }

  true_pred = test_sim$pred
  s_true = true_pred$data$distr

  # True model results
  true_model_res = suppressWarnings(
    data.table(
      learner_id = "true",
      mise = mise(s_true, s_true, times),
      miae = miae(s_true, s_true, times),
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

    # flexsurv learners → set prediction grid
    if (inherits(learner, "LearnerSurvFlexreg")) {
      learner$param_set$set_values(times = times)
    }

    pred = learner$predict(test_task)

    # Cox / KM / RSF → inteprolate to the sparse grid
    if (!inherits(learner, "LearnerSurvFlexreg")) {
      pred$data$distr = survdistr::interp(
        x = pred$data$distr,
        eval_times = times,
        method = "const_surv",
        output = "surv",
        add_times = TRUE,
        check = FALSE
      )
    }

    out[[learner_id]] = suppressWarnings(
      data.table(
        learner_id = learner_id,
        mise = mise(s_true, pred$data$distr, times),
        miae = miae(s_true, pred$data$distr, times),
        rcll = pred$score(rcll)
      )
    )
  }

  combined = rbind(true_model_res, rbindlist(out), fill = TRUE)
  combined[, c("task_id", "rsmp_id", "prop", "n_times") := list(cens_id, rsmp_id, prop, length(times))]
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
