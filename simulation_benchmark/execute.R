#' Under misspecification, which scoring rule is most discriminative?
#' `Experiment`: simulate data from a known DGP, evaluate different scoring rules on 
#' predictions from the true model and various misspecified models (Cox, RSF, KM,
#' flexible parametric) across different test set sizes and prediction time grid 
#' sizes.
#' 
#' Run: `Rscript simulation_benchmark/execute.R`
library(mlr3)
library(mlr3proba)
library(mlr3extralearners)
library(mlr3pipelines)
library(survdistr)
library(data.table)
library(future)
library(future.apply)
library(progressr)
source("simulation_benchmark/simulate.R")
source("weighted_RCLL.R") # estimates both f from S and f_C from S_C (C estimated via Kaplan-Meier)

plan("multicore", workers = 14)

# Enable progress bars
options(progressr.enable = TRUE)
handlers(global = TRUE)
handlers("progress")

set.seed(42)

# Training Tasks ----
if (!file.exists("simulation_benchmark/tasks.rds") && 
    !file.exists("simulation_benchmark/trained_low.rds") && 
    !file.exists("simulation_benchmark/trained_high.rds")) {
  message("Execute `gen_train_tasks_and_models.R` to generate train tasks and train the models.")
}

# Benchmark Grid ----
learner_ids = readRDS("simulation_benchmark/learner_ids.rds")
tasks = readRDS("simulation_benchmark/tasks.rds")
task_id = names(tasks) # low and high censoring tasks
n_rsmps = 100 # number of Monte Carlo repetitions (sampling `n_test` obs for prediction)
rsmp_id = seq_len(n_rsmps) # 100 test sets
n_test_sizes = c(10, 25, 50, 100, 250, 500, 1000) # number of test observations (sampled from DGP for prediction)
n_times_grid = c(10, 25, 50, 100, 250, 500, 1000) # number of prediction time points
sim_grid = data.table::CJ(task_id, rsmp_id, n_test_sizes, n_times_grid)
sim_grid[, config_id := .I]

# Measures ----
rcll = msr("surv.rcll")
cindex = msr("surv.cindex")
sbs_1 = msr("surv.graf", integrated = FALSE, times = 1) # early
sbs_5 = msr("surv.graf", integrated = FALSE, times = 5) # ~median 
sbs_9 = msr("surv.graf", integrated = FALSE, times = 9) # late
# 2 ISBS variants (per task): change the upper limit of integration
# 1) up to admin censoring time (tau = 9.99)
# 2) up to 90% quantile of observed event times
q90_low = quantile(tasks$low$unique_event_times(), 0.9) # ~ 5.3
q90_high = quantile(tasks$high$unique_event_times(), 0.9) # ~ 3.4
tau = 9.99 # close to admin censoring time

# we keep evaluation grid for SBS fixed, as we want to test its sensitivity to
# the granularity of the prediction time grid
isbs = msr("surv.graf", integrated = TRUE, times = seq(0.05, tau, length.out = 300))
isbs_q90_low = msr("surv.graf", integrated = TRUE, times = seq(0.05, q90_low, length.out = 300))
isbs_q90_high = msr("surv.graf", integrated = TRUE, times = seq(0.05, q90_high, length.out = 300))

# measure distance between predictions and true S(t) (mean integrated squared error)
mise = function(S_true, S_pred, times) {
  dt_vec = diff(times)
  sq_diff = (S_true - S_pred)^2
  mean(rowSums(0.5 * (sq_diff[, -ncol(sq_diff)] + sq_diff[, -1]) * dt_vec))
}

# alternative S(t) distance: mean integrated absolute error
miae = function(S_true, S_pred, times) {
  dt_vec = diff(times)
  abs_diff = abs(S_true - S_pred)
  mean(rowSums(0.5 * (abs_diff[, -ncol(abs_diff)] + abs_diff[, -1]) * dt_vec))
}

# Eval function ----
eval_config = function(row, p) {
  # to ensure that the functions inside the sourced file are available in each parallel worker
  source("weighted_RCLL.R")
  wrcll = msr("surv.wrcll")
  cens    = row$task_id
  rsmp_id = row$rsmp_id
  n       = row$n_test_sizes
  n_times = row$n_times_grid

  # Notify progress via `p = progressr::progressor()`
  p(sprintf("Task: %s, RSMP-id: %s, N_test: %i, Time grid: %i", cens, rsmp_id, n, n_times))

  times = seq(0, 9.99, length.out = n_times)
  train_task = tasks[[cens]]

  # reproducibility
  set.seed(row$config_id)

  # simulate test data
  test_sim = simulate(
    n = n,
    b0 = 1.15, b1 = 0.15, b2 = -0.55, b12 = -0.75,
    sigma = c(0.5, 1.5),
    lambdaC = if (cens == "low") 0.075 else 0.45,
    times = times
  )
  test_task = test_sim$task
  true_pred = test_sim$pred
  s_true = true_pred$data$distr

  # Evaluate true model's predictions
  isbs_q90 = if (cens == "low") isbs_q90_low else isbs_q90_high

  true_model_res = suppressWarnings(
    data.table(
      learner_id = "true",
      mise = mise(S_true = s_true, S_pred = s_true, times = times),
      miae = miae(S_true = s_true, S_pred = s_true, times = times),
      cindex = true_pred$score(cindex),
      rcll = true_pred$score(rcll), # true S(t), estimated f(t) via linear interpolation of S(t)
      true_rcll = test_sim$rcll, # true S(t) + true f(t)
      wrcll = true_pred$score(wrcll, task = train_task),
      sbs_1 = true_pred$score(sbs_1, task = train_task, train_set = train_task$row_ids),
      sbs_5 = true_pred$score(sbs_5, task = train_task, train_set = train_task$row_ids),
      sbs_9 = true_pred$score(sbs_9, task = train_task, train_set = train_task$row_ids),
      isbs = true_pred$score(isbs, task = train_task, train_set = train_task$row_ids),
      isbs_q90 = true_pred$score(isbs_q90, task = train_task, train_set = train_task$row_ids)
    )
  )

  # Evaluate learner predictions
  out = list()
  for (learner_id in learner_ids) {
    if (learner_id == "RSF") {
      # load trained RSF model separately (as it is large in size)
      learner = readRDS(sprintf("simulation_benchmark/trained_%s_rsf.rds", cens))
    } else {
      learner = readRDS(sprintf("simulation_benchmark/trained_%s.rds", cens))[[learner_id]]
    }

    # flexsurv learners → set prediction grid
    if (inherits(learner, "LearnerSurvFlexreg")) {
      learner$param_set$set_values(times = times)
    }

    pred = learner$predict(test_task)

    # Cox / KM / RSF → hack: discretize predicted time grid 
    # via constant interpolation
    if (!inherits(learner, "LearnerSurvFlexreg")) {
      pred$data$distr = survdistr::mat_interp(
        x = pred$data$distr,
        eval_times = times,
        constant = TRUE
      )
    }

    out[[learner_id]] = suppressWarnings(
      data.table(
        learner_id = learner_id,
        mise = mise(S_true = s_true, S_pred = pred$data$distr, times = times),
        miae = miae(S_true = s_true, S_pred = pred$data$distr, times = times),
        cindex = pred$score(cindex),
        rcll = pred$score(rcll), # pred S(t), estimated f(t) via linear interpolation of S(t)
        wrcll = pred$score(wrcll, task = train_task),
        sbs_1 = pred$score(sbs_1, task = train_task, train_set = train_task$row_ids),
        sbs_5 = pred$score(sbs_5, task = train_task, train_set = train_task$row_ids),
        sbs_9 = pred$score(sbs_9, task = train_task, train_set = train_task$row_ids),
        isbs = pred$score(isbs, task = train_task, train_set = train_task$row_ids),
        isbs_q90 = pred$score(isbs_q90, task = train_task, train_set = train_task$row_ids)
      )
    )
  }

  lrn_res = rbindlist(out)
  combined_dt = rbind(true_model_res, lrn_res, fill = TRUE)
  # add the following columns and reorder
  new_cols = c("task_id", "rsmp_id", "n_test", "n_times")
  combined_dt[, (new_cols) := list(cens, rsmp_id, n, n_times)]
  setcolorder(combined_dt, c(new_cols, setdiff(names(combined_dt), new_cols)))

  combined_dt
}

execute_sim = function(bm_grid) {
  row_seq = sim_grid$config_id
  p = progressr::progressor(along = row_seq) # progress tracking

  results = future.apply::future_lapply(
    row_seq,
    function(i) eval_config(bm_grid[i, ], p),
    future.seed = TRUE
  )

  rbindlist(results)
}

# options(future.globals.maxSize = 1500 * 1024^2) # 1.5 GB (avoid hitting RAM limits)
with_progress({
  results_dt = execute_sim(sim_grid)
})

saveRDS(results_dt, "results/sim_bm.rds")
