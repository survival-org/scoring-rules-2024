#' Under misspecification, which scoring rule is most discriminative?
#' `Experiment`: simulate data from a known DGP, evaluate different scoring rules on
#' predictions from the true model and various misspecified models (Cox, RSF, KM,
#' flexible parametric) across large sample sizes (n_test = 1000) and dense
#' prediction time grid (unique event times from the training task).
#' All learners use the same grid.
#'
#' Run: `Rscript simulation_benchmark/execute_dense.R`
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
plan("multisession", workers = 50)
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

# Fixed test set size
n_test = 1000

# Grid: only task and replicate
sim_grid = data.table::CJ(task_id = names(tasks), rsmp_id = rsmp_id)
sim_grid[, config_id := .I]

# Measures
# RCLL => right-censored negative log-likelihood based on predicted S(t) and
# estimated f(t) via linear interpolation of S(t)
rcll = msr("surv.rcll")
cindex = msr("surv.cindex")
sbs = msr("surv.graf", integrated = FALSE, times = 5) # chosen time is a bit after the crossing
# ISBS: choose upper limit of integration as the 90% percentile of observed times (different per task)
q90_low = unname(quantile(tasks$low$times(), 0.9)) # ~ 6.19
q90_high = unname(quantile(tasks$high$times(), 0.9)) # ~ 3.35
# we keep evaluation grid for SBS fixed, as we want to test its sensitivity to the choice of the grid
isbs_q90_low = msr("surv.graf", integrated = TRUE, times = seq(0, q90_low, length.out = 100))
isbs_q90_high = msr("surv.graf", integrated = TRUE, times = seq(0, q90_high, length.out = 100))

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

# Eval function
eval_config = function(row, p) {
  config_id = row$config_id
  cens_id = row$task_id
  rsmp_id = row$rsmp_id

  p(sprintf("Config-id: %i, Task: %s, RSMP-id: %i, n_test = %i (fixed)",
            config_id, cens_id, rsmp_id, n_test))

  train_task = tasks[[cens_id]]

  # Dense prediction grid: unique event times from training task
  times = train_task$unique_event_times()

  # Reproducible seed based on task and replicate
  seed = 100000 * rsmp_id + ifelse(cens_id == "low", 1, 2)
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

  # Evaluate TRUE model's predictions
  isbs = if (cens_id == "low") isbs_q90_low else isbs_q90_high

  true_model_res = suppressWarnings(
    data.table(
      learner_id = "true",
      mise = mise(s_true, s_true, times),
      miae = miae(s_true, s_true, times),
      cindex = true_pred$score(cindex),
      rcll = true_pred$score(rcll),
      true_rcll = test_sim$rcll,
      sbs = true_pred$score(sbs, task = train_task, train_set = train_task$row_ids),
      isbs = true_pred$score(isbs, task = train_task, train_set = train_task$row_ids)
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

    # Cox / KM / RSF → use only the unique event times (intepolate here just to
    # be on the safe side, by default these learners will have the prediction grid
    # on the total train time points, events + censoring)
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
        cindex = pred$score(cindex),
        rcll = pred$score(rcll),
        sbs = pred$score(sbs, task = train_task, train_set = train_task$row_ids),
        isbs = pred$score(isbs, task = train_task, train_set = train_task$row_ids)
      )
    )
  }

  combined = rbind(true_model_res, rbindlist(out), fill = TRUE)
  combined[, c("task_id", "rsmp_id", "n_test") := list(cens_id, rsmp_id, n_test)]
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

saveRDS(results_dt, "results/sim_bm_dense.rds")
