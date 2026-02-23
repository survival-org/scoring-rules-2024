# Under misspecification, which scoring rule is most discriminative?
library(mlr3)
library(mlr3proba)
library(mlr3extralearners)
library(mlr3pipelines)
library(survdistr)
library(data.table)
library(future)
library(future.apply)
library(progressr)
source("simulate.R")
source("weighted_RCLL.R") # estimates both f from S and f_C from S_C (C estimated via Kaplan-Meier)

plan("multicore", workers = 14)

# Enable progress bars
options(progressr.enable = TRUE)
handlers(global = TRUE)
handlers("progress")

set.seed(42)

# Training Tasks ---- 
## 2 tasks: (low, high) censoring
n_train = 1000  # training sample size (fixed)

# ~25% total censoring (split ~50:50 between x2-groups)
sim_low = simulate(
  n = n_train,
  b0 = 1.15, b1 = 0.15, b2 = -0.55, b12 = -0.75,
  sigma = c(0.5, 1.5),
  lambdaC = 0.075
)

if (FALSE) {
  sim_low$task$kaplan(strata="x2") |> plot()
  sim_low$task$cens_prop()

  # per group
  print(mean(sim_low$data$status == 0 & sim_low$data$x2 == 0))
  print(mean(sim_low$data$status == 0 & sim_low$data$x2 == 1))
  # administrative censoring
  print(mean(sim_low$data$time == 10))
  # proportion: admin censoring / total censoring
  print(mean(sim_low$data$time == 10) / (mean(sim_low$data$status == 0)))
}

# ~65% total censoring (split ~50:50 between x2-groups), 
sim_high = simulate(
  n = n_train,
  b0 = 1.15, b1 = 0.15, b2 = -0.55, b12 = -0.75,
  sigma = c(0.5, 1.5),
  lambdaC = 0.45
)

sim_low$task$id  = "low"
sim_high$task$id = "high"

tasks = list(
  low  = sim_low$task,
  high = sim_high$task
)

# Learners ----
## 10 learners in total
learners = list(
  LogNorm_int_shape_x2 = lrn("surv.flexreg", id = "LogNorm_int_shape_x2",
    formula = survival::Surv(time, status) ~ x1 + x2 + x1:x2,
    anc = list(sdlog = ~ x2),
    dist = "lognormal"
  ),
  LogNorm_noint_shape_x2 = lrn("surv.flexreg", id = "LogNorm_noint_shape_x2",
    formula = survival::Surv(time, status) ~ x1 + x2,
    anc = list(sdlog = ~ x2),
    dist = "lognormal"
  ),
  LogNorm_int_noshape = lrn("surv.flexreg", id = "LogNorm_int_noshape",
    formula = survival::Surv(time, status) ~ x1 + x2 + x1:x2,
    dist = "lognormal"
  ),
  LogNorm_noint_noshape = lrn("surv.flexreg", id = "LogNorm_noint_noshape",
    formula = survival::Surv(time, status) ~ x1 + x2,
    dist = "lognormal"
  ),
  Weib_int_shape_x2 = lrn("surv.flexreg", id = "Weib_int_shape_x2",
    formula = survival::Surv(time, status) ~ x1 + x2 + x1:x2,
    anc = list(shape = ~ x2),
    dist = "weibull"
  ),
  LogLog_int_shape_x2 = lrn("surv.flexreg", id = "LogLog_int_shape_x2",
    formula = survival::Surv(time, status) ~ x1 + x2 + x1:x2,
    anc = list(shape = ~ x2),
    dist = "llogis"
  ),
  CoxPH = lrn("surv.coxph", id = "CoxPH"),
  CoxPH_int = po("modelmatrix", formula =  ~ -1 + x1 + x2) %>>%
              lrn("surv.coxph") |>
              as_learner(),
  KM = lrn("surv.kaplan", id = "KM"),
  RSF = lrn("surv.ranger", id = "RSF", time.interest = 1000, num.threads = 4)
)
learners$CoxPH_int$id = "CoxPH_int"

# Train once each task ----
bm_grid = benchmark_grid(
  tasks = tasks,
  learners = learners,
  resamplings = rsmp("insample")
)
bm = benchmark(bm_grid, store_models = TRUE)

bm_tbl = as.data.table(bm)
bm_tbl[, model := bm_tbl$learner]

# list of lists (1st level tasks (low, high), 2nd level: trained models)
trained_models = lapply(
  split(bm_tbl, bm_tbl$task_id),
  function(dt) {
    setNames(dt$model, dt$learner_id)
  }
)

# Benchmark Grid ----
task_id = names(trained_models) # low, high censoring tasks
n_rsmps = 100 # number of Monte Carlo repetitions (sampling `n_test` obs for prediction)
rsmp_id = seq_len(n_rsmps) # 100 test sets
n_test_sizes = c(10, 25, 50, 100, 250, 500, 1000) # number of test observations (sampled from DGP for prediction)
n_times_grid = c(10, 25, 50, 100, 250, 500, 1000) # number of prediction time points
bench_grid = data.table::CJ(task_id, rsmp_id, n_test_sizes, n_time_grid)
bench_grid[, config_id := .I]

# Measures ----
rcll = msr("surv.rcll")
cindex = msr("surv.cindex")
wrcll = msr("surv.wrcll")
sbs_1 = msr("surv.graf", integrated = FALSE, times = 1) # early
sbs_5 = msr("surv.graf", integrated = FALSE, times = 5) # ~median 
sbs_9 = msr("surv.graf", integrated = FALSE, times = 9) # late
# 2 ISBS variants (per task): change the upper limit of integration
# 1) up to admin censoring time (tau = 9.99)
# 2) up to 90% quantile of observed event times
q90_low = quantile(tasks$low$unique_event_times(), 0.9) # ~ 5.3
q90_high = quantile(tasks$high$unique_event_times(), 0.9) # ~ 3.4
tau = 9.99 # close to admin censoring time

isbs = msr("surv.graf", integrated = TRUE, times = seq(0.05, tau, length.out = 100))
isbs_q90_low = msr("surv.graf", integrated = TRUE, times = seq(0.05, q90_low, length.out = 100))
isbs_q90_high = msr("surv.graf", integrated = TRUE, times = seq(0.05, q90_high, length.out = 100))

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
eval_config = function(row) {
  source("weighted_RCLL.R") # to ensure the function is available in each parallel worker
  cens    = row$task_id
  rsmp_id = row$rsmp_id
  n       = row$n_test_sizes
  n_times = row$n_time_grid
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
  for (learner_id in names(trained_models[[cens]])) {
    learner = trained_models[[cens]][[learner_id]]$clone(deep = TRUE)

    # flexsurv learners → set prediction grid
    if (inherits(learner, "LearnerSurvFlexreg")) {
      learner$param_set$set_values(times = times)
    }

    pred = learner$predict(test_task)

    # Cox / KM / RSF → hack: discretize prediction time grid 
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
        rcll = pred$score(rcll), # true S, estimated f!
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

options(future.globals.maxSize = 1500 * 1024^2) # 1.5 GB (to not hit some limits)
bench_list = split(bench_grid, bench_grid$config_id)
results = future.apply::future_lapply(
  bench_list[1:5],
  eval_config,
  future.seed = TRUE
)

results_dt = rbindlist(results)
saveRDS(results_dt, "results/sens_bench_results.rds")

## Setting 1 - Notes
# RCLL: 
# - Distributional assumption not as important as interaction and correct scaling.
# - Lower Cens has worse RCLL than higher cens. Structurally the same observation though
# - with increasing misspecification in the log normal, decreasing RCLL
