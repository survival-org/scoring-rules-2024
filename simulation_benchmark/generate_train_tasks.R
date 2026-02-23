#' Generate training tasks and train learners for `sensitivity_benchmark.R`
library(mlr3)
library(mlr3proba)
source("simulate.R")

if (file.exists("data/train_tasks.rds")) {
  tasks = readRDS("data/tasks.rds")
  trained_models_low = readRDS("data/trained_models_low.rds")
  trained_models_high = readRDS("data/trained_models_high.rds")
} else {
  source("train_learners.R") # trains all learners on both tasks (low, high censoring) and saves the trained models in "trained_models.rds"
}
saveRDS(trained_models, "trained_models.rds")
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
  RSF = lrn("surv.ranger", id = "RSF", time.interest = 1000)
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

saveRDS(trained_models$low,  "trained_low.rds")
saveRDS(trained_models$high, "trained_high.rds")
