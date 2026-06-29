source("weighted_RCLL.R") # add RCLL*
library(mlr3extralearners)
library(mlr3proba)
library(data.table)
library(mlr3misc)
library(progressr)

# get all survival tasks in mlr3proba
keys = data.table::as.data.table(mlr_tasks)[task_type == "surv"][["key"]]
tasks = lapply(keys, function(key) {
  tsk(key)
})
# remove actg dataset
index = which(mlr3misc::map(tasks, `[[`, "id") == "actg")
tasks[[index]] = NULL
# stratify by event status for resampling
for (i in seq_along(tasks)) {
  task = tasks[[i]]
  task$set_col_roles(cols = task$target_names[[2]], add_to = "stratum")
}

# logging
lgr::get_logger("mlr3")$set_threshold("warn")

# Progress bars
options(progressr.enable = TRUE)
handlers(global = TRUE)
handlers("progress")

# parallelization
future::plan("multicore", workers = 10)

# define learners
learners = list(
  lrn("surv.kaplan", id = "KM"), # baseline
  lrn("surv.coxph", id = "Cox"), # semi-parametric
  lrn("surv.ranger", id = "RSF"), # non-parametric
  lrn("surv.flexreg", id = "LogNorm", dist = "lognormal"), # parametric
  lrn("surv.flexreg", id = "Weibull", dist = "weibull") # parametric
)

# encapsulate learners to fallback to Kaplan-Meier if they error
learners = lapply(learners, function(learner) {
  learner$encapsulate("evaluate", fallback = lrn("surv.kaplan"))
  learner
})

# run benchmark
set.seed(42)
bm_grid = benchmark_grid(
  tasks = tasks,
  learners = learners,
  resamplings = rsmp("repeated_cv", repeats = 10, folds = 5)
)
bm = benchmark(bm_grid)

# calculate scores
measures = c(
  msr("surv.rcll", id = "RCLL"),
  msr("surv.wrcll", id = "RCLL*", weighted = TRUE), # weighted RCLL (RCLL* in the paper)
  msr("surv.cindex", id = "C-index"),
  msr("surv.dcalib", id = "D-calib"),
  msr("surv.graf", id = "ISBS")
)
res = bm$score(measures)

# store results
res = data.table::as.data.table(res) |>
  dplyr::select(task_id, learner_id, iteration, RCLL, `RCLL*`, `C-index`, `D-calib`, `ISBS`)
saveRDS(res, file = "results/real_data_bm.rds")
