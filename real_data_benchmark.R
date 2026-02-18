source("weighted_RCLL.R") # add RCLL*
library(mlr3extralearners)
library(progressr)

# get all survival tasks in mlr3proba
keys = as.data.table(mlr_tasks)[task_type == "surv"][["key"]]
tasks = lapply(keys, function(key) {
  tsk(key)
})

# logging
lgr::get_logger("mlr3")$set_threshold("warn")

# Progress bars
options(progressr.enable = TRUE)
handlers(global = TRUE)
handlers("progress")

# parallelization
future::plan("multisession", workers = 15)

# conduct benchmark (<=2 min on a modern laptop)
set.seed(42)
bm_grid = benchmark_grid(
  tasks = tasks,
  learners = lrns(c("surv.kaplan", "surv.ranger", "surv.coxph")),
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
res = bm$aggregate(measures)

# store results
res = data.table::as.data.table(res) |>
  dplyr::select(task_id, learner_id, RCLL, `RCLL*`, `C-index`, `D-calib`, `ISBS`)
saveRDS(res, file = "results/real_data_bm.rds")
