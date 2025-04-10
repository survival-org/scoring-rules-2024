#' Compare ISBS proper vs improper on some real-world datasets
#'
#' Run: `Rscript run_bench.R`
suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(mlr3proba)
  library(mlr3extralearners)
  library(progressr)
  library(future.apply)
})

#' see `prepare_tasks.R` for converting datasets from saved files to `mlr3` tasks
task_tbl = readRDS(file = "task_tbl.rds")

# logging
lgr::get_logger("mlr3")$set_threshold("warn")

# Progress bars
options(progressr.enable = TRUE)
handlers(global = TRUE)
handlers("progress")

# parallelization
future::plan("multicore", workers = 10)

# how many times to train/test split each dataset?
n_rsmps = 5

# models
learners = c("cox", "aft", "rsf")

# Construct the grid for rsmp_id for each dataset
bm_grid = lapply(seq_len(nrow(task_tbl)), function(i) {
  task_id = task_tbl[i, ]$id
  expand.grid(
    task_id = task_id,
    lrn_id = learners,
    rsmp_id = 1:n_rsmps,
    stringsAsFactors = FALSE
  )
}) |> bind_rows()

with_progress({
  row_seq = seq_len(nrow(bm_grid))
  p = progressor(along = row_seq)

  bench_res = future.apply::future_lapply(row_seq, function(i) {
    set.seed(i)

    # Task
    task_id = bm_grid[i, "task_id"]
    task = task_tbl |> filter(id == task_id) |> pull(task)
    task = task[[1]]$clone()
    # Resampling iteration id
    rsmp_id = bm_grid[i, "rsmp_id"]
    # Learner
    lrn_id = bm_grid[i, "lrn_id"]

    # print message
    p(sprintf("Task: %s, Learner: %s, Resampling Iteration: %s", task_id, lrn_id, rsmp_id))

    # create mlr3 learner (CoxPH, AFT-Weibull or RSF)
    if (lrn_id == "cox") {
      learner = lrn("surv.coxph")
    } else if (lrn_id == "aft") {
      learner = lrn("surv.parametric", form = "aft", dist = "weibull", discrete = TRUE)
    } else if (lrn_id == "rsf") {
      learner = lrn("surv.ranger", num.trees = 500, splitrule = "logrank")
      learner$timeout = c(train = 25, predict = 5)
    }

    # add encapsulation for capturing errors
    suppressWarnings(learner$encapsulate(method = "evaluate", fallback = lrn("surv.kaplan")))

    # split dataset to train/test sets
    task$set_col_roles(cols = "status", add_to = "stratum")
    part = partition(task, ratio = 0.8) # by default stratified

    # train learner
    learner$train(task, row_ids = part$train)

    # if training had no errors
    if (length(learner$errors) == 0) {
      # Define metrics: ISBS (improper) and RISBS (proper)
      # Use 80% quantile of outcome times in the train set as time horizon
      train_times = task$times(rows = part$train)
      t_max = unname(quantile(train_times, probs = 0.8))
      graf_improper = msr("surv.graf", proper = FALSE, id = "graf.improper",
                          t_max = t_max, eps = 0.01)
      graf_proper = msr("surv.graf", proper = TRUE, id = "graf.proper",
                        t_max = t_max, eps = 0.01)

      # Preprocessing step: keep test observations with time points strictly smaller
      # than the max time point in the training data to minimize effects of extrapolation
      #test_times = task$times(rows = part$test)
      #part$test = part$test[test_times < max(train_times)]

      # keep censoring status of the test observations
      test_times = task$times(rows = part$test)
      test_status = task$status(rows = part$test)
      # test_status = test_status[test_times <= t_max] # with t <= t_max

      # predict on the test set
      pred = learner$predict(task, row_ids = part$test)

      # PROPER score + obs-wise scores
      risbs = unname(pred$score(graf_proper, task = task, train_set = part$train))
      risbs_scores = graf_proper$scores

      # IMPROPER score + obs-wise scores
      isbs = unname(pred$score(graf_improper, task = task, train_set = part$train))
      isbs_scores = graf_improper$scores

      tibble::tibble(
        task_id = task_id,
        lrn_id = lrn_id,
        rsmp_id = rsmp_id,
        test_status = list(test_status),
        risbs = risbs,
        risbs_scores = list(risbs_scores),
        isbs = isbs,
        isbs_scores = list(isbs_scores)
      )
    } else {

    }
  }, future.seed = TRUE) |> bind_rows()
})

cat("------------\nSaving results\n")
saveRDS(bench_res, file = "bench_res5.rds")
