#' Compare ISBS proper vs improper on some real-world datasets
suppressPackageStartupMessages(library(tidyverse))
library(mlr3proba)
library(mlr3extralearners)
library(progressr)
library(future.apply)

#' see `prepare_tasks.R`
task_tbl = readRDS(file = "task_tbl.rds")

# Progress bars
options(progressr.enable = TRUE)
handlers(global = TRUE)
handlers("progress")

# use all available CPUs
future::plan("multicore")

# how many times to split to train and test set each dataset?
n_rsmps = 100

# Integrated Survival Brier Score (improper) and re-weighted version (proper)
graf_improper = msr("surv.graf", proper = FALSE, id = "graf.improper")
graf_proper   = msr("surv.graf", proper = TRUE,  id = "graf.proper")

with_progress({
  row_seq = seq_len(nrow(task_tbl))
  p = progressor(along = row_seq)

  bench_res = future.apply::future_lapply(row_seq, function(i) {
    set.seed(i)
    task = task_tbl[i, ]$task[[1]]
    p(sprintf("task = %s", task$id))

    # define Kaplan-Meier, CoxPH and AFT-Weibull models
    kaplan = lrn("surv.kaplan")
    cox = lrn("surv.coxph")
    aft = lrn("surv.parametric", form = "aft", dist = "weibull", discrete = TRUE)

    # add encapsulation for capturing errors
    kaplan$encapsulate = c(train = "evaluate", predict = "evaluate")
    cox$encapsulate    = c(train = "evaluate", predict = "evaluate")
    aft$encapsulate    = c(train = "evaluate", predict = "evaluate")

    # split each dataset a number of times to train/test sets (repeated holdout)
    lapply(1:n_rsmps, function(rsmp_num) {
      part = partition(task, ratio = 0.8) # by default stratified

      # train
      kaplan$train(task, row_ids = part$train)
      cox$train(task, row_ids = part$train)
      aft$train(task, row_ids = part$train)

      # evaluate graf proper and improper on the test set
      # using various models, but check if training succeeded first

      # Kaplan-Meier
      if (length(kaplan$errors) == 0) {
        # predict
        pred_kaplan = kaplan$predict(task, row_ids = part$test)

        # calculate graf scores
        km_proper = pred_kaplan$score(graf_proper, task = task, train_set = part$train)
        km_improper = pred_kaplan$score(graf_improper, task = task, train_set = part$train)
        km_proper_scores = graf_proper$scores
        km_improper_scores = graf_improper$scores
      } else {
        km_proper = NA
        km_improper = NA
        km_proper_scores = NA
        km_improper_scores = NA
      }

      # Cox
      if (length(cox$errors) == 0) {
        # predict
        pred_cox = cox$predict(task, row_ids = part$test)

        # calculate graf scores
        cox_proper = pred_cox$score(graf_proper, task = task, train_set = part$train)
        cox_improper = pred_cox$score(graf_improper, task = task, train_set = part$train)
        cox_proper_scores = graf_proper$scores
        cox_improper_scores = graf_improper$scores
      } else {
        cox_proper = NA
        cox_improper = NA
        cox_proper_scores = NA
        cox_improper_scores = NA
      }

      # AFT
      if (length(aft$errors) == 0) {
        # predict
        pred_aft = aft$predict(task, row_ids = part$test)

        # calculate graf scores
        aft_proper = pred_aft$score(graf_proper, task = task, train_set = part$train)
        aft_improper = pred_aft$score(graf_improper, task = task, train_set = part$train)
        aft_proper_scores = graf_proper$scores
        aft_improper_scores = graf_improper$scores
      } else {
        aft_proper = NA
        aft_improper = NA
        aft_proper_scores = NA
        aft_improper_scores = NA
      }

      tibble(
        task_id = task$id,
        aft = list(aft),
        # KM
        km_proper = km_proper,
        km_improper = km_improper,
        km_proper_scores = list(km_proper_scores),
        km_improper_scores = list(km_improper_scores),
        # Cox
        cox_proper = cox_proper,
        cox_improper = cox_improper,
        cox_proper_scores = list(cox_proper_scores),
        cox_improper_scores = list(cox_improper_scores),
        # AFT
        aft_proper = aft_proper,
        aft_improper = aft_improper,
        aft_proper_scores = list(aft_proper_scores),
        aft_improper_scores = list(aft_improper_scores)
      )
    }) |> bind_rows()
  }, future.seed = TRUE) |> bind_rows()
})

saveRDS(bench_res, file = "bench_res.rds")
