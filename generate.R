#' Generate the simulated data
suppressPackageStartupMessages(library(tidyverse))
library(survival)
suppressPackageStartupMessages(library(coxed))
library(mlr3proba)
library(mlr3extralearners)
library(progressr)
suppressPackageStartupMessages(library(pracma))
suppressPackageStartupMessages(library(future.apply))

#' Parameters of interest:
#' 1) independent (random) vs dependent censoring (`cens_dep`)
#' 2) PH vs non-PH (time-varying coefficients) data (`prop_haz`)
#' 3) %censoring => 10,20,..,80 (`cens_prop`)
#' 4) Number of observations: 100 - 1000 (`n_obs`)

#' Fixed parameters (for each combo of 1-4 above):
#' 4) Time horizon (max event time): **365 days** - 1 year (`time_horizon`)
#' 5) Number of datasets to generate: 1000 `n_dfs`
#' 6) Number of covariates (random select): 3-10 `n_vars` (low-dim setting)

n_dfs = 100 # Number of datasets to generate PER combo of the (1)-(4) parameters
sim_grid = expand.grid(
  cens_dep = c(TRUE, FALSE),
  prop_haz = c(TRUE, FALSE),
  cens_prop = seq(from = 0.1, to = 0.8, by = 0.1),
  n_obs = seq(from = 100, to = 1000, by = 100),
  time_horizon = 365,
  n_dfs = n_dfs
) |> as_tibble()

# How many different categories of datasets to generate (per category)?
n_rows = nrow(sim_grid)
print(paste0("Number of dataset categories: ", n_rows))

# censored obs at last time point?
cens_last = FALSE

# Integrated Survival Brier Score (improper) and re-weighted version (proper)
graf_improper = msr("surv.graf", proper = FALSE, id = "graf.improper")
graf_proper   = msr("surv.graf", proper = TRUE,  id = "graf.proper")

# Progress bars
options(progressr.enable = TRUE)
handlers(global = TRUE)
handlers("progress")

# use all available CPUs
future::plan("multicore")

with_progress({
  row_seq = seq_len(n_rows)
  p = progressor(along = row_seq)

  res = future.apply::future_lapply(row_seq, function(i) {
    p(sprintf("dataset-id = %g", i))

    set.seed(i)
    sim_params = sim_grid[i,]
    n_obs = sim_params$n_obs
    n_dfs = sim_params$n_dfs
    time_horizon = sim_params$time_horizon
    prop_haz = sim_params$prop_haz
    cens_prop = sim_params$cens_prop
    cens_dep = sim_params$cens_dep

    # generate data
    index = 1
    data_list = list()
    while (TRUE) {
      if (index == n_dfs + 1) break

      # how many X variables in the data
      n_vars = sample(x = 3:10, size = 1)

      # generate dataset
      simdata = coxed::sim.survdata(
        N = n_obs,
        T = time_horizon,
        xvars = n_vars,
        num.data.frames = 1,
        type = ifelse(prop_haz, "none", "tvbeta"),
        censor = cens_prop,
        censor.cond = cens_dep
      )

      # censor all observations at t_max
      if (cens_last) {
        times = simdata$data$y
        status = simdata$data$failed

        t_max = max(times)
        indx = which(times == t_max)
        simdata$data$failed[indx] = 0
      }

      # convert to survival mlr3 task
      task = mlr3proba::as_task_surv(x = simdata$data, time = "y",
        event = "failed", type = "right", id = "coxed.sim.surv")

      # check PH assumption
      cox = lrn("surv.coxph")$train(task)
      ok = (length(cox$errors) == 0) &
           (length(cox$warnings) == 0)
      # rare case, somehow cox model didn't converge, train didn't succeed, etc
      if (!ok) next

      # p_value < 0.05 indicates PH violation
      zph_test = survival::cox.zph(fit = cox$model)
      p_value = zph_test$table["GLOBAL", "p"]

      # keep data only if it agrees with `prop_haz`
      ok = (p_value < 0.05 & !prop_haz) || (p_value > 0.05 & prop_haz)

      if (ok) {
        # train/test split
        part = partition(task, ratio = 0.7) # by default stratified

        # define Kaplan-Meier, CoxPH and AFT-Weibull models
        kaplan = lrn("surv.kaplan")
        cox = lrn("surv.coxph")
        aft = lrn("surv.parametric", form = "aft", dist = "weibull", discrete = TRUE)

        # add encapsulation for capturing errors
        kaplan$encapsulate = c(train = "evaluate", predict = "evaluate")
        cox$encapsulate    = c(train = "evaluate", predict = "evaluate")
        aft$encapsulate    = c(train = "evaluate", predict = "evaluate")

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

        #' Convert the simulated distr predictions to `mlr3proba::PredictionSurv()`
        #' `ind.survive` is a survival matrix (N x T)
        surv = simdata$ind.survive
        colnames(surv) = 1:ncol(surv) # all integer time points

        #' survival probabilities must by non-strictly decreasing
        #' checked by `survivalmodels:::assert_surv_matrix()`
        p = try(mlr3proba::.surv_return(surv = surv), silent = TRUE)

        if (inherits(p, "try-error")) {
          sim_proper = NA
          sim_improper = NA
          sim_proper_scores = NA
          sim_improper_scores = NA
        } else {
          # keep only the simulated predictions of the test set
          pred_sim = mlr3proba::PredictionSurv$new(
            row_ids = part$test, truth = task$truth(rows = part$test),
            distr = p$distr[part$test, ], crank = p$crank[part$test]
          )
          sim_proper = pred_sim$score(graf_proper, task = task, train_set = part$train)
          sim_improper = pred_sim$score(graf_improper, task = task, train_set = part$train)
          sim_proper_scores = graf_proper$scores
          sim_improper_scores = graf_improper$scores
        }

        data_list[[index]] = tibble(
          # general info about the simulated data
          n_obs = n_obs,
          n_vars = n_vars,
          prop_haz = prop_haz,
          cens_prop = cens_prop,
          cens_dep = cens_dep,
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
          aft_improper_scores = list(aft_improper_scores),
          # Simulated distr predictions
          sim_proper = sim_proper,
          sim_improper = sim_improper,
          sim_proper_scores = list(sim_proper_scores),
          sim_improper_scores = list(sim_improper_scores),
          # task and simulated survival data (takes too much space to keep)
          # task = list(task),
          # part = list(part),
          # sim_surv = list(simdata$ind.survive)
        )
        index = index + 1
      }
    }

    # return tibble
    data_list |> bind_rows()
  }, future.seed = TRUE) |> bind_rows()
})

stopifnot(nrow(res) == n_rows * n_dfs)

print("Save results")
saveRDS(res, file = "res.rds")
