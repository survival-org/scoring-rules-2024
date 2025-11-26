source("helper.R") # interpolation functions
library(mlr3proba) # 0.8.5 version
library(mlr3extralearners)
library(paradox)
library(progressr)

# Define new mlr3(proba) measure
MeasureSurvWeightedRCLL = R6::R6Class("MeasureSurvWeightedRCLL",
  inherit = MeasureSurv,
  public = list(
    initialize = function() {
      param_set = paradox::ps(
        weighted = p_lgl(default = TRUE),
        eps = p_dbl(0, 1, default = 1e-6)
      )

      param_set$set_values(weighted = TRUE, eps = 1e-6)

      super$initialize(
        id = "surv.wrcll",
        param_set = param_set,
        minimize = TRUE,
        predict_type = "distr",
        label = "Weighted Right-Censored Log Loss",
        man = "mlr3proba::mlr_measures_surv.wrcll",
        range = c(0, Inf)
      )

      invisible(self)
    }
  ),

  private = list(
    .score = function(prediction, task, train_set, ...) {
      pv = self$param_set$values

      if (pv$weighted && is.null(task)) {
        # we need all the data to estimate S_C(t) and f_C(t) via KM
        stop("'task' is required for weighted score for estimating the censoring distribution")
      }

      truth = prediction$truth
      n_obs = length(truth)
      test_times = truth[, 1L]
      test_status = truth[, 2L]
      eps = pv$eps

      # get survival matrix
      surv_mat = prediction$data$distr
      times = as.numeric(colnames(surv_mat))

      out = vapply(seq_len(n_obs), function(obs_index) {
        outcome_time = test_times[obs_index] # event time or censoring time

        # predicted survival curve for observation
        surv_pred = list(surv = surv_mat[obs_index, ], time = times)

        if (test_status[obs_index] == 1) { # event => use f(t)
          interp_pdf(surv_pred, outcome_time)
        } else { # censored => use S(t)
          interp_surv(surv_pred, outcome_time)
        }
      }, numeric(1))

      out = -log(pmax(eps, out))

      # if weighted, apply IPCW
      if (pv$weighted) {
        ghat = task$kaplan(reverse = TRUE)

        out = vapply(seq_len(n_obs), function(obs_index) {
          outcome_time = test_times[obs_index] # event time or censoring time

          # event => divide by censoring S(t)
          if (test_status[obs_index] == 1) {
            survCt = max(eps, interp_surv(ghat, outcome_time))
            out[obs_index] / survCt
          } else {
            # censored => divide by censoring density f(t)
            densityCt = max(eps, interp_pdf(ghat, outcome_time))
            out[obs_index] / densityCt
          }
        }, numeric(1))
      }

      mean(out)
    }
  )
)

# Add measure for easy use with `msr()`
mlr_measures = utils::getFromNamespace("mlr_measures", ns = "mlr3")
mlr_measures$add("surv.wrcll", MeasureSurvWeightedRCLL)

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
