source("helper.R") # interpolation functions
library(mlr3proba)
library(mlr3extralearners)
library(paradox)
library(progressr)

# Define mlr3(proba) measure
MeasureSurvReweightedRCLL = R6::R6Class("MeasureSurvReweightedRCLL",
  inherit = MeasureSurv,
  public = list(
    initialize = function() {
      param_set = paradox::ps(
        proper = p_lgl(default = TRUE),
        eps = p_dbl(0, 1, default = 1e-6)
      )

      param_set$set_values(proper = TRUE, eps = 1e-6)

      super$initialize(
        id = "surv.rrcll",
        param_set = param_set,
        minimize = TRUE,
        predict_type = "distr",
        packages = "distr6",
        label = "Reweighted Right-Censored Log Loss",
        man = "mlr3proba::mlr_measures_surv.rcll",
        range = c(0, Inf)
      )

      invisible(self)
    }
  ),

  private = list(
    .score = function(prediction, task, train_set, ...) {
      pv = self$param_set$values

      if (pv$proper && is.null(task)) {
        # we need all the data to estimate S_C(t) and f_C(t) via KM
        stop("'task' is required for proper score for estimating the censoring distribution")
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

      # if proper, apply IPCW
      if (pv$proper) {
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
mlr_measures$add("surv.rrcll", MeasureSurvReweightedRCLL)

# (normal) RCLL
rcll = msr("surv.rrcll", id = "RCLL", proper = FALSE)
# (re-weighted) rRCLL/RCLL**
rrcll = msr("surv.rrcll", id = "RRCLL", proper = TRUE)

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
future::plan("multicore", workers = 15)

# conduct benchmark
set.seed(42)
bm_grid = benchmark_grid(
  tasks = tasks,
  learners = lrns(c("surv.kaplan", "surv.ranger", "surv.coxph")),
  resamplings = rsmp("repeated_cv", folds = 5)
)
bm = benchmark(bm_grid)

# calculate scores
res = bm$aggregate(c(rcll, rrcll, msr("surv.cindex"), msr("surv.dcalib"), msr("surv.graf")))

# store results
res =
  data.table::as.data.table(res) |>
  dplyr::select(task_id, learner_id, RCLL, RRCLL, surv.cindex, surv.dcalib, surv.graf)
saveRDS(res, file = "results/rrcll_bm.rds")
