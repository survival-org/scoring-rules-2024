source("helper.R") # interpolation functions
library(mlr3)
library(mlr3proba) # 0.8.5 version
library(paradox)

# Define new mlr3(proba) measure: RCLL*. RCLL is a subcase (weighted = FALSE).
# RCLL* = -[δ * log(f(t))/S_C(t) + (1-δ) * S(t)/f_C(t)]
# (f,S) => predicted distribution
# (f_C,S_C) => censoring distribution
# δ = 1 => event, δ = 0 => censored
MeasureSurvWeightedRCLL = R6::R6Class("MeasureSurvWeightedRCLL",
  inherit = mlr3proba::MeasureSurv,
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
        # we use the WHOLE dataset (train+test set) to estimate S_C(t) and f_C(t) via KM
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
# check measure is okay to use 
stopifnot(!is.null(msr("surv.wrcll")))
