#' Investigation of ISBS Gaming (NOT included in paper)
#' 
#' Attempted to show ISBS can be "gamed" via degenerate step-function S(t) = I(t < t0)
#' 
#' FAILED: Outperforms KM/RSF on ISBS only on `rats` dataset (-13% vs KM, -12% vs RSF)
#' Worse everywhere else; always C-index=0.5, D-calib∼10⁹, RCLL correctly punishes
#' 
#' Conclusion: this "ISBS gaming" doesn't work broadly
#' 
#' Full results: results/degenerate_model_bm.rds (9 tasks × 50 resamplings)
#' Execute: `Rscript degenerate_model.R`
library(R6)
library(mlr3proba)
library(mlr3extralearners)
library(mlr3tuning)
library(paradox)
library(progressr)

# define degenerate model as an mlr3proba LearnerSurv
LearnerSurvDegenerate = R6Class("LearnerSurvDegenerate",
  inherit = LearnerSurv,
  public = list(
    initialize = function() {
      super$initialize(
        id = "surv.degenerate",
        predict_types = c("distr"),
        feature_types = c("logical", "integer", "numeric", "character", "factor", "ordered"),
        properties = "missings",
        packages = c("survival", "distr6"),
        label = "Degenerate Estimator",
        param_set = ps(
           quantile = p_dbl(0, 1)
        )
      )
    }
  ),

  private = list(
    .train = function(task) {
      list(time = task$truth()[,1L]) # store observed times
    },

    .predict = function(task) {
      quantile_ps = self$param_set$values$quantile
      times = sort(unique(self$model$time))
      surv = matrix(nrow = task$nrow, ncol = length(times))

      q_t = quantile(seq.int(length(times)), quantile_ps)[[1]]
      
      # same S for all test observations, sharp drop from 1 to 0 at q_t
      surv[, 1:floor(q_t)] = 1
      surv[, ceiling(q_t):ncol(surv)] = 0
      colnames(surv) = times
      mlr3proba::surv_return(times = times, surv = surv)
    }
  )
)

l = LearnerSurvDegenerate$new()
l$param_set$values$quantile = to_tune(0, 1)

# ISBS
isbs = msr("surv.graf")

# Tune the quantile parameter of the degenerate model
at = auto_tuner(
  tuner = tnr("grid_search", resolution = 50),
  learner = l,
  resampling = rsmp("holdout"),
  measure = isbs
)

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

# conduct benchmark (<=20 min on a modern laptop)
set.seed(42)
bm_grid = benchmark_grid(
  tasks = tasks,
  learners = c(lrns(c("surv.coxph", "surv.kaplan", "surv.ranger")), at),
  resamplings = rsmp("repeated_cv", repeats = 10, folds = 5)
)
bm = benchmark(bm_grid)

# calculate scores
measures = c(
  msr("surv.cindex", id = "C-index"),
  msr("surv.dcalib", id = "D-calib"),
  msr("surv.rcll", id = "RCLL"),
  msr("surv.graf", id = "ISBS")
)
res = bm$aggregate(measures)

res = data.table::as.data.table(res) |>
  dplyr::select(task_id, learner_id, RCLL, `C-index`, `D-calib`, `ISBS`)
saveRDS(res, file = "results/degenerate_model_bm.rds")

# plot results nicely
if (FALSE) {
  res = readRDS(file = "results/degenerate_model_bm.rds")

  res |>
  dplyr::mutate(
    learner_id = dplyr::recode(
      learner_id,
      "surv.kaplan" = "KM",
      "surv.ranger" = "RSF",
      "surv.coxph"  = "CoxPH",
      "surv.degenerate.tuned" = "Degenerate-Model"
    )
  )|>
  DT::datatable(rownames = FALSE, filter = "top", options = list(searching = TRUE)) |>
  DT::formatRound(columns = 3:7, digits = 3) |>
  DT::formatSignif(columns = 5, digits = 3) # scientific notation for D-calib
}

