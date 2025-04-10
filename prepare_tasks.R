#' Convert datasets to tasks, record if PH assumption is satisfied or not and
#' other censoring-related info
#'
#' Run: `Rscript prepare_tasks.R`
suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(mlr3misc)
  library(mlr3proba)
})

# logging
lgr::get_logger("mlr3")$set_threshold("warn")

paths = list.files(path = "data", full.names = TRUE)

task_tbl = lapply(paths, function(path) {
  # read data, convert to survival task
  data = readRDS(file = path)
  task_id = gsub(".*/([^/]+)\\.[^.]+$", "\\1", path)

  task = mlr3proba::as_task_surv(
    x = data,
    time = "time",
    event = "status",
    id = task_id
  )

  # general info about the task
  n_obs = task$nrow
  n_times = length(task$unique_times())
  n_vars = length(task$feature_names)
  n_factors = sum(task$feature_types$type == "factor")
  n_numeric = n_vars - n_factors

  # Proportion of censored observations in the dataset
  cens_prop = task$cens_prop()
  # Proportion of censored observations that are censored administratively
  # Default: censored at or after the last outcome time (0.99 quantile)
  admin_cens_prop = task$admin_cens_prop()
  # Proportion of significant coefficients that predict censoring status (proxy
  # to dependent censoring)
  dep_cens_prop = task$dep_cens_prop()
  # p_value < 0.05 indicates PH violation
  p_value = task$prop_haz()
  prop_haz = ifelse(p_value < 0.05, FALSE, TRUE)

  tibble::tibble(
    task = list(task),
    id = task_id,
    n_obs = n_obs,
    n_vars = n_vars,
    n_times = n_times,
    n_factors = n_factors,
    n_numeric = n_numeric,
    cens_prop = cens_prop,
    admin_cens_prop = admin_cens_prop,
    dep_cens_prop = dep_cens_prop,
    prop_haz = prop_haz
  )
}) |> bind_rows()

cat("Saving task table\n")
saveRDS(task_tbl, file = "task_tbl.rds")
