#' Convert datasets to tasks, record if PH assumption is satisfied or not and
#' other censoring-related info
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(mlr3misc))
library(mlr3proba)

paths = list.files(path = "data", full.names = TRUE)

task_tbl = lapply(paths, function(path) {
  # read data, convert to survival task
  data = readRDS(file = path)
  task_id = gsub(".*/([^/]+)\\.[^.]+$", "\\1", path)
  task = mlr3proba::as_task_surv(x = data, time = "time",
    event = "status", type = "right", id = task_id)

  # general info about the task
  n_obs = task$nrow
  n_vars = length(task$feature_names)
  n_factors = sum(task$feature_types$type == "factor")
  n_numeric = n_vars - n_factors

  # censoring info
  truth = task$truth()
  times  = truth[, 1]
  status = truth[, 2]
  cens_prop = sum(status == 0)/length(status)
  cens_times = times[status == 0]
  cens_tmax = max(cens_times) # max censoring time
  # percentage of censored observations that are censored administratively
  # i.e. at the last censoring time
  admin_cens_prop = sum(cens_times == cens_tmax)/length(cens_times)

  # Dependent vs random censoring
  # Check proportion of significant coefficients to predict censoring status
  coef_tbl = glm(status ~ .,
    data = task$data(cols = c(task$feature_names, "status")),
    family = "binomial") |>
    broom::tidy()
  dep_cens_prop = sum(p.adjust(coef_tbl$p.value, method = "holm") < 0.05)/nrow(coef_tbl)

  # PH vs non-PH
  cox = lrn("surv.coxph")$train(task)
  zph_test = try(survival::cox.zph(fit = cox$model), silent = TRUE)
  if (inherits(zph_test, "try-error")) {
    prop_haz = NA
  } else {
    # < 0.05 indicates PH violation
    p_value = zph_test$table["GLOBAL", "p"]
    prop_haz = ifelse(p_value < 0.05, FALSE, TRUE)
  }

  tibble(
    task = list(task),
    task_id = task_id,
    n_obs = n_obs,
    n_vars = n_vars,
    n_factors = n_factors,
    n_numeric = n_numeric,
    cens_prop = cens_prop,
    admin_cens_prop = admin_cens_prop,
    dep_cens_prop = dep_cens_prop,
    prop_haz = prop_haz
  )
}) |> bind_rows()

saveRDS(task_tbl, file = "task_tbl.rds")
