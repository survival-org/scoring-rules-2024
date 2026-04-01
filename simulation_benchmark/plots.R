library(data.table)
library(ggplot2)

# Plot score diffs ----
# load results
res = readRDS("results/sim_bm.rds")

# Generate a faceted mutli-plot per measure, score difference vs. MIAE for
# different learners, faceted by n_test and n_times
plot_score_grid = function(
  res,
  y = "miae",
  measure_id = "rcll",
  task_id_val = "low",
  n_test_vals = c(10, 25, 50, 100, 250, 500, 1000),
  n_times_vals = c(10, 25, 50, 100, 250, 500, 1000),
  thr = 20,
  drop_learners = c("CoxPH")
) {
  # subset relevant data
  dt = res[
    task_id == task_id_val &
    n_test %in% n_test_vals &
    n_times %in% n_times_vals
  ]

  # Remove CoxPH and potentially other learners
  dt = dt[!learner_id %in% drop_learners]

  # --- SPECIAL HANDLING FOR RCLL ---
  if (measure_id == "rcll") {
    # full oracle (true_rcll) per replicate
    true_dt = dt[learner_id == "true",
      .(rsmp_id, n_test, n_times, true_val = true_rcll)
    ]

    # merge full oracle
    dt = merge(dt, true_dt, by = c("rsmp_id", "n_test", "n_times"))

    # compute diff vs FULL oracle
    dt[, diff := rcll - true_val]

    # change label for true model and keep it in the plot
    dt[learner_id == "true", learner_id := "Oracle (S true, f interp)"]
  } else if (measure_id == "cindex") {
    # convert to error = 1 - C
    true_dt = dt[learner_id == "true",
      .(rsmp_id, n_test, n_times, true_val = 1 - cindex)
    ]

    dt = merge(dt, true_dt, by = c("rsmp_id", "n_test", "n_times"))

    dt[, diff := (1 - cindex) - true_val]

    # remove true model (same logic as before)
    dt = dt[learner_id != "true"]
  } else {
    # default: compare to "true" learner for SBS/ISBS/wRCLL (S true, f via linear interpolation of S for wRCLL)
    true_dt = dt[learner_id == "true",
      .(rsmp_id, n_test, n_times, true_val = get(measure_id))
    ]

    dt = merge(dt, true_dt, by = c("rsmp_id", "n_test", "n_times"))

    # compute diff
    dt[, diff := get(measure_id) - true_val]

    dt = dt[learner_id != "true"]
  }

  # wRCLL outlier handling
  if (measure_id == "wrcll") {
    dt = dt[abs(diff) <= thr]
  }

  # factor ordering for grid layout
  dt[, n_test := factor(n_test, levels = sort(n_test_vals))]
  dt[, n_times := factor(n_times, levels = sort(n_times_vals))]

  # learner ordering
  learner_order = c(
    "Oracle (S true, f interp)", # full oracle for RCLL, half-oracle for others
    "LogNorm_int_shape_x2",
    "LogLog_int_shape_x2",
    "Weib_int_shape_x2",
    "LogNorm_noint_shape_x2",
    "LogNorm_int_noshape",
    "LogNorm_noint_noshape",
    "CoxPH_int",
    "KM",
    "RSF"
  )
  dt[, learner_id := factor(learner_id, levels = learner_order)]

  # colors
  cols = c(
    `Oracle (S true, f interp)` = "#000000",
    LogNorm_int_shape_x2   = "#dfce0cff",
    LogNorm_noint_shape_x2 = "#2171B5",
    LogNorm_int_noshape    = "#6BAED6",
    LogNorm_noint_noshape  = "#C6DBEF",
    LogLog_int_shape_x2    = "#FFA500",
    Weib_int_shape_x2      = "#756BB1",
    CoxPH_int              = "#238B45",
    KM                     = "#7F7F7F",
    RSF                    = "#D62728"
  )

  # label for measure
  measure_labels = c(
    rcll     = "RCLL",
    wrcll    = "Weighted RCLL",
    cindex   = "1 − C-index (error)",
    sbs_1    = "SBS (t = 1)",
    sbs_5    = "SBS (t = 5)",
    sbs_9    = "SBS (t = 9)",
    isbs     = "ISBS (τ = 9.99)",
    isbs_q90 = "ISBS (τ = q90)"
  )

  ggplot(dt, aes(x = .data[[y]], y = diff, color = learner_id)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(alpha = 0.5, size = 1.5) +
    scale_color_manual(values = cols) +
    facet_grid(
      rows = vars(n_times),
      cols = vars(n_test),
      labeller = labeller(
        n_test = function(x) paste0("n_test = ", x),
        n_times = function(x) paste0("n_times = ", x)
      )
    ) +
    labs(
      title = sprintf(
        "%s censoring DGP — %s",
        paste0(toupper(substring(task_id_val, 1, 1)), substring(task_id_val, 2)),
        measure_labels[[measure_id]]
      ),
      x = ifelse(y == "mise",
        "Mean Integrated Squared Error (distance to true S)",
        "Mean Integrated Absolute Error (distance to true S)"),
      y = "Score difference (pred − true)",
      color = "Learner"
    ) +
    theme_bw(base_size = 16, base_family = "Arial") +
    guides(color = guide_legend(override.aes = list(size = 4)))
}

## Low-cens DGP ----
# learner ids to remove from some figures (keeping only the closer to the oracle)
learner_ids = c("CoxPH", "CoxPH_int", "LogNorm_noint_noshape",
                "LogNorm_noint_shape_x2", "LogNorm_int_noshape", "RSF", "KM")
# RCLL proper
plot_score_grid(res, measure_id = "rcll", task_id_val = "low",
                n_test_vals = c(10, 50, 100, 500), n_times_vals = c(10, 100, 500))

plot_score_grid(res, measure_id = "rcll", task_id_val = "low",
                n_test_vals = c(10, 50, 100, 500), n_times_vals = c(10, 100, 500),
                drop_learners = learner_ids)

# RCLL* is not proper!!! (negative Pred - True)
plot_score_grid(res, measure_id = "wrcll", task_id_val = "low", n_test_vals = c(10, 50, 100),
                n_times_vals = c(10, 100, 500))
# C-index is not proper! (negative Pred - True)
plot_score_grid(res, measure_id = "cindex", task_id_val = "low", n_test_vals = c(10, 50, 100),
                n_times_vals = c(10, 100, 500))
# SBS
plot_score_grid(res, measure_id = "sbs_5", task_id_val = "low", n_test_vals = c(10, 50, 100),
                n_times_vals = c(10, 100, 500))
plot_score_grid(res, measure_id = "sbs_1", task_id_val = "low", n_test_vals = c(25, 50, 100),
                n_times_vals = c(25, 100, 500))
plot_score_grid(res, measure_id = "sbs_9", task_id_val = "low", n_test_vals = c(10, 50, 100),
                n_times_vals = c(10, 100, 500))
# ISBS (τ = 9.99)
plot_score_grid(res, measure_id = "isbs", task_id_val = "low", n_test_vals = c(10, 50, 100, 500),
                n_times_vals = c(10, 100, 500))
# ISBS (τ = q90)
plot_score_grid(res, measure_id = "isbs_q90", task_id_val = "low", n_test_vals = c(10, 50, 100, 500),
                n_times_vals = c(10, 100, 500))

## High-cens DGP ----
# RCLL* is not proper!!! (negative Pred - True)
plot_score_grid(res, measure_id = "wrcll", task_id_val = "high", n_test_vals = c(10, 50, 100),
                n_times_vals = c(10, 100, 500))
# C-index is not proper! (negative Pred - True)
plot_score_grid(res, measure_id = "cindex", task_id_val = "high", n_test_vals = c(10, 50, 100),
                n_times_vals = c(10, 100, 500))
# SBS
plot_score_grid(res, measure_id = "sbs_5", task_id_val = "high", n_test_vals = c(10, 50, 100),
                n_times_vals = c(10, 100, 500))
plot_score_grid(res, measure_id = "sbs_1", task_id_val = "high", n_test_vals = c(25, 50, 100),
                n_times_vals = c(25, 100, 500))
plot_score_grid(res, measure_id = "sbs_9", task_id_val = "high", n_test_vals = c(10, 50, 100),
                n_times_vals = c(10, 100, 500))
# ISBS (τ = 9.99)
plot_score_grid(res, measure_id = "isbs", task_id_val = "high", n_test_vals = c(10, 50, 100),
                n_times_vals = c(10, 100, 500))
# ISBS (τ = q90)
plot_score_grid(res, measure_id = "isbs_q90", task_id_val = "high", n_test_vals = c(10, 50, 100),
                n_times_vals = c(10, 100, 500))

# Plot all measures ----
# Plot score differences vs. MIAE/MISE for different measures, compared to the true model
# not very accurate for RCLL, C-index
plot_score_diff = function(
  res,
  y = "miae",
  task_id_val = "low",
  n_test_val = 1000,
  n_times_val = 10,
  thr = 10,
  drop_learners = c("CoxPH")
) {
  # Subset data
  dt = res[
    task_id == task_id_val &
    n_test == n_test_val &
    n_times == n_times_val
  ]

  # Remove CoxPH and potentially other learners
  dt = dt[!learner_id %in% drop_learners]

  measures = c("rcll", "wrcll", "cindex", "sbs_1", "sbs_5", "sbs_9", "isbs", "isbs_q90")

  # --- TRUE REFERENCES ---
  true_dt = dt[learner_id == "true",
    .(
      rsmp_id,
      rcll_true_full = true_rcll,
      rcll_true_interp = rcll,
      cindex_true = cindex,
      wrcll_true = wrcll,
      sbs_1_true = sbs_1,
      sbs_5_true = sbs_5,
      sbs_9_true = sbs_9,
      isbs_true = isbs,
      isbs_q90_true = isbs_q90
    )
  ]

  dt = merge(dt, true_dt, by = "rsmp_id")

  # --- COMPUTE DIFFERENCES ---

  # RCLL: compare to FULL oracle
  dt[, diff_rcll := rcll - rcll_true_full]

  # rest of measures: compare to half-oracle (S true, f via linear interpolation of S)
  dt[, diff_wrcll := wrcll - wrcll_true]
  dt[, diff_sbs_1 := sbs_1 - sbs_1_true]
  dt[, diff_sbs_5 := sbs_5 - sbs_5_true]
  dt[, diff_sbs_9 := sbs_9 - sbs_9_true]
  dt[, diff_isbs := isbs - isbs_true]
  dt[, diff_isbs_q90 := isbs_q90 - isbs_q90_true]

  # C-index → error scale
  dt[, diff_cindex := (1 - cindex) - (1 - cindex_true)]

  # --- RELABEL TRUE MODEL FOR RCLL ---
  dt[learner_id == "true", learner_id := "Oracle (S true, f interp)"]

  # --- LONG FORMAT ---
  long_dt = melt(
    dt,
    id.vars = c("learner_id", "mise", "miae", "rsmp_id"),
    measure.vars = paste0("diff_", measures),
    variable.name = "measure",
    value.name = "diff"
  )

  long_dt[, measure := sub("diff_", "", measure)]

  # Remove interpolated oracle for all EXCEPT RCLL
  long_dt = long_dt[
    !(learner_id == "Oracle (S true, f interp)" & measure != "rcll")
  ]

  # Remove wRCLL outliers
  plot_dt = long_dt[!(measure == "wrcll" & abs(diff) > thr)]

  # --- COLORS ---
  cols = c(
    `Oracle (S true, f interp)` = "#000000",
    LogNorm_int_shape_x2   = "#dfce0cff",
    LogNorm_noint_shape_x2 = "#2171B5",
    LogNorm_int_noshape    = "#6BAED6",
    LogNorm_noint_noshape  = "#C6DBEF",
    LogLog_int_shape_x2    = "#FFA500",
    Weib_int_shape_x2      = "#756BB1",
    CoxPH_int              = "#238B45",
    KM                     = "#7F7F7F",
    RSF                    = "#D62728"
  )

  # --- LABELS ---
  measure_labels = c(
    rcll     = "RCLL (vs True Likelihood)",
    wrcll    = "Weighted RCLL",
    cindex   = "1 − C-index (error)",
    sbs_1    = "SBS (t = 1)",
    sbs_5    = "SBS (t = 5)",
    sbs_9    = "SBS (t = 9)",
    isbs     = "ISBS (τ = 9.99)",
    isbs_q90 = "ISBS (τ = q90)"
  )

  plot_dt[, measure_label := factor(
    measure,
    levels = names(measure_labels),
    labels = measure_labels
  )]

  # --- ORDER ---
  learner_order = c(
    "Oracle (S true, f interp)",
    "LogNorm_int_shape_x2",
    "LogLog_int_shape_x2",
    "Weib_int_shape_x2",
    "LogNorm_noint_shape_x2",
    "LogNorm_int_noshape",
    "LogNorm_noint_noshape",
    "CoxPH_int",
    "KM",
    "RSF"
  )

  plot_dt[, learner_id := factor(learner_id, levels = learner_order)]

  # --- PLOT ---
  ggplot(plot_dt, aes(x = .data[[y]], y = diff, color = learner_id)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = cols) +
    facet_wrap(~ measure_label, scales = "free_y") +
    labs(
      title = sprintf(
        "%s censoring DGP, Test set size = %.0f, Prediction time grid size = %.0f",
        paste0(toupper(substring(task_id_val, 1, 1)), substring(task_id_val, 2)),
        n_test_val, n_times_val
      ),
      x = ifelse(
        y == "mise",
        "Mean Integrated Squared Error (distance to true S)",
        "Mean Integrated Absolute Error (distance to true S)"
      ),
      y = "Score difference (pred − true)",
      color = "Learner"
    ) +
    theme_bw(base_size = 16, base_family = "Arial") +
    guides(color = guide_legend(override.aes = list(size = 4)))
}

plot_score_diff(res, task_id_val = "low", n_test_val = 1000, n_times_val = 500)
