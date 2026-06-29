# Plot score differences vs. MIAE/MISE for one measure and multiple test set sizes
# Score diffs are (predicted - true)
plot_score_diff_by_ntest = function(
  res,
  measure_id = "rcll",# "rcll", "cindex", "sbs", "isbs"
  y = "miae", # "mise" or "miae"
  task_id_val = "low", # "low" or "high"
  n_test_vals = c(10, 25, 50, 100, 250, 500), # test sizes to facet
  drop_learners = c("CoxPH"),
  rcll_truth = "interp", # "true" or "interp"
  legend_position = "right" # "right", "left", "top", "bottom", "none"
) {
  # subset to the chosen task and test sizes
  dt = res[task_id == task_id_val & n_test %in% n_test_vals]
  # Remove unwanted learners
  dt = dt[!learner_id %in% drop_learners]
  # Rename learner_ids
  label_map = c(
    "LogNorm_int_shape_x2"   = "Oracle",
    "LogLog_int_shape_x2"    = "LLogis",
    "Weib_int_shape_x2"      = "Weibull",
    "LogNorm_noint_shape_x2" = "LogNorm_scale",
    "LogNorm_int_noshape"    = "LogNorm_int",
    "LogNorm_noint_noshape"  = "LogNorm",
    "CoxPH_int"              = "Cox_int",
    "RSF"                    = "RSF",
    "KM"                     = "KM"
  )

  dt[, learner_id := ifelse(
    learner_id %in% names(label_map),
    label_map[learner_id],
    learner_id
  )]

  # True references
  true_dt = dt[learner_id == "true",
    .(rsmp_id,
      n_test,
      rcll_true_full = true_rcll,
      rcll_true_interp = rcll,
      true_cindex = cindex,
      true_sbs = sbs,
      true_isbs = isbs)
  ]

  dt = merge(dt, true_dt, by = c("rsmp_id", "n_test"))

  # Differences
  if (measure_id == "rcll") {
    if (rcll_truth == "interp") {
      dt[, diff := rcll - rcll_true_interp]
    } else {
      dt[, diff := rcll - rcll_true_full]
    }
  } else if (measure_id == "cindex") {
    dt[, diff := true_cindex - cindex]
  } else if (measure_id == "sbs") {
    dt[, diff := sbs - true_sbs]
  } else if (measure_id == "isbs") {
    dt[, diff := isbs - true_isbs]
  }

  # Remove true model
  dt = dt[learner_id != "true"]

  # Ordering
  dt[, n_test := factor(n_test, levels = sort(n_test_vals))]

  learner_order = c(
    "Oracle",
    "LLogis",
    "Weibull",
    "LogNorm_scale",
    "LogNorm_int",
    "LogNorm",
    "Cox_int",
    "RSF",
    "KM"
  )
  present_learners = unique(dt$learner_id)
  learner_order = learner_order[learner_order %in% present_learners]
  dt[, learner_id := factor(learner_id, levels = learner_order)]

  cols = c(
    Oracle        = "#dfce0cff",
    LLogis        = "#FFA500",
    Weibull       = "#756BB1",
    LogNorm_scale = "#2171B5",
    LogNorm_int   = "#6BAED6",
    LogNorm       = "#C6DBEF",
    Cox_int       = "#238B45",
    RSF           = "#D62728",
    KM            = "#7F7F7F"
  )
  cols = cols[names(cols) %in% present_learners]

  measure_labels = c(
    rcll   = if (rcll_truth == "interp") "RCLL" else "RCLL (difference to true likelihood)",
    sbs    = "SBS (t = 5)",
    isbs   = "ISBS (τ* = q90)",
    cindex = "1 − C‑index (error)"
  )

  ggplot(dt, aes(x = .data[[y]], y = diff, color = learner_id)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = cols) +
    facet_grid(
      . ~ n_test,
      labeller = labeller(
        n_test = function(x) paste0("n['test']==", x),
        .default = label_parsed
      )
    ) +
    labs(
      title = sprintf(
        "%s censoring task — %s",
        paste0(toupper(substring(task_id_val, 1, 1)), substring(task_id_val, 2)),
        measure_labels[[measure_id]]
      ),
      x = ifelse(y == "mise", "MISE (distance to true S)", "MIAE (distance to true S)"),
      y = expression(L[pred] - ~L[true]),
      color = "Model"
    ) +
    theme_bw(base_size = 14, base_family = "Arial") +
    theme(
      strip.background = element_rect(fill = "grey95"),
      strip.text = element_text(face = "bold"),
      legend.position = if (legend_position == "none") "none" else legend_position,
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12)
    ) +
    guides(color = guide_legend(override.aes = list(size = 5)))
}

res = readRDS("results/sim_bm_dense.rds")

# low censoring task
p1 = plot_score_diff_by_ntest(res, measure_id = "rcll", task_id_val = "low", legend_position = "none")
p2 = plot_score_diff_by_ntest(res, measure_id = "isbs", task_id_val = "low", legend_position = "bottom") +
  guides(color = guide_legend(override.aes = list(size = 5, alpha = 1), nrow = 1, byrow = TRUE)) +
  theme(legend.direction = "horizontal")
p1 / p2
ggsave(filename = "sensitivity_small_n_low_cens.png", dpi = 600, height = 8, width = 13, units = "in")
ggsave(filename = "sensitivity_small_n_low_cens.pdf", dpi = 600, height = 8, width = 13, units = "in", device = cairo_pdf)

# high censoring task
p3 = plot_score_diff_by_ntest(res, measure_id = "rcll", task_id_val = "high", legend_position = "none")
p4 = plot_score_diff_by_ntest(res, measure_id = "isbs", task_id_val = "high", legend_position = "bottom") +
  guides(color = guide_legend(override.aes = list(size = 5, alpha = 1), nrow = 1, byrow = TRUE)) +
  theme(legend.direction = "horizontal")
p3 / p4
ggsave(filename = "sensitivity_small_n_high_cens.png", dpi = 600, height = 8, width = 13, units = "in")
ggsave(filename = "sensitivity_small_n_high_cens.pdf", dpi = 600, height = 8, width = 13, units = "in", device = cairo_pdf)

plot_rcll_diff_by_grid = function(
  res,
  y = "miae", # "mise" or "miae"
  task_id_val = "low", # "low" or "high"
  n_test_val = 500, # test set size to display
  prop_vals = c(0.02, 0.05, 0.1, 0.2, 0.5, 0.8, 1),
  drop_learners = c("CoxPH"),
  legend_position = "right" # "right", "left", "top", "bottom", "none"
) {
  # subset to the chosen task and test size
  dt = res[task_id == task_id_val & prop %in% prop_vals & n_test == n_test_val]
  # Remove unwanted learners
  dt = dt[!learner_id %in% drop_learners]
  # Rename learner_ids
  label_map = c(
    "LogNorm_int_shape_x2"   = "Oracle",
    "LogLog_int_shape_x2"    = "LLogis",
    "Weib_int_shape_x2"      = "Weibull",
    "LogNorm_noint_shape_x2" = "LogNorm_scale",
    "LogNorm_int_noshape"    = "LogNorm_int",
    "LogNorm_noint_noshape"  = "LogNorm",
    "CoxPH_int"              = "Cox_int",
    "RSF"                    = "RSF",
    "KM"                     = "KM"
  )

  dt[, learner_id := ifelse(
    learner_id %in% names(label_map),
    label_map[learner_id],
    learner_id
  )]

  # Extract true references (per replicate + prop)
  true_dt = dt[learner_id == "true",
    .(rsmp_id,
      prop,
      rcll_true_full = true_rcll, # true S + true f
      rcll_true_interp = rcll # true S + interpolated f
    )
  ]

  dt = merge(dt, true_dt, by = c("rsmp_id", "prop"))

  # Compute differences relative to TRUE likelihood
  dt[, diff := rcll - rcll_true_full]

  # Rename true model with interpolated density
  dt[learner_id == "true", learner_id := "Truth (true S, interp f)"]

  # Factor ordering
  # Create facet label: B = n_times (percentage)
  dt[, facet_label := sprintf("B == %d~'('*%d*'%%'*')'", n_times, round(prop * 100))]
  dt[, facet_label := factor(facet_label, levels = dt[order(prop), unique(facet_label)])]

  learner_order = c(
    "Truth (true S, interp f)",
    "Oracle",
    "LLogis",
    "Weibull",
    "LogNorm_scale",
    "LogNorm_int",
    "LogNorm",
    "Cox_int",
    "RSF",
    "KM"
  )
  present_learners = unique(dt$learner_id)
  learner_order = learner_order[learner_order %in% present_learners]
  dt[, learner_id := factor(learner_id, levels = learner_order)]

  # Colors
  cols = c(
    `Truth (true S, interp f)` = "#000000",
    Oracle        = "#dfce0cff",
    LLogis        = "#FFA500",
    Weibull       = "#756BB1",
    LogNorm_scale = "#2171B5",
    LogNorm_int   = "#6BAED6",
    LogNorm       = "#C6DBEF",
    Cox_int       = "#238B45",
    RSF           = "#D62728",
    KM            = "#7F7F7F"
  )
  cols = cols[names(cols) %in% present_learners]

  # Plot
  ggplot(dt, aes(x = .data[[y]], y = diff, color = learner_id)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = cols) +
    facet_grid(
      . ~ facet_label,
      labeller = label_parsed
    ) +
    labs(
      title = sprintf(
        "%s censoring task (Test set size = %d)",
        paste0(toupper(substring(task_id_val, 1, 1)), substring(task_id_val, 2)),
        n_test_val
      ),
      x = ifelse(y == "mise", "MISE (distance to true S)", "MIAE (distance to true S)"),
      y = expression(L[pred] - L[true]),
      color = "Model"
    ) +
    theme_bw(base_size = 14, base_family = "Arial") +
    theme(
      strip.background = element_rect(fill = "grey95"),
      strip.text = element_text(face = "bold"),
      legend.position = if (legend_position == "none") "none" else legend_position,
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12)
    ) +
    guides(color = guide_legend(override.aes = list(size = 5, alpha = 1)))
}

res = readRDS("results/rcll_sim.rds")

# n = 500 results
p1 = plot_rcll_diff_by_grid(res, task_id_val = "low", n_test_val = 500, legend_position = "none")
# ggsave(filename = "sensitivity_RCLL_grid1.png", dpi = 600, height = 8, width = 10, units = "in")
p2 = plot_rcll_diff_by_grid(res, task_id_val = "high", n_test_val = 500, legend_position = "bottom") +
  guides(color = guide_legend(override.aes = list(size = 5), nrow = 1, byrow = TRUE)) +
  theme(legend.direction = "horizontal")
# ggsave(filename = "sensitivity_RCLL_grid2.png", dpi = 600, height = 8, width = 10, units = "in")
p1 / p2

# n = 250 results
p1 = plot_rcll_diff_by_grid(res, task_id_val = "low", n_test_val = 250, legend_position = "none")
p2 = plot_rcll_diff_by_grid(res, task_id_val = "high", n_test_val = 250, legend_position = "bottom") +
  guides(color = guide_legend(override.aes = list(size = 5, alpha = 1), nrow = 1, byrow = TRUE)) +
  theme(legend.direction = "horizontal")
p1 / p2
ggsave(filename = "sensitivity_RCLL_grid.png", dpi = 600, height = 8, width = 13, units = "in")
ggsave(filename = "sensitivity_RCLL_grid.pdf", dpi = 600, height = 8, width = 13, units = "in", device = cairo_pdf)

# n = 100 results
p1 = plot_rcll_diff_by_grid(res, task_id_val = "low", n_test_val = 100)
p2 = plot_rcll_diff_by_grid(res, task_id_val = "high", n_test_val = 100)
p1 / p2

# n = 50 results
p1 = plot_rcll_diff_by_grid(res, task_id_val = "low", n_test_val = 50)
p2 = plot_rcll_diff_by_grid(res, task_id_val = "high", n_test_val = 50)
p1 / p2

# reduce learners to focus on the ones near the true model
learner_ids = c("CoxPH", "CoxPH_int", "LogNorm_noint_noshape",
                "LogNorm_noint_shape_x2", "LogNorm_int_noshape", "RSF", "KM")
# n = 500 results
p1 = plot_rcll_diff_by_grid(res, task_id_val = "low", n_test_val = 500, drop_learners = learner_ids)
p2 = plot_rcll_diff_by_grid(res, task_id_val = "high", n_test_val = 500, drop_learners = learner_ids)
p1 / p2
ggsave(filename = "sensitivity_RCLL_grid_few_learners.png", dpi = 600, height = 6, width = 12, units = "in")

# n = 250 results
p1 = plot_rcll_diff_by_grid(res, task_id_val = "low", n_test_val = 250,
  drop_learners = learner_ids, legend_position = "none")
p2 = plot_rcll_diff_by_grid(res, task_id_val = "high", n_test_val = 250,
  drop_learners = learner_ids, legend_position = "bottom") +
  guides(color = guide_legend(override.aes = list(size = 5), nrow = 1, byrow = TRUE)) +
  theme(legend.direction = "horizontal")
p1 / p2
ggsave(filename = "sensitivity_RCLL_grid_few_lrns.pdf", dpi = 600, height = 8,
 width = 13, units = "in", device = cairo_pdf)

# TABLE C2
# some RCLL stats
res_expanded = rbindlist(list(
  res,
  res[learner_id == "true",
      .(learner_id = "Truth_trueSf",   # new label
        mise,
        miae,
        rcll = true_rcll,              # use exact likelihood here
        true_rcll,
        task_id,
        rsmp_id,
        prop,
        n_times,
        n_test)]
), fill = TRUE)

# Rename existing "true" to reflect interpolation-based version
res_expanded[learner_id == "true", learner_id := "Truth_interpF"]

res_q = res_expanded[task_id == "low" & n_test == 250,
  .(
    mean = mean(rcll, na.rm = TRUE),
    sd = sd(rcll, na.rm = TRUE),
    w25 = quantile(rcll, 0.25, na.rm = TRUE),
    w75 = quantile(rcll, 0.75, na.rm = TRUE)
  ),
  by = .(learner_id, n_times)
]

res_q[n_times %in% c(16, 38, 76, 152, 379, 758) & learner_id == "Truth_trueSf"]
