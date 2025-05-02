# Extra analyses + plots (that didn't make it to the main HTML report)
library(ggplot2)
library(dplyr)
library(tidyr)

# Data ----
res = readRDS("results/res_sims10000_distrs1000_0.rds")
# res = readRDS("results/res_sims10000_distrs1000_1.rds") # S constant interp

res = res |>
  select(!matches("shape|scale")) |> # remove columns
  mutate(cens_bin = cut(prop_cens, breaks = seq(0, 1, by = 0.2), include.lowest = TRUE)) |>
  mutate(tv_dist_bin = cut(tv_dist, breaks = seq(0, 1, by = 0.25), include.lowest = TRUE)) |>
  select(!c("prop_cens", "tv_dist"))
res |> count(n, name = "sim")

# Violation stats ----
# Define a violation threshold
threshold = 1e-04

all_stats = res |>
  group_by(n) |>
  summarize(
    total = n(),
    SBS_median_n_violations = sum(SBS_median_diff > threshold),
    SBS_median_violation_rate = mean(SBS_median_diff > threshold),
    SBS_median_diff_mean = mean(SBS_median_diff[SBS_median_diff > threshold]),
    SBS_median_se_mean = mean(SBS_median_se[SBS_median_diff > threshold]),
    SBS_q10_n_violations = sum(SBS_q10_diff > threshold),
    SBS_q10_violation_rate = mean(SBS_q10_diff > threshold),
    SBS_q10_diff_mean = mean(SBS_q10_diff[SBS_q10_diff > threshold]),
    SBS_q90_n_violations = sum(SBS_q90_diff > threshold),
    SBS_q90_violation_rate = mean(SBS_q90_diff > threshold),
    SBS_q90_diff_mean = mean(SBS_q90_diff[SBS_q90_diff > threshold]),
    ISBS_n_violations = sum(ISBS_diff > threshold),
    ISBS_violation_rate = mean(ISBS_diff > threshold),
    ISBS_diff_mean = mean(ISBS_diff[ISBS_diff > threshold]),
    RCLL_n_violations = sum(RCLL_diff > threshold),
    RCLL_violation_rate = mean(RCLL_diff > threshold),
    rRCLL_n_violations = sum(rRCLL_diff > threshold),
    rRCLL_violation_rate = mean(rRCLL_diff > threshold),
    .groups = "drop"
  )

# no (r)RCLL violations:
all_stats |> select(n | contains("RCLL"))
# SBS violations:
all_stats |> select(n | contains("SBS"))
all_stats |> select(n | contains("ISBS"))

# check se for violations ((I)SBS scores)
measures = c("SBS_median", "SBS_q10", "SBS_q90", "ISBS")

# Compute proportions per measure
se_res = lapply(measures, function(measure) {
  diff_col = paste0(measure, "_diff")
  se_col = paste0(measure, "_se")

  res |>
    filter(.data[[diff_col]] > threshold) |>
    mutate(se_diff = .data[[se_col]] >= .data[[diff_col]]) |>
    group_by(n) |>
    summarise(prop_se_gte_diff = mean(se_diff), .groups = "drop") |>
    mutate(measure = measure)
}) |>
bind_rows() |>
pivot_wider(names_from = measure, values_from = prop_se_gte_diff)

se_res

# look at SBS violations grouped by (n, %cens)
stats_cens =
  res |>
  group_by(n, cens_bin) |>
  summarize(
    total = n(),
    SBS_median_n_violations = sum(SBS_median_diff > threshold),
    SBS_median_violation_rate = mean(SBS_median_diff > threshold),
    SBS_median_diff_mean = mean(SBS_median_diff[SBS_median_diff > threshold]),
    SBS_q10_n_violations = sum(SBS_q10_diff > threshold),
    SBS_q10_violation_rate = mean(SBS_q10_diff > threshold),
    SBS_q10_diff_mean = mean(SBS_q10_diff[SBS_q10_diff > threshold]),
    SBS_q90_n_violations = sum(SBS_q90_diff > threshold),
    SBS_q90_violation_rate = mean(SBS_q90_diff > threshold),
    SBS_q90_diff_mean = mean(SBS_q90_diff[SBS_q90_diff > threshold]),
    ISBS_n_violations = sum(ISBS_diff > threshold),
    ISBS_violation_rate = mean(ISBS_diff > threshold),
    ISBS_diff_mean = mean(ISBS_diff[ISBS_diff > threshold]),
    .groups = "drop"
  )
stats_cens
# table(stats_cens$cens_bin)

# look at SBS violations grouped by (n, distance of (Y:true, Y:pred) distributions)
stats_tv_dist =
  res |>
  group_by(n, tv_dist_bin) |>
  summarize(
    total = n(),
    SBS_median_n_violations = sum(SBS_median_diff > threshold),
    SBS_median_violation_rate = mean(SBS_median_diff > threshold),
    SBS_median_diff_mean = mean(SBS_median_diff[SBS_median_diff > threshold]),
    SBS_q10_n_violations = sum(SBS_q10_diff > threshold),
    SBS_q10_violation_rate = mean(SBS_q10_diff > threshold),
    SBS_q10_diff_mean = mean(SBS_q10_diff[SBS_q10_diff > threshold]),
    SBS_q90_n_violations = sum(SBS_q90_diff > threshold),
    SBS_q90_violation_rate = mean(SBS_q90_diff > threshold),
    SBS_q90_diff_mean = mean(SBS_q90_diff[SBS_q90_diff > threshold]),
    ISBS_n_violations = sum(ISBS_diff > threshold),
    ISBS_violation_rate = mean(ISBS_diff > threshold),
    ISBS_diff_mean = mean(ISBS_diff[ISBS_diff > threshold]),
    .groups = "drop"
  )
stats_tv_dist

# (I)SBS: Line plots ----
## Cens Prop ----
stats_cens |>
  ggplot(aes(x = n, y = SBS_median_violation_rate, color = cens_bin, group = cens_bin)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.001, linetype = "dashed", color = "black") +
  annotate("text", x = min(stats_cens$n), y = 0.00115,
           label = "0.1% threshold = 1:1000",
           hjust = 0, vjust = -0.5, size = 3) +
  labs(
    title = "SBS Violation Rates by Sample Size and Censorship Proportion",
    subtitle = expression("" * tau^"*" == Q[50](T[obs]) * ""),
    x = "Sample Size (n)",
    y = "SBS Violation Rate (%)",
    color = "Censorship Bin"
  ) +
  scale_x_log10() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_minimal() +
  scale_color_brewer(palette = "Set1")

stats_cens |>
  ggplot(aes(x = n, y = SBS_q10_violation_rate, color = cens_bin, group = cens_bin)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.001, linetype = "dashed", color = "black") +
  annotate("text", x = min(stats_cens$n), y = 0.00115,
           label = "0.1% threshold = 1:1000",
           hjust = 0, vjust = -0.5, size = 3) +
  labs(
    title = "SBS Violation Rates by Sample Size and Censorship Proportion",
    subtitle = expression("" * tau^"*" == Q[10](T[obs]) * ""),
    x = "Sample Size (n)",
    y = "SBS Violation Rate (%)",
    color = "Censorship Bin"
  ) +
  scale_x_log10() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_minimal() +
  scale_color_brewer(palette = "Set1")

stats_cens |>
  ggplot(aes(x = n, y = SBS_q90_violation_rate, color = cens_bin, group = cens_bin)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.001, linetype = "dashed", color = "black") +
  annotate("text", x = min(stats_cens$n), y = 0.00115,
           label = "0.1% threshold = 1:1000",
           hjust = 0, vjust = -0.5, size = 3) +
  labs(
    title = "SBS Violation Rates by Sample Size and Censorship Proportion",
    subtitle = expression("" * tau^"*" == Q[90](T[obs]) * ""),
    x = "Sample Size (n)",
    y = "SBS Violation Rate (%)",
    color = "Censorship Bin"
  ) +
  scale_x_log10() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_minimal() +
  scale_color_brewer(palette = "Set1")

stats_cens |>
  ggplot(aes(x = n, y = ISBS_violation_rate, color = cens_bin, group = cens_bin)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.001, linetype = "dashed", color = "black") +
  annotate("text", x = min(stats_cens$n), y = 0.00115,
           label = "0.1% threshold = 1:1000",
           hjust = 0, vjust = -0.5, size = 3) +
  labs(
    title = "ISBS Violation Rates by Sample Size and Censorship Proportion",
    x = "Sample Size (n)",
    y = "ISBS Violation Rate (%)",
    color = "Censorship Bin"
  ) +
  scale_x_log10() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_minimal() +
  scale_color_brewer(palette = "Set1")

## TV distance ----
stats_tv_dist |>
  ggplot(aes(x = n, y = SBS_median_violation_rate, color = tv_dist_bin, group = tv_dist_bin)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.001, linetype = "dashed", color = "black") +
  annotate("text", x = min(stats_tv_dist$n), y = 0.00115,
           label = "0.1% threshold = 1:1000",
           hjust = 0, vjust = -0.5, size = 3) +
  labs(
    title = expression("SBS Violation Rates by Sample Size and Total Variational Distance (Y, "*hat(Y)*")"),
    subtitle = expression("" * tau^"*" == Q[50](T[obs]) * ""),
    x = "Sample Size (n)",
    y = "SBS Violation Rate (%)",
    color = "Total Variational\nDistance"
  ) +
  scale_x_log10() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_minimal() +
  scale_color_brewer(palette = "Set1")

stats_tv_dist |>
  ggplot(aes(x = n, y = ISBS_violation_rate, color = tv_dist_bin, group = tv_dist_bin)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.001, linetype = "dashed", color = "black") +
  annotate("text", x = min(stats_tv_dist$n), y = 0.00115,
           label = "0.1% threshold = 1:1000",
           hjust = 0, vjust = -0.5, size = 3) +
  labs(
    title = expression("ISBS Violation Rates by Sample Size and Total Variational Distance (Y, "*hat(Y)*")"),
    x = "Sample Size (n)",
    y = "ISBS Violation Rate (%)",
    color = "Total Variational\nDistance"
  ) +
  scale_x_log10() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_minimal() +
  scale_color_brewer(palette = "Set1")

# SBS: Barplots ----
## Cens Prop ----
stats_cens |>
  ggplot(aes(x = cens_bin, y = SBS_median_violation_rate, fill = cens_bin)) +
  geom_col() +
  facet_wrap(~ n, labeller = labeller(n = function(x) paste0("n = ", x))) +
  labs(
    title = "SBS Violation Rates by Censoring Proportion",
    subtitle = expression("" * tau^"*" == Q[50](T[obs]) * ""),
    x = "Censoring Proportion",
    y = "SBS Violation Rate (%)"
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_minimal() +
  scale_fill_brewer(palette = "Blues")

stats_cens |>
  ggplot(aes(x = cens_bin, y = ISBS_violation_rate, fill = cens_bin)) +
  geom_col() +
  facet_wrap(~ n, labeller = labeller(n = function(x) paste0("n = ", x))) +
  labs(
    title = "ISBS Violation Rates by Censoring Proportion",
    subtitle = expression("" * tau^"*" == Q[10](T[obs]) * ""),
    x = "Censoring Proportion",
    y = "ISBS Violation Rate (%)"
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_minimal() +
  scale_fill_brewer(palette = "Blues")

## TV distance ----
stats_tv_dist |>
  ggplot(aes(x = tv_dist_bin, y = SBS_median_violation_rate, fill = tv_dist_bin)) +
  geom_col() +
  facet_wrap(~ n, ncol = 2, labeller = labeller(n = function(x) paste0("n = ", x))) +
  labs(
    title = expression("SBS Violation Rates by Total Variational Distance (Y, "*hat(Y)*")"),
    subtitle = expression("" * tau^"*" == Q[50](T[obs]) * ""),
    x = "Censoring Proportion",
    y = "SBS Violation Rate (%)",
    fill = "Total Variational\nDistance"
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_minimal() +
  scale_fill_brewer(palette = "Greens")
