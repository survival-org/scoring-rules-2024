# plotting the true survival curves by x2 group for the two tasks
library(survival)
library(ggplot2)
library(patchwork)

tasks = readRDS("simulation_benchmark/tasks.rds")

plot_km_task = function(task, title, show_censor = FALSE, extend_to_last_time = TRUE) {
  d = task$data()

  # summary quantities
  cens_total = mean(d$status == 0)
  cens_x20 = mean(d$status == 0 & d$x2 == 0)
  cens_x21 = mean(d$status == 0 & d$x2 == 1)
  n_events = sum(d$status == 1)
  admin_cens = mean(d$time == 10)
  #prop_admin = admin_cens / cens_total

  subtitle =
    paste0(
      sprintf("Total censoring: %.1f%% (x2 = 0: %.1f%%, x2 = 1: %.1f%%)\n",
        100 * cens_total, 100 * cens_x20, 100 * cens_x21
      ),
      sprintf("Number of events: %d\n", n_events),
      sprintf("Administrative censoring: %.1f%%", 100*admin_cens)
    )

  # Kaplan-Meier fit
  fit = survival::survfit(
    Surv(time, status) ~ factor(x2),
    data = d
  )

  s = summary(fit)
  df = data.frame(
    time = s$time,
    surv = s$surv,
    lower = s$lower,
    upper = s$upper,
    strata = s$strata,
    n_censor = s$n.censor
  )

  # ensure S(0) = 1 for each strata
  by_strata = split(df, df$strata)
  s0_rows = lapply(by_strata, function(dd) {
    if (any(dd$time == 0)) {
      NULL
    } else {
      data.frame(
        time = 0,
        surv = 1,
        lower = 1,
        upper = 1,
        strata = dd$strata[1],
        n_censor = 0
      )
    }
  })
  s0_rows = do.call(rbind, s0_rows)
  if (!is.null(s0_rows)) {
    df = rbind(s0_rows, df)
  }

  if (extend_to_last_time) {
    max_time = max(d$time)
    by_strata = split(df, df$strata)
    ext = lapply(by_strata, function(dd) {
      last_row = dd[nrow(dd), , drop = FALSE]
      if (max_time > last_row$time) {
        last_row$time = max_time
        last_row
      } else {
        NULL
      }
    })
    ext = do.call(rbind, ext)
    if (!is.null(ext)) {
      df = rbind(df, ext)
    }
  }

  censor_df = df[df$n_censor > 0, , drop = FALSE]

  p =
    ggplot(df, aes(x = time, y = surv, color = strata, linetype = strata)) +
    geom_ribbon(
      aes(ymin = lower, ymax = upper, group = strata),
      fill = "grey80",
      alpha = 0.4,
      color = NA
    ) +
    geom_step() +
    scale_color_manual(
      values = c("black", "red"),
      labels = c(
        expression(x[2] == 0 ~ paste("(", sigma, "  = ", 0.5, ")")),
        expression(x[2] == 1 ~ paste("(", sigma, "  = ", 1.5, ")"))
      )
    ) +
    scale_linetype_manual(
      values = c(1, 1),
      labels = c(
        expression(x[2] == 0 ~ paste("(", sigma, "  = ", 0.5, ")")),
        expression(x[2] == 1 ~ paste("(", sigma, "  = ", 1.5, ")"))
      )
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Time",
      y = "Survival probability",
      color = NULL,
      linetype = NULL
    ) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = c(0.98, 0.98),
      legend.justification = c(1, 1),
      legend.background = element_rect(fill = "white", color = "grey80"),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 16),
      legend.key.size = grid::unit(1, "lines")
    )

  if (show_censor && nrow(censor_df) > 0) {
    p = p + geom_point(
      data = censor_df,
      aes(x = time, y = surv, color = strata),
      shape = 3,
      stroke = 0.9,
      size = 2,
      show.legend = FALSE
    )
  }

  p
}

p_low  = plot_km_task(tasks$low, title = "Low censoring task", show_censor = TRUE)
p_low
p_high = plot_km_task(tasks$high, title = "High censoring task", show_censor = TRUE)
p_high

(p_low | p_high)
ggsave("kaplan_meier_plots_DPG.png", device = png, dpi = 600, width = 11, height = 6)
