# S(t)/f(t) ESTIMATION/INTERPOLATION FUNCTIONS (via Kaplan-Meier)

#' Linearly interpolate (and extrapolate) a Kaplan-Meier survival curve
#' @param fit survfit object
#' @param new_times vector of times (unordered, possibly duplicated)
#' @return interpolated S(t) values
interp_surv = function(fit, new_times) {
  # remove constant-interpolated values from Kaplan-Meier S(t) fit
  keep = !duplicated(fit$surv)
  surv = fit$surv[keep] # decreasing
  times = fit$time[keep] # ordered

  # Edge case: constant survival
  if (length(unique(surv)) == 1) {
    return(rep(surv[1], length(new_times)))
  }

  # Normal interpolation (at least two S(t) values here)
  int_surv = stats::approx(
    x = times, y = surv, xout = new_times,
    yleft = 1, method = "linear", rule = 2
  )$y

  # Extrapolate manually if needed
  min_time = min(times)
  max_time = max(times)

  # Precompute slopes for extrapolation
  slope_min = (surv[1L] - 1) / min_time
  slope_max = (surv[length(surv)] - surv[length(surv) - 1L]) / (max_time - times[length(times) - 1L])

  idx_min = new_times < min_time
  idx_max = new_times > max_time

  # Linear extrapolation considering that S(t = 0) = 1
  if (any(idx_min) && surv[1L] < 1) {
    int_surv[idx_min] = 1 + slope_min * new_times[idx_min]
  }

  # Linear extrapolation using the last time interval
  if (any(idx_max) && surv[length(surv)] > 0) {
    extrapolated_value = surv[length(surv)] + slope_max * (new_times[idx_max] - max_time)
    int_surv[idx_max] = pmax(0, extrapolated_value) # force S >= 0
  }

  int_surv
}

#' PDF estimation from Kaplan-Meier survival curve
#' @param fit survfit object
#' @param new_times numeric vector (unordered, duplicated allowed)
#' @return numeric vector of density values f(t)
interp_pdf = function(fit, new_times) {
  # keep all unique sorted times (KM and requested) for pdf
  utimes = sort(unique(c(fit$time, new_times)))

  # Create a mapping of `new_times` to `utimes`
  indx = match(new_times, utimes)

  # Interpolated survival function
  surv = interp_surv(fit, utimes)

  # CDF = 1 - S
  cdf = 1 - surv

  # Numerical derivative: f = dF/dt
  dt = diff(utimes)
  dF = diff(cdf)

  # Density (finite difference)
  dens = dF / dt

  # For timepoints exactly at utimes, align left
  dens_full = c(dens[1], dens) # replicate first slope for first point

  # return density at `new_times`, clip any negatives to 0
  pmax(dens_full[indx], 0)
}
