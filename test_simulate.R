# Example use of the simulate function to generate non-PH data and 
# evaluate predictions from a correctly specified model 
# (Log AFT parametric) and a misspecified model (Cox PH).
source("simulate.R")

# simulate non-PH data for training
set.seed(42)
sim = simulate(n = 1000, b0 = 1.65, b1 = 0.5, b2 = -0.5, b12 = 2, sigma = c(0.5, 1.5))
  
task = sim$task
task$cens_prop() # ~ 42%
task$prop_haz() # very non-PH!
task$kaplan(strata = "x2") |> plot()

library(mlr3extralearners)
library(survdistr)
# true model
m1 = lrn("surv.flexreg", id = "LogNorm_int_shape_x2",
         formula = survival::Surv(time, status) ~ x1 + x2 + x1:x2,
         anc = list(sdlog = ~ x2),
         dist = "lognormal",
         times = times) # get survival predictions exactly on the `times`
# misspecified model
m2 = lrn("surv.coxph")

# train models
m1$train(task)
m2$train(task)

# generate test data (same DGP)
times = seq(0, 9.99, length.out = 500) # select specific grid for prediction times
#times = sort(runif(1000, min = 0, max = 10)) # random prediction times
test_sim = simulate(n = 50, b0 = 1.65, b1 = 0.5, b2 = -0.5, b12 = 2,
                    sigma = c(0.5, 1.5), times = times)
test_task = test_sim$task
true_pred = test_sim$pred
# get true survival matrix
true_surv = true_pred$data$distr
# true_surv[1:3, 1:3] # example S_Y(t) values

# predictions
p1 = m1$predict(test_task)
p2 = m2$predict(test_task)

# measure distance between predictions (interpolate to match grid of time points)
p1_surv = survdistr::mat_interp(x = p1$data$distr, eval_times = times, constant = TRUE)
p2_surv = survdistr::mat_interp(x = p2$data$distr, eval_times = times, constant = TRUE)
# check: same time points for S(t)
stopifnot(all(dim(true_surv) == dim(p1_surv)))
stopifnot(all(dim(p1_surv) == dim(p2_surv)))
# trick: force granularity of predicted S(t)
p1$data$distr = p1_surv
p2$data$distr = p2_surv

mise = function(S_true, S_pred, times) {
  dt_vec = diff(times)
  sq_diff = (S_true - S_pred)^2
  mean(rowSums(0.5 * (sq_diff[, -ncol(sq_diff)] + sq_diff[, -1]) * dt_vec))
}
mise(true_surv, p1_surv, times)
mise(true_surv, p2_surv, times)

# alternative S(t) distance:
miae = function(S_true, S_pred, times) {
  dt_vec = diff(times)
  abs_diff = abs(S_true - S_pred)
  mean(rowSums(0.5 * (abs_diff[, -ncol(abs_diff)] + abs_diff[, -1]) * dt_vec))
}

miae(true_surv, p1_surv, times)
miae(true_surv, p2_surv, times)

# C-index
true_pred$score() # highly biased due to large censoring
p1$score()
p2$score()

# RCLL
test_sim$rcll # true RCLL
rcll = msr("surv.rcll")
true_pred$score(rcll) # estimates f(t) from S(t) via linear interpolation
p1$score(rcll)
p2$score(rcll)

# RCLL*
source(file = "weighted_RCLL.R") # estimates both f from S and f_C from S_C (estimated itself from KM)
rcll2 = msr("surv.wrcll")
true_pred$score(rcll2, task = task)
p1$score(rcll2, task = task)
p2$score(rcll2, task = task)

# SBS (evaluated at a specific time point)
eval_time = 5
sbs = msr("surv.graf", integrated = FALSE, times = eval_time)
true_pred$score(sbs, task = task, train_set = task$row_ids)
p1$score(sbs, task = task, train_set = task$row_ids)
p2$score(sbs, task = task, train_set = task$row_ids)

# ISBS
q5 = quantile(times, 0.05)
q80 = quantile(times, 0.80)
eval_times = seq(q5, q80, length.out = 50) # 5th - 80th percentile of times
isbs = msr("surv.graf", times = eval_times) # specific times to integrate
# isbs = msr("surv.graf", p_max = 0.8) # test set times
true_pred$score(isbs, task = task, train_set = task$row_ids)
p1$score(isbs, task = task, train_set = task$row_ids)
p2$score(isbs, task = task, train_set = task$row_ids)

# MULTIPLE TEST SETS
B = 100 # number of test replications
n_test = 10 # n_test

#  predicted S(t) resolution (time grid) => MATTERS FOR ESTIMATION OF f(t)!
times = seq(0, 9.99, length.out = 1000)

scores = matrix(NA_real_, nrow = B, ncol = 8)
colnames(scores) = c("true", "true_S_est_f", "p1", "p2",
                     "p1_pred_grid_conS", "p2_pred_grid_conS",
                     "p1_pred_grid_linS", "p2_pred_grid_linS")

for (b in seq_len(B)) {
  message(b)
  test_sim = simulate(
    n = n_test,
    b0 = 1.65,
    b1 = 0.5,
    b2 = -0.5,
    b12 = 2,
    sigma = c(0.5, 1.5),
    times = times
  )

  test_task = test_sim$task
  true_pred = test_sim$pred

  p1 = m1$predict(test_task)
  p2 = m2$predict(test_task)

  scores[b, "true"] = test_sim$rcll
  scores[b, "true_S_est_f"] = true_pred$score(rcll)
  scores[b, "p1"] = p1$score(rcll)
  scores[b, "p2"] = p2$score(rcll)

  p1_distr = p1$data$distr
  p2_distr = p2$data$distr
  p1$data$distr = survdistr::mat_interp(x = p1_distr, eval_times = times, constant = TRUE)
  p2$data$distr = survdistr::mat_interp(x = p2_distr, eval_times = times, constant = TRUE)

  scores[b, "p1_pred_grid_conS"] = p1$score(rcll)
  scores[b, "p2_pred_grid_conS"] = p2$score(rcll)

  p1$data$distr = survdistr::mat_interp(x = p1_distr, eval_times = times, constant = FALSE)
  p2$data$distr = survdistr::mat_interp(x = p2_distr, eval_times = times, constant = FALSE)

  scores[b, "p1_pred_grid_linS"] = p1$score(rcll)
  scores[b, "p2_pred_grid_linS"] = p2$score(rcll)
}

sort(colMeans(scores))
