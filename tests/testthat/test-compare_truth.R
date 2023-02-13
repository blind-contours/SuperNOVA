# context("Ensure truth for moderate sample size is covered by SuperNOVA' CIs ")
# library(dplyr)
#
# data_shift1 <- simulate_data(n_obs = 2000, shift_var_index = 1)
# data_sample <- data_shift1$data
# effect <- data_shift1$effect
#
# w <- data_sample[, c("W1", "W2")]
# a <- data_sample[, c("M1", "M2", "M3")]
# y <- data_sample$Y
#
# deltas <- list("M1" = 1, "M2" = 0, "M3" = 0)
#
# ptm <- proc.time()
#
# sim_results <- SuperNOVA(
#   w = w,
#   a = a,
#   y = y,
#   delta = deltas,
#   n_folds = 2,
#   num_cores = 6,
#   family = "continuous",
#   quantile_thresh = 0,
#   seed = 294580
# )
#
# proc.time() - ptm
#
# exp_1_results <- sim_results$`Indiv Shift Results`$M1
#
# cov_ind <- between(effect, exp_1_results$`Lower CI`[nrow(exp_1_results)],
#                   exp_1_results$`Upper CI`[nrow(exp_1_results)])
#
# test_that("Resulting CIs include ground truth under one shift",
#           {expect_equal(TRUE, cov_ind)})
#
#
