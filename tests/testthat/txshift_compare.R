library(txshift)

data_sample <- data[sample(nrow(data), 2000), ]

w <- data_sample[, c("W1", "W2")]
a <- data_sample[, c("M1", "M2", "M3")]
y <- data_sample$Y

txshift_results <- txshift(
  W = w,
  A = a[, 1],
  Y = y,
  delta = 1,
  estimator = "tmle",
  g_exp_fit_args = list(
    fit_type = "sl",
    sl_learners_density = Lrnr_density_hse$new(Lrnr_hal9001$new())
  ),
  Q_fit_args = list(fit_type = "glm", glm_formula = "Y ~ .^2")
)

txshift_results$psi - mean(data_sample$Y)
