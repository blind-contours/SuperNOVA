# Load necessary packages
library(devtools)
library(dplyr)
library(magrittr)
library(data.table)

# Load SuperNOVA functions
load_all()

# Test function for SuperNOVA
test_that("Ensure initial estimates from Q are pretty close to truth in discrete a case", {
  # Generate simulated data
  full_data <- simulate_mediation_data(n_obs = 100000, delta = 1)
  total_effect_truth <- full_data$ate_a1_quant
  p0_data <- full_data$data

  # Set seed for reproducibility
  seed <- sample(1:10000, 1)
  set.seed(seed)

  # Sample data
  sample_size <- 1000
  data_sim <- p0_data %>%
    dplyr::slice_sample(n = sample_size)

  # Prepare input data
  w <- data_sim[, c("w_1", "w_2")]
  a <- data_sim[, "a_1_quant"]
  z <- data_sim[, "z_1_quant"]
  y <- data_sim[, "y_quant"]

  # Create data tables
  at <- data.table(w, a = a, z = z, y = y)
  av <- data.table(w, a = a, z = z, y = y)

  # Test est_Q_w_shifted_mediation function
  mu_learner <- create_sls()$mu_learner

  qn_estim_no_med <- est_Q_w_shifted_mediation(
    exposure = "a",
    mediator = "z",
    delta = 1,
    mu_learner = mu_learner,
    covars = c("w_1", "w_2", "a"),
    av = av,
    at = at,
    upper_bound = max(a),
    lower_bound = min(a)
  )

  upshift_diff <- mean(qn_estim_no_med$av_predictions$upshift) - mean(y)

  # Check that total_effect_truth and upshift_diff are not very different
  tolerance <- 1
  expect_equal(upshift_diff, total_effect_truth, tolerance = tolerance)
})
