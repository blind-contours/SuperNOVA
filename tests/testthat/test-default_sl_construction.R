library(sl3)
# Test Script

# Create default Super Learner estimators
sls <- create_sls()

# Create a test set
test_data <- data.frame(
  x = runif(1000),
  y = runif(1000)
)

# Check that the correct number of estimators is created
expect_equal(length(sls), 6)

# Check that the estimators are of the correct type
expect_true(class(sls$pi_learner)[1] == "Lrnr_sl")
expect_true(class(sls$zeta_learner)[1] == "Stack")
expect_true(class(sls$mu_learner)[1] == "Stack")
expect_true(class(sls$e_learner)[1] == "Stack")
expect_true(class(sls$g_learner)[1] == "Stack")
expect_true(class(sls$em_learner)[1] == "Stack")
