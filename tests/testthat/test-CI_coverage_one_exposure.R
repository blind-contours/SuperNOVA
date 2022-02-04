context("Make sure CI folds coverage and pooled results cover ground-truth for one exposure")
library(data.table)
library(sl3)
set.seed(172943)

# Example based on the data-generating mechanism presented in the simulation
n <- 500
W <- data.frame(W1 = runif(n), W2 = rbinom(n, 1, 0.7))
A <- rpois(n, lambda = exp(3 + .3 * log(W$W1) - 0.2 * exp(W$W1) * W$W2))
Y <- rbinom(
  n, 1,
  plogis(-1 + 0.05 * A - 0.02 * A * W$W2 + 0.2 * A * tan(W$W1^2) -
    0.02 * W$W1 * W$W2 + 0.1 * A * W$W1 * W$W2)
)
delta_shift <- 2

fitY.0 <- glm(
  Y ~ A + A:W2 + A:I(tan(W1^2)) + W1:W2 + A:W1:W2,
  family = binomial, data = data.frame(A, W)
)

Qn.0 <- function(A = A, W = W) {
  predict(
    fitY.0,
    newdata = data.frame(A, W, row.names = NULL),
    type = "response"
  )
}

Qn_ext_fitted <- as.data.table(
  lapply(c(0, delta_shift), function(shift_value) {
    Qn_out <- Qn.0(A = A + shift_value, W = W)
  })
)
setnames(Qn_ext_fitted, c("noshift", "upshift"))

sl_density_lrnr <- make_density_superlearner()

Lrnr_earth_1 <- Lrnr_earth$new(linpreds = FALSE, degree = 1)
Lrnr_earth_2 <- Lrnr_earth$new(linpreds = FALSE, degree = 2)
Lrnr_earth_3 <- Lrnr_earth$new(linpreds = FALSE, degree = 2, pmethod = "none")


learners <- c(
  Lrnr_earth_1,
  Lrnr_earth_2,
  Lrnr_earth_3
)

names(learners) <- c(
  "full earth 1",
  "full earth 2",
  "full earth 3"
)


Exposures_stack <- make_learner(Stack, learners)

Lrnr_glm_basic <- Lrnr_glm$new()
Lrnr_mean_base <- Lrnr_mean$new()

Lrnr_ridge <- Lrnr_glmnet$new(alpha = 0)
Lrnr_lasso <- Lrnr_glmnet$new(alpha = 1)

learners <- c(
  Lrnr_earth_1,
  Lrnr_earth_2,
  Lrnr_earth_3,
  Lrnr_glm_basic,
  Lrnr_mean_base,
  Lrnr_ridge,
  Lrnr_lasso
)

names(learners) <- c(
  "full earth 1",
  "full earth 2",
  "full earth 3",
  "Lrnr_glm",
  "Lrnr_mean",
  "Lrnr_ridge",
  "Lrnr_lasso"
)

Covariate_stack <- make_learner(Stack, learners)

mean_lrnr <- Lrnr_mean$new()
fglm_lrnr <- Lrnr_glm_fast$new()
rf_lrnr <- Lrnr_ranger$new()
lasso_learner <- Lrnr_glmnet$new(alpha = 1)
ridge_learner <- Lrnr_glmnet$new(alpha = 0)
lrn_polspline <- Lrnr_polspline$new()
lrn_ranger100 <- make_learner(Lrnr_ranger, num.trees = 100)
# hal_lrnr <- Lrnr_hal9001$new(max_degree = 3, n_folds = 3)

Outcome_stack <- make_learner(
  Stack, mean_lrnr, fglm_lrnr, rf_lrnr, lasso_learner, ridge_learner, lrn_polspline, lrn_ranger100
)

# fit TMLE
tmle_results <- SuperNOVA(
  W = W,
  V = NULL,
  A = A,
  Y = Y,
  delta = delta_shift,
  LOD_val = 0,
  estimator = "tmle",
  fluctuation = "standard",
  max_iter = 10,
  Density_stack = sl_density_lrnr,
  Exposures_stack = Exposures_stack,
  Covariate_stack = Covariate_stack,
  Outcome_stack = Outcome_stack,
  n_folds = 5,
  family = "continuous",
  quantile_thresh = 0,
  verbose = TRUE,
  parallel = TRUE
)

meta_results <- compute_meta_results(SuperNOVA_results = tmle_results, parameter = "Indiv Shift")

pooled_tmle_results <- meta_results$`Pooled Results`[[1]]
truth <- mean(Qn_ext_fitted$upshift - Qn_ext_fitted$noshift)

all_CI_between <- between(truth, pooled_tmle_results$`Lower CI`, pooled_tmle_results$`Upper CI`)
# test for reasonable equality between estimators
test_that("All CI bands contain truth", {
  expect_true(all(all_CI_between))
})
