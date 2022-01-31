context("Make sure CI folds coverage and pooled results cover ground-truth for two exposures")
library(data.table)
set.seed(172943)

# Example based on the data-generating mechanism presented in the simulation
n <- 1200
W <- data.frame(W1 = runif(n), W2 = rbinom(n, 1, 0.7))
A_1 <- rpois(n, lambda = exp(3 + .3 * log(W$W1) - 0.2 * exp(W$W1) * W$W2))
A_2 <- rpois(n, lambda = exp(4 + .1 * log(W$W1) - 0.2 * exp(W$W1) * W$W2))

A <- data.frame(A_1 = A_1, A_2 = A_2)
Y <- 0.01 *(A_1 * A_2 - 0.01 * W$W1 + -0.04 * W$W2)

delta_shift <- 2

fitY.0 <- glm(
  Y ~  A_1*A_2 + W1*W2,
  family = gaussian, data = data.frame(A, W)
)

Qn.0 <- function(A_1 = A_1, A_2 = A_2,  W = W) {
  predict(
    fitY.0,
    newdata = data.frame(A_1, A_2, W, row.names = NULL),
    type = "response"
  )
}

shift_preds <-  function(shift_value) {
    no_shift_Qn_out <- Qn.0(A_1 = A_1, A_2 = A_2, W = W)
    A1_Qn_out <- Qn.0(A_1 = A_1 + shift_value, A_2 = A_2, W = W)
    A2_Qn_out <- Qn.0(A_1 = A_1,  A_2 = A_2 + shift_value, W = W)
    A1A2_Qn_out <- Qn.0(A_1 = A_1 + shift_value,  A_2 = A_2 + shift_value, W = W)

    results <- as.data.table(cbind(no_shift_Qn_out, A1_Qn_out, A2_Qn_out, A1A2_Qn_out ))

    return(results)
}

shift_results <- shift_preds(delta_shift)

setnames(shift_results, c("noshift", "A1 upshift", "A2 upshift", "A1 A2 upshift"))

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
  Stack, mean_lrnr, fglm_lrnr, rf_lrnr, lasso_learner, ridge_learner, lrn_polspline, lrn_ranger100, Lrnr_earth_1,Lrnr_earth_2, Lrnr_earth_3
)

# fit TMLE
tmle_results <- SuperNOVA(W = W,
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
                          n_folds = 4,
                          family = "continuous",
                          quantile_thresh = 0,
                          verbose = TRUE,
                          parallel = TRUE)

meta_results <- compute_meta_results(SuperNOVA_results = tmle_results, parameter = "Joint Shift")

psi_tmle_results <- meta_results$`Pooled Results`$`Psi-A_1A_2`
A1_tmle_results <- meta_results$`Pooled Results`$`A_1-A_1A_2`
A2_tmle_results <- meta_results$`Pooled Results`$`A_2-A_1A_2`
joint_tmle_results <- meta_results$`Pooled Results`$`A_1A_2-A_1A_2`
no_shift_tmle_results <- meta_results$`Pooled Results`$`No Shift-A_1A_2`

psi_truth <- mean(shift_results$`A1 A2 upshift`) - mean(shift_results$`A1 upshift`) - mean(shift_results$`A2 upshift`)
A1_shift_truth <- mean(shift_results$`A1 upshift`)
A2_shift_truth <- mean(shift_results$`A2 upshift`)
joint_shift_truth <- mean(shift_results$`A1 A2 upshift`)
no_shift_truth <- mean(shift_results$noshift)

A1_pooled_CI_between <- between(A1_shift_truth, A1_tmle_results$`Lower CI`[nrow(A1_tmle_results)] , A1_tmle_results$`Upper CI`[nrow(A1_tmle_results)])
A2_pooled_CI_between <- between(A2_shift_truth, A2_tmle_results$`Lower CI`[nrow(A2_tmle_results)] , A2_tmle_results$`Upper CI`[nrow(A2_tmle_results)])
joint_pooled_CI_between <- between(joint_shift_truth, joint_tmle_results$`Lower CI`[nrow(joint_tmle_results)] , joint_tmle_results$`Upper CI`[nrow(joint_tmle_results)])

# test for reasonable equality between estimators
test_that("All CI bands contain truth", {
  expect_true(all(all_CI_between))
})


