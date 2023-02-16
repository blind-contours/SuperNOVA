library(data.table)
library(sl3)
library(tmle3)
devtools::install_github("nhejazi/medshift@fix-tmle3")
library(medshift)

# produces a simple data set based on ca causal model with mediation
make_simple_mediation_data <- function(n_obs = 1000) {
  # baseline covariate -- simple, binary
  W <- rbinom(n_obs, 1, prob = 0.50)

  # create treatment based on baseline W
  A <- as.numeric(rbinom(n_obs, 1, prob = W / 4 + 0.1))

  # single mediator to affect the outcome
  z1_prob <- 1 - plogis((A^2 + W) / (A + W^3 + 0.5))
  Z <- rbinom(n_obs, 1, prob = z1_prob)

  # create outcome as a linear function of A, W + white noise
  Y <- Z + A - 0.1 * W + rnorm(n_obs, mean = 0, sd = 0.25)

  # full data structure
  data <- as.data.table(cbind(Y, Z, A, W))
  setnames(data, c("Y", "Z", "A", "W"))
  return(data)
}

# set seed and simulate example data
set.seed(75681)
example_data <- make_simple_mediation_data()

# instantiate learners
g_learner <- e_learner <-
  Lrnr_cv$new(Lrnr_glm$new(family = stats::binomial()), full_fit = TRUE)
m_learner <- phi_learner <-
  Lrnr_cv$new(Lrnr_glm$new(family = stats::gaussian()), full_fit = TRUE)

# Test different estimators
reweighted_medshift <- medshift(
  W = example_data$W, A = example_data$A,
  Z = example_data$Z, Y = example_data$Y,
  delta = 3, estimator = "reweighted",
  estimator_args = list(cv_folds = 3)
)

onestep_medshift <- medshift(
  W = example_data$W, A = example_data$A,
  Z = example_data$Z, Y = example_data$Y,
  delta = 3, estimator = "onestep",
  estimator_args = list(cv_folds = 3)
)

sub_medshift <- medshift(
  W = example_data$W, A = example_data$A,
  Z = example_data$Z, Y = example_data$Y,
  delta = 3, estimator = "substitution",
  estimator_args = list(cv_folds = 3)
)

# get 'tmle3_Task' not found error here
tmle_medshift <- medshift(
  W = example_data$W, A = example_data$A,
  Z = example_data$Z, Y = example_data$Y,
  delta = 3, estimator = "tmle",
  g_learners = g_learner, e_learner = e_learner,
  m_learner = m_learner, phi_learners = phi_learner,
  estimator_args = list(
    cv_folds = 3,
    max_iter = 10,
    step_size = 1e-6
  )
)
