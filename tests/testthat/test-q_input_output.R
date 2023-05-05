library(data.table)

# Create dummy data
exposure <- c("exposure1")
delta <- 0.5
mu_learner <- c("lrn_glm")
covars <- c("cov1", "cov2")

av <- data.table(
  "exposure1" = rnorm(100),
  "cov1" = rnorm(100),
  "cov2" = rnorm(100),
  "y" = rnorm(100)
)

at <- data.table(
  "exposure1" = rnorm(100),
  "cov1" = rnorm(100),
  "cov2" = rnorm(100),
  "y" = rnorm(100)
)

# Test function
mu_learner <- make_learner(Stack, c(Lrnr_glm$new(), Lrnr_mean$new()))

out_list <- indiv_stoch_shift_est_Q(
  exposure,
  delta,
  mu_learner = mu_learner,
  covars = c(covars, exposure),
  av,
  at,
  lower_bound = -Inf,
  upper_bound = Inf
)

# Test that the function returns a list with two elements
if (length(out_list) != 2) {
  stop("Function did not return a list with two elements")
}

# Test that the list elements are data tables
if (class(out_list$q_at)[1] != "data.table") {
  stop("q_at element of list is not a data table")
}

if (class(out_list$q_av)[1] != "data.table") {
  stop("q_av element of list is not a data table")
}

# Test that the data tables have the correct number of columns
if (ncol(out_list$q_at) != 4) {
  stop("q_at data table does not have four columns")
}

if (ncol(out_list$q_av) != 4) {
  stop("q_av data table does not have four columns")
}
