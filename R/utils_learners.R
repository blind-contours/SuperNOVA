#' Title
#'
#' @param n_obs
#' @param mu
#' @param sigma_mod
#'
#' @return
#' @export
#' @importFrom MASS mvrnorm
#' @importFrom bindata rmvbin

#' @examples
make_density_superlearner <- function() {
  make_metalearner <- function() {
    metalearner <- sl3::make_learner(
      Lrnr_nnls
    )
    return(metalearner)
  }

  make_continuous_superlearner <- function() {
    learners <- c(
      sl3::Lrnr_glm$new(),
      sl3::Lrnr_mean$new(),
      sl3::Lrnr_glmnet$new(alpha = 1),
      sl3::Lrnr_glmnet$new(alpha = 0),
      sl3::Lrnr_glmnet$new(alpha = .5),
      sl3::Lrnr_ranger$new(num.trees = 100),
      sl3::Lrnr_xgboost$new(),
      sl3::Lrnr_xgboost$new(nrounds = 50)
      # sl3::Lrnr_hal9001$new(max_degree = 2, reduce_basis = sqrt(1 / 1500))
    )


    learner_stack <- sl3::make_learner(sl3::Stack, learners)

    continuous_sl <- Lrnr_sl$new(
      learners = learner_stack,
      metalearner = make_metalearner()
    )

    return(continuous_sl)
  }

  fglm_lrnr <- Lrnr_glm_fast$new()
  # hal_lrnr <- Lrnr_hal9001$new(max_degree = 2, reduce_basis = sqrt(1 / 1500))
  # haldensify_lrnr <- Lrnr_haldensify$new(
  #   n_bins = c(10, 20),
  #   lambda_seq = exp(seq(-1, -10, length = 200))
  # )
  mean_lrnr <- make_continuous_superlearner()

  # semiparametric density estimator based on homoscedastic errors (HOSE)
  hose_lrnr <- Lrnr_density_semiparametric$new(mean_learner = mean_lrnr)

  # semiparametric density estimator based on heteroscedastic errors (HESE)
  hese_lrnr <- Lrnr_density_semiparametric$new(
    mean_learner = mean_lrnr,
    var_learner = make_continuous_superlearner()
  )

  sl_dens_lrnr <- Lrnr_sl$new(
    learners = list(hose_lrnr, hese_lrnr),
    metalearner = Lrnr_solnp_density$new()
  )

  return(sl_dens_lrnr)
}
