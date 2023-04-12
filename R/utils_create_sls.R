#' @title Create default Super Learner estimators for the data adaptive
#' and nuisance parameters used in `SuperNOVA`
#' @description If Super Learners are not passed to the stack arguments
#' this function is used to create some default ensemble machine
#' learning estimators for each parameter. The default estimators are
#' fast but also flexible.
#'
#' @import sl3
#' @return List of ensemble estimators
#' @export

create_sls <- function() {
  # Create default pi estimator ---------------------------
  make_metalearner <- function() {
    metalearner <- sl3::make_learner(
      Lrnr_nnls
    )
    return(metalearner)
  }

  make_continuous_superlearner <- function() {
    learners <- c(
      sl3::Lrnr_glm$new(),
      # sl3::Lrnr_glmnet$new(alpha = 1),
      # sl3::Lrnr_glmnet$new(alpha = 0),
      # sl3::Lrnr_glmnet$new(alpha = .5),
      # sl3::Lrnr_ranger$new(num.trees = 100),
      # sl3::Lrnr_ranger$new(num.trees = 500),
      # sl3::Lrnr_xgboost$new(nrounds = 50),
      # sl3::Lrnr_xgboost$new(nrounds = 200),
      sl3::Lrnr_earth$new(degree = 2),
      sl3::Lrnr_earth$new(degree = 3)
    )


    learner_stack <- sl3::make_learner(sl3::Stack, learners)

    continuous_sl <- Lrnr_sl$new(
      learners = learner_stack,
      metalearner = make_metalearner()
    )

    return(continuous_sl)
  }

  mean_lrnr <- make_continuous_superlearner()


  # semiparametric density estimator based on homoscedastic errors (HOSE)
  hose_lrnr <- Lrnr_density_semiparametric$new(mean_learner = mean_lrnr)

  # semiparametric density estimator based on heteroscedastic errors (HESE)
  hese_lrnr <- Lrnr_density_semiparametric$new(
    mean_learner = mean_lrnr,
    var_learner = make_continuous_superlearner()
  )

  pi_learner <- Lrnr_sl$new(
    learners = list(hose_lrnr, hese_lrnr),
    metalearner = Lrnr_solnp_density$new()
  )

  # Create default zeta estimator ---------------------------

  lrnr_earth_1 <- Lrnr_earth$new(linpreds = FALSE, degree = 1)
  lrnr_earth_2 <- Lrnr_earth$new(linpreds = FALSE, degree = 2)
  lrnr_earth_3 <- Lrnr_earth$new(
    linpreds = FALSE, degree = 2,
    pmethod = "none"
  )
  lrnr_earth_4 <- Lrnr_earth$new(
    linpreds = FALSE, degree = 2,
    thresh = 0.0005
  )
  lrnr_earth_5 <- Lrnr_earth$new(linpreds = FALSE, degree = 2)
  lrnr_earth_6 <- Lrnr_earth$new(
    linpreds = FALSE, degree = 2,
    fast.k = 10
  )
  lrnr_earth_7 <- Lrnr_earth$new(
    linpreds = FALSE, degree = 2,
    fast.k = 5
  )
  lrnr_earth_8 <- Lrnr_earth$new(
    linpreds = FALSE, degree = 2,
    fast.beta = 0
  )
  lrnr_earth_9 <- Lrnr_earth$new(
    linpreds = FALSE, degree = 2,
    fast.k = 10, pmethod = "cv", nfold = 5
  )
  lrnr_earth_10 <- Lrnr_earth$new(
    linpreds = TRUE, degree = 1
  )

  lrnr_earth_11 <- Lrnr_earth$new(
    linpreds = TRUE, degree = 2,
    fast.k = 0
  )

  lrnr_earth_12 <- Lrnr_earth$new(
    degree = 3
  )

  learners <- c(
    lrnr_earth_1,
    lrnr_earth_2,
    lrnr_earth_3,
    lrnr_earth_4,
    lrnr_earth_5,
    lrnr_earth_6,
    lrnr_earth_7,
    lrnr_earth_8,
    lrnr_earth_9,
    lrnr_earth_10,
    lrnr_earth_11,
    lrnr_earth_12
  )

  names(learners) <- c(
    "full earth 1",
    "full earth 2",
    "full earth 3",
    "full earth 4",
    "full earth 5",
    "full earth 6",
    "full earth 7",
    "full earth 8",
    "full earth 9",
    "full earth 10",
    "full earth 11",
    "full earth 12"
  )

  zeta_learner <- make_learner(Stack, learners)

  # Create default mu estimator ---------------------------

  lrnr_glm_basic <- Lrnr_glm$new()
  lrnr_ridge <- Lrnr_glmnet$new(alpha = 0)
  # lrnr_lasso <- Lrnr_glmnet$new(alpha = 1)
  # lrnr_ranger_100 <- make_learner(Lrnr_ranger, num.trees = 100)
  # lrnr_xgboost_df <- make_learner(Lrnr_xgboost)
  # lrnr_xgboost_50 <- make_learner(Lrnr_xgboost, nrounds = 50)
  # lrnr_xgboost_100 <- make_learner(Lrnr_xgboost, nrounds = 100)
  # lrnr_xgboost_200 <- make_learner(Lrnr_xgboost, nrounds = 200)
  # lrnr_xgboost_300 <- make_learner(Lrnr_xgboost, nrounds = 300)

  learners <- c(
    lrnr_glm_basic,
    # lrnr_ridge
    # lrnr_lasso,
    # lrnr_ranger_100,
    lrnr_xgboost_df,
    # lrnr_xgboost_50,
    # lrnr_xgboost_100,
    # lrnr_xgboost_200,
    # lrnr_xgboost_300
  )


  mu_learner <- make_learner(Stack, learners)
  e_learner <- make_learner(Stack, learners)
  g_learner <- make_learner(Stack, learners)


  ## categorical exposure learners
  # lrnr_glm_multivar <- Lrnr_multivariate$new()
  lrnr_ridge <- Lrnr_glmnet$new(alpha = 0)
  lrnr_lasso <- Lrnr_glmnet$new(alpha = 1)
  lrnr_ranger_100 <- make_learner(Lrnr_ranger, num.trees = 100)
  lrnr_xgboost_df <- make_learner(Lrnr_xgboost)
  lrnr_xgboost_50 <- make_learner(Lrnr_xgboost, nrounds = 50)
  lrnr_xgboost_100 <- make_learner(Lrnr_xgboost, nrounds = 100)
  lrnr_xgboost_200 <- make_learner(Lrnr_xgboost, nrounds = 200)
  lrnr_xgboost_300 <- make_learner(Lrnr_xgboost, nrounds = 300)
  lrnr_polspline <- Lrnr_polspline$new()

  learners <- c(
    lrnr_ranger_100,
    lrnr_polspline,
    lrnr_xgboost_100
  )

  quant_stack <- make_learner(Stack, learners)
  discrete_sl_metalrn <- sl3::Lrnr_cv_selector$new(sl3::loss_loglik_multinomial)

  quant_lrnr <- sl3::Lrnr_sl$new(
    learners = quant_stack,
    metalearner = discrete_sl_metalrn,
  )

  return(list(
    "pi_learner" = pi_learner,
    "zeta_learner" = zeta_learner,
    "mu_learner" = mu_learner,
    "e_learner" = e_learner,
    "g_learner" = g_learner,
    "quant_learner" = quant_lrnr
  ))
}
