#' Targeted Minimum Loss Estimate of Counterfactual Mean of Stochastic Shift
#' Intervention
#'
#' @details Invokes the procedure to construct a targeted minimum loss estimate
#'  (TMLE) of the counterfactual mean under a modified treatment policy.
#'
#' @param data_internal A \code{data.table} constructed internally by a call to
#'  expo_shift. This contains most of the data for computing the
#'  targeted minimum loss (TML) estimator.
#' @param delta A \code{numeric} value indicating the shift in the treatment to
#'  be used in defining the target parameter. This is defined with respect to
#'  the scale of the treatment (A).
#' @param Qn_scaled An object providing the value of the outcome evaluated after
#'  imposing a shift in the treatment. This object is passed in after being
#'  constructed by a call to the internal function for Q estimation.
#' @param Hn An object providing values of the auxiliary ("clever")
#'  covariate, constructed from the treatment mechanism and required for
#'  targeted minimum loss-based estimation. This object object should be passed
#'  in after being constructed by a call to clever covariate construction.
#' @param fluctuation The method to be used in the submodel fluctuation step
#'  (targeting step) to compute the TML estimator. The choices are "standard"
#'  and "weighted" for where to place the auxiliary covariate in the logistic
#'  tilting regression.
#' @param eif_reg_type Whether a flexible nonparametric function ought to be
#'  used in the dimension-reduced nuisance regression of the targeting step for
#'  the censored data case. By default, the method used is a nonparametric
#'  regression based on the Highly Adaptive Lasso.
#'  Set this to \code{"glm"} to instead use a simple linear regression model.
#'  In this step, the efficient influence function (EIF) is regressed against
#'  covariates contributing to the censoring mechanism (i.e., EIF ~ V | C = 1).
#' @param y Outcome variable
#'
#' @importFrom assertthat assert_that
#' @importFrom data.table as.data.table setnames
#' @importFrom stringr str_detect
#' @importFrom Rdpack reprompt
#' @export
#'
#' @return S3 object of class \code{txshift} containing the results of the
#'  procedure to compute a TML estimate of the treatment shift parameter.
tmle_exposhift <- function(data_internal,
                           delta,
                           Qn_scaled,
                           Qn_unscaled,
                           Hn,
                           fluctuation = c("standard", "weighted"),
                           eif_reg_type = c("hal", "glm"),
                           y,
                           estimator = "tmle"
                           ) {
  # initialize counter
  n_steps <- 0

  # fit logistic regression to fluctuate along the sub-model
  fitted_fluc_mod <- fit_fluctuation(
    y = y,
    Qn_scaled = Qn_scaled,
    Hn = Hn,
    method = fluctuation
  )

  # compute TML estimate and EIF for the treatment shift parameter
  tmle_eif_out <- eif(
    y = y,
    qn = Qn_scaled,
    qn_unscaled = Qn_unscaled,
    hn = Hn,
    estimator = estimator,
    fluc_mod_out = fitted_fluc_mod,
    data
  )

  # create output object
  exposure_shift_out <- unlist(
    list(
      psi = tmle_eif_out[["psi"]],
      var = tmle_eif_out[["var"]],
      se = tmle_eif_out[["se"]],
      ci = tmle_eif_out[["CI"]],
      p_value = tmle_eif_out[["p_value"]],
      eif = list(tmle_eif_out[["eif"]]),
      .iter_res = NULL,
      n_iter = n_steps,
      estimator = estimator,
      .outcome = list(data_internal$y),
      .delta = list(delta),
      qn_shift_star = list(fitted_fluc_mod$qn_shift_star),
      qn_noshift_star = list(fitted_fluc_mod$qn_noshift_star),
      noshift_psi = tmle_eif_out[["no shift psi"]],
      noshift_var = tmle_eif_out[["no shift var"]],
      noshift_se = tmle_eif_out[["no shift se"]],
      noshift_CI = tmle_eif_out[["no shift CI"]],
      p_value = tmle_eif_out[["p_value"]],
      noshift_eif = list(tmle_eif_out[["no shift eif"]])
    ),
    recursive = FALSE
  )


  # S3-ify and return output object
  class(exposure_shift_out) <- "SuperNOVA"
  return(exposure_shift_out)
}

###############################################################################

#' Fit One-Dimensional Fluctuation Model for Updating Initial Estimates
#'
#' @details Procedure for fitting a one-dimensional fluctuation model to update
#'  the initial estimates of the outcome regression based on the auxiliary
#'  covariate. These updated estimates are subsequently used to construct the
#'  TML estimator of the counterfactual mean under a modified treatment policy.
#'
#' @param y A \code{numeric} vector corresponding to an outcome variable.
#' @param Qn_scaled An object providing the value of the outcome evaluate
#'  after inducing a shift in the exposure. This object should be passed in
#'  after being constructed by the Q estimators.
#' @param Hn An object providing values of the auxiliary ("clever") covariate,
#'  constructed from the treatment mechanism and required for targeted minimum
#'  loss estimation. This object object should be passed in after being
#'  constructed by a call to clever covariate construction
#' @param ipc_weights A \code{numeric} vector that gives inverse probability of
#'  censoring weights for each observation. These are generated by invoking the
#'  routines for estimating the censoring mechanism.
#' @param method A \code{character} giving the type of regression to be used in
#'  traversing the fluctuation sub-model. The available choices are "weighted"
#'  and "standard". Consult the literature for details on the differences.
#' @param flucmod_tol A \code{numeric} indicating the largest value to be
#'  tolerated in the fluctuation model for the targeted minimum loss estimator.
#'
#' @importFrom stats qlogis glm fitted predict as.formula coef
#' @importFrom data.table as.data.table setnames
#' @export
#' @return A \code{list} containing the fluctuation model (a \code{glm} object)
#'  produced by logistic regression, a \code{character} vector indicating the
#'  type of fluctuation (whether the auxiliary covariates was used as a weight
#'  or included directly in the model formula), the updated estimates of the
#'  outcome regression under the shifted value of the exposure, and the updated
#'  estimates of the outcome regression under the natural value of exposure.
fit_fluctuation <- function(y,
                            Qn_scaled,
                            Hn,
                            ipc_weights = rep(1, length(y)),
                            method = c("standard", "weighted"),
                            flucmod_tol = 50) {
  y_star <- scale_to_unit(
    vals = y
  )

  # bound precision for use of logit transform
  Qn_scaled_bounded <- data.table::as.data.table(apply(
    Qn_scaled, 2,
    bound_precision
  ))

  # extract Q and obtain logit transform
  Qn_noshift_logit <- stats::qlogis(Qn_scaled_bounded$noshift)

  # fit the fluctuation regression in one of two ways
  if (method == "standard") {
    # note that \epsilon_n will be the coefficient of the covariate Hn
    suppressWarnings(
      fluctuation_model <- stats::glm(
        formula = stats::as.formula(
          "y_star ~ -1 + offset(logit_Qn) + Hn"
        ),
        data = data.table::as.data.table(list(
          y_star = y_star,
          logit_Qn = Qn_noshift_logit,
          Hn = Hn$noshift
        )),
        weights = ipc_weights,
        family = "quasibinomial",
        start = 0
      )
    )
    coefs_fluc <- stats::coef(fluctuation_model)

    # check convergence of fluctuation model and sanity of estimates
    if (!fluctuation_model$converged || abs(max(coefs_fluc)) > flucmod_tol) {
      suppressWarnings(
        fluctuation_model <- stats::glm(
          formula = stats::as.formula("y_star ~ -1 + offset(logit_Qn) + Hn"),
          data = data.table::as.data.table(list(
            y_star = y_star,
            logit_Qn = Qn_noshift_logit,
            Hn = Hn$noshift
          )),
          weights = ipc_weights,
          family = "binomial"
        )
      )
      # if the fluctuation model hasn't converged or is unstable, simply set
      # the coefficients to disable updating, i.e., coef(Hn) := 0
      if (!fluctuation_model$converged | abs(max(coefs_fluc)) > flucmod_tol) {
        fluctuation_model$coefficients <- 0
      }
    }
  } else if (method == "weighted") {
    # note that epsilon_n will be the intercept term here
    suppressWarnings(
      fluctuation_model <- stats::glm(
        formula = stats::as.formula("y_star ~ offset(logit_Qn)"),
        data = data.table::as.data.table(list(
          y_star = y_star,
          logit_Qn = Qn_noshift_logit
        )),
        weights = as.numeric(Hn$noshift * ipc_weights),
        family = "binomial",
        start = 0
      )
    )
    coefs_fluc <- stats::coef(fluctuation_model)

    # check covergence of fluctuation model and sanity of estimates
    if (!fluctuation_model$converged || abs(max(coefs_fluc)) > flucmod_tol) {
      suppressWarnings(
        fluctuation_model <- stats::glm(
          formula = stats::as.formula("y_star ~ offset(logit_Qn)"),
          data = data.table::as.data.table(list(
            y_star = y_star,
            logit_Qn = Qn_noshift_logit
          )),
          weights = as.numeric(Hn$noshift * ipc_weights),
          family = "binomial"
        )
      )
      # if the updated fluctuation model hasn't converged or is unstable,
      # simply set the coefficient to zero to disable updating
      if (!fluctuation_model$converged | abs(max(coefs_fluc)) > flucmod_tol) {
        fluctuation_model$coefficients <- 0
      }
    }
  }

  # get fitted values from fluctuation model
  Qn_noshift_star_unit <- unname(stats::fitted(fluctuation_model))
  Qn_noshift_star <- scale_to_original(
    scaled_vals = Qn_noshift_star_unit,
    max_orig = max(y),
    min_orig = min(y)
  )

  # need to logit transform Qn(d(A,W),W)
  Qn_shift_logit <- stats::qlogis(Qn_scaled_bounded$upshift)

  # get Qn_star for the SHIFTED data
  if (method == "standard") {
    Qn_shift_star_data <- data.table::as.data.table(list(
      logit_Qn = Qn_shift_logit,
      Hn = Hn$shift
    ))

    # predict from fluctuation model on Q(d(A,W),W) and Hn(d(A,W))
    Qn_shift_star_unit <- unname(stats::predict(
      object = fluctuation_model,
      newdata = Qn_shift_star_data,
      type = "response"
    ))
  } else if (method == "weighted") {
    Qn_shift_star_data <- data.table::as.data.table(Qn_shift_logit)
    data.table::setnames(Qn_shift_star_data, "logit_Qn")

    # predict from fluctuation model on Q(d(A,W),W)
    Qn_shift_star_unit <- unname(stats::predict(
      object = fluctuation_model,
      newdata = Qn_shift_star_data,
      type = "response"
    ))
  }

  Qn_shift_star <- scale_to_original(
    scaled_vals = Qn_shift_star_unit,
    max_orig = max(y),
    min_orig = min(y)
  )

  # return the fit model object
  out <- list(
    fluc_fit = fluctuation_model,
    covar_method = method,
    qn_shift_star = as.numeric(Qn_shift_star),
    qn_noshift_star = as.numeric(Qn_noshift_star)
  )
  return(out)
}
