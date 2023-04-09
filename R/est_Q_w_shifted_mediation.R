#' Estimate the Outcome Mechanism with Shifted A and Z
#'
#' @details Compute the outcome regression for the observed data, including
#'  with the shift imposed by the intervention. This returns the outcome
#'  regression for the observed data (at A) and under the counterfactual shift
#'  shift (at A + delta).
#'
#' @param exposure A \code{character} vector of exposures to be shifted.
#' @param mediator The mediator variable
#' @param covars A \code{character} vector covariates to adjust for.
#' @param at The training data
#' @param delta A \code{numeric} indicating the magnitude of the shift to be
#'  computed for the exposure \code{A}. This is passed to the internal
#'  \code{\link{shift_additive}} and is currently limited to additive shifts.
#' @param mu_learner Object containing a set of instantiated learners from the
#'  \pkg{sl3}, to be used in fitting an ensemble model.
#' @param av A \code{dataframe} of validation data specific to the fold
#' @param at A \code{dataframe} of training data specific to the fold
#' @param zn_estim Density estimates of the mediator under various shifts
#'
#' @importFrom stats glm as.formula predict
#' @importFrom data.table as.data.table setnames copy set
#' @importFrom stringr str_detect
#' @importFrom assertthat assert_that
#' @import sl3
#' @export
#' @return A \code{data.table} with two columns, containing estimates of the
#'  outcome mechanism at the natural value of the exposure Q(A, W) and an
#'  upshift of the exposure Q(A + delta, W).

est_Q_w_shifted_mediation <- function(exposure,
                                      mediator,
                                      delta,
                                      mu_learner,
                                      covars,
                                      av,
                                      at,
                                      upper_bound = upper_bound,
                                      lower_bound = lower_bound) {
  future::plan(future::sequential, gc = TRUE)

  # # scale the outcome for logit transform
  # y_star_av <- scale_to_unit(vals = av$y)
  # y_star_at <- scale_to_unit(vals = at$y)
  #
  # av$y <- y_star_av
  # at$y <- y_star_at

  # need a data set with the exposure stochastically shifted DOWNWARDS A-delta
  # do this for both AV and AT as AT is used in mediation

  av_a_upshifted <- data.table::copy(av)
  at_a_upshifted <- data.table::copy(at)

  data.table::set(av_a_upshifted,
    j = exposure,
    value = shift_additive(
      a = subset(av, select = exposure),
      delta = delta,
      lower_bound = lower_bound,
      upper_bound = upper_bound
    )
  )

  data.table::set(at_a_upshifted,
    j = exposure,
    value = shift_additive(
      a = subset(at, select = exposure),
      delta = delta,
      lower_bound = lower_bound,
      upper_bound = upper_bound
    )
  )

  # need a data set with the exposure stochastically shifted DOWNWARDS A-delta
  av_a_downshifted <- data.table::copy(av)
  data.table::set(av_a_downshifted, j = exposure, value = shift_additive(
    a = subset(av, select = exposure),
    delta = -delta,
    lower_bound = lower_bound,
    upper_bound = upper_bound
  ))

  # need a data set with the exposure stochastically shifted UPWARDS A+2delta
  av_a_upupshifted <- data.table::copy(av)
  data.table::set(av_a_upupshifted, j = exposure, value = shift_additive(
    a = subset(av, select = exposure),
    delta = 2 * delta,
    lower_bound = lower_bound,
    upper_bound = upper_bound
  ))


  # Outcome mechanism
  sl <- sl3::Lrnr_sl$new(
    learners = mu_learner,
    metalearner = make_learner(sl3::Lrnr_nnls)
  )

  at_task_noshift <- sl3::sl3_Task$new(
    data = at,
    covariates = covars,
    outcome = "y",
    outcome_type = "continuous"
  )

  at_task_upshift <- sl3::sl3_Task$new(
    data = at_a_upshifted,
    covariates = covars,
    outcome = "y",
    outcome_type = "continuous"
  )

  av_task_noshift <- sl3::sl3_Task$new(
    data = av,
    covariates = covars,
    outcome = "y",
    outcome_type = "continuous"
  )

  av_task_a_upshift <- sl3::sl3_Task$new(
    data = av_a_upshifted,
    covariates = covars,
    outcome = "y",
    outcome_type = "continuous"
  )

  av_task_a_upupshift <- sl3::sl3_Task$new(
    data = av_a_upupshifted,
    covariates = covars,
    outcome = "y",
    outcome_type = "continuous"
  )

  av_task_a_downshift <- sl3::sl3_Task$new(
    data = av_a_downshifted,
    covariates = covars,
    outcome = "y",
    outcome_type = "continuous"
  )

  sl_fit <- sl$train(at_task_noshift)

  # fit new Super Learner to the natural (no shift) data and predict
  at_q_pred <- sl_fit$predict(at_task_noshift)

  av_q_pred <- sl_fit$predict(av_task_noshift)
  av_pred_a_upshift <- sl_fit$predict(av_task_a_upshift)
  av_pred_a_downshift <- sl_fit$predict(av_task_a_downshift)
  av_pred_a_upupshift <- sl_fit$predict(av_task_a_upupshift)
  at_pred_a_upshift <- sl_fit$predict(at_task_upshift)
  at_pred_a_noshift <- sl_fit$predict()


  # create output data frame and return result
  av_a_shifted <- data.table::as.data.table(cbind(
    av_q_pred,
    av_pred_a_upshift,
    av_pred_a_upupshift,
    av_pred_a_downshift
  ))
  data.table::setnames(av_a_shifted, c("noshift", "upshift", "upupshift", "downshift"))

  # create output data frame and return result
  at_a_shifted <- data.table::as.data.table(cbind(
    at_pred_a_noshift,
    at_pred_a_upshift
  ))
  data.table::setnames(at_a_shifted, c("noshift", "upshift"))


  return(list("av_predictions" = av_a_shifted, "at_predictions" = at_a_shifted, "model" = sl_fit))
}
