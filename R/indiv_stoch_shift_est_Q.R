#' Estimate the Outcome Mechanism
#'
#' @details Compute the outcome regression for the observed data, including
#'  with the shift imposed by the intervention. This returns the outcome
#'  regression for the observed data (at A) and under the counterfactual shift
#'  shift (at A + delta).
#'
#' @param exposure A \code{character} vector of exposures to be shifted.
#' @param covars A \code{character} vector covariates to adjust for.
#' @param delta A \code{numeric} indicating the magnitude of the shift to be
#'  computed for the exposure \code{A}. This is passed to the internal
#'  \code{\link{shift_additive}} and is currently limited to additive shifts.
#' @param mu_learner Object containing a set of instantiated learners from the
#'  \pkg{sl3}, to be used in fitting an ensemble model.
#' @param av A \code{dataframe} of validation data specific to the fold
#' @param at A \code{dataframe} of training data specific to the fold
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

indiv_stoch_shift_est_Q <- function(exposure,
                                    delta,
                                    mu_learner,
                                    covars,
                                    av,
                                    at) {
  future::plan(future::sequential, gc = TRUE)

  # scale the outcome for logit transform
  y_star_av <- scale_to_unit(vals = av$y)
  y_star_at <- scale_to_unit(vals = at$y)

  av$y <- y_star_av
  at$y <- y_star_at

  # need a data set with the exposure stochastically shifted DOWNWARDS A-delta
  # do this for both AV and AT as AT is used in mediation

  av_downshifted <- data.table::copy(av)

  data.table::set(av_downshifted,
    j = exposure,
    value = shift_additive(
      a = subset(av, select = exposure), delta = -delta
    )
  )

  at_downshifted <- data.table::copy(at)
  data.table::set(at_downshifted,
    j = exposure,
    value = shift_additive(
      a = subset(at, select = exposure), delta = -delta
    )
  )

  # need a data set with the exposure stochastically shifted UPWARDS A+delta
  # do this for both AV and AT as AT is used in mediation

  av_upshifted <- data.table::copy(av)
  data.table::set(av_upshifted,
    j = exposure,
    value = shift_additive(
      a = subset(av, select = exposure), delta = delta
    )
  )

  at_upshifted <- data.table::copy(at)
  data.table::set(at_upshifted,
    j = exposure,
    value = shift_additive(
      a = subset(at, select = exposure), delta = delta
    )
  )

  # need a data set with the exposure stochastically shifted UPWARDS A+2delta
  # do this for both AV and AT as AT is used in mediation

  av_upupshifted <- data.table::copy(av)
  data.table::set(av_upupshifted,
    j = exposure,
    value = shift_additive(
      a = subset(av, select = exposure), delta = 2 * delta
    )
  )

  at_upupshifted <- data.table::copy(at)
  data.table::set(at_upupshifted,
    j = exposure,
    value = shift_additive(
      a = subset(at, select = exposure), delta = 2 * delta
    )
  )

  # Outcome mechanism
  sl <- sl3::Lrnr_sl$new(
    learners = mu_learner,
    metalearner = sl3::Lrnr_nnls$new()
  )

  at_task_noshift <- sl3::sl3_Task$new(
    data = at,
    covariates = covars,
    outcome = "y",
    outcome_type = "quasibinomial"
  )

  av_task_noshift <- sl3::sl3_Task$new(
    data = av,
    covariates = covars,
    outcome = "y",
    outcome_type = "quasibinomial"
  )

  at_task_upshift <- sl3::sl3_Task$new(
    data = at_upshifted,
    covariates = covars,
    outcome = "y",
    outcome_type = "quasibinomial"
  )

  av_task_upshift <- sl3::sl3_Task$new(
    data = av_upshifted,
    covariates = covars,
    outcome = "y",
    outcome_type = "quasibinomial"
  )

  at_task_upupshift <- sl3::sl3_Task$new(
    data = at_upupshifted,
    covariates = covars,
    outcome = "y",
    outcome_type = "quasibinomial"
  )

  av_task_upupshift <- sl3::sl3_Task$new(
    data = av_upupshifted,
    covariates = covars,
    outcome = "y",
    outcome_type = "quasibinomial"
  )

  at_task_downshift <- sl3::sl3_Task$new(
    data = at_downshifted,
    covariates = covars,
    outcome = "y",
    outcome_type = "quasibinomial"
  )

  av_task_downshift <- sl3::sl3_Task$new(
    data = av_downshifted,
    covariates = covars,
    outcome = "y",
    outcome_type = "quasibinomial"
  )

  sl_fit <- sl$train(at_task_noshift)

  # fit new Super Learner to the natural (no shift) data and predict
  at_q_pred <- bound_precision(sl_fit$predict(at_task_noshift))
  av_q_pred <- bound_precision(sl_fit$predict(av_task_noshift))

  # predict with Super Learner from unshifted data on the shifted data
  at_pred_upshifted <- bound_precision(sl_fit$predict(at_task_upshift))
  av_pred_upshifted <- bound_precision(sl_fit$predict(av_task_upshift))

  # predict with Super Learner from unshifted data on the shifted data
  at_pred_upupshifted <- bound_precision(sl_fit$predict(at_task_upupshift))
  av_pred_upupshifted <- bound_precision(sl_fit$predict(av_task_upupshift))

  # predict with Super Learner from unshifted data on the shifted data
  at_pred_downshifted <- bound_precision(sl_fit$predict(at_task_downshift))
  av_pred_downshifted <- bound_precision(sl_fit$predict(av_task_downshift))


  # create output data frame and return result
  av_out <- data.table::as.data.table(cbind(
    av_q_pred,
    av_pred_upshifted,
    av_pred_upupshifted,
    av_pred_downshifted
  ))
  data.table::setnames(av_out, c("noshift", "upshift", "upupshift", "downshift"))

  # create output data frame and return result
  at_out <- data.table::as.data.table(cbind(
    at_q_pred,
    at_pred_upshifted,
    at_pred_upupshifted,
    at_pred_downshifted
  ))
  data.table::setnames(at_out, c("noshift", "upshift", "upupshift", "downshift"))

  out_list <- list("q_at" = at_out, "q_av" = av_out)

  return(out_list)
}
