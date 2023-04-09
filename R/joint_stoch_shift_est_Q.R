#' Estimate the Outcome Mechanism
#'
#' @details Compute the outcome regression for the observed data, including
#'  with the shift imposed by the intervention. This returns the outcome
#'  regression for the observed data (at A) and under the counterfactual shift
#'  shift (at A + delta).
#'
#' @param exposures A \code{character} vector of exposures to be shifted.
#' @param covars A \code{character} vector covariates to adjust for.
#' @param deltas A \code{numeric} indicating the magnitude of the shift to be
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
#' @export
#' @return A \code{data.table} with two columns, containing estimates of the
#'  outcome mechanism at the natural value of the exposure Q(A, W) and an
#'  upshift of the exposure Q(A + delta, W).

joint_stoch_shift_est_Q <- function(exposures,
                                    deltas,
                                    mu_learner,
                                    covars,
                                    av,
                                    at,
                                    upper_bound,
                                    lower_bound) {
  future::plan(future::sequential, gc = TRUE)

  # scale the outcome for logit transform
  y_star_av <- scale_to_unit(vals = av$y)
  y_star_at <- scale_to_unit(vals = at$y)

  av$y <- y_star_av
  at$y <- y_star_at

  results <- list()

  for (i in 1:length(exposures)) {
    exposure <- exposures[[i]]
    delta <- deltas[[i]]


    # need a data set with the exposure stochastically shifted DOWNWARDS A-delta
    av_downshifted <- data.table::copy(av)
    data.table::set(av_downshifted, j = exposure, value = shift_additive(
      a = subset(av, select = exposure),
      delta = -delta,
      upper_bound = upper_bound,
      lower_bound = lower_bound
    ))

    # need a data set with the exposure stochastically shifted UPWARDS A+delta
    av_upshifted <- data.table::copy(av)
    data.table::set(av_upshifted, j = exposure, value = shift_additive(
      a = subset(av, select = exposure),
      delta = delta,
      upper_bound = upper_bound,
      lower_bound = lower_bound
    ))

    # need a data set with the exposure stochastically shifted UPWARDS A+2delta
    av_upupshifted <- data.table::copy(av)
    data.table::set(av_upupshifted, j = exposure, value = shift_additive(
      a = subset(av, select = exposure),
      delta = 2 * delta,
      upper_bound = upper_bound,
      lower_bound = lower_bound
    ))

    sl <- Lrnr_sl$new(
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

    av_task_upshift <- sl3::sl3_Task$new(
      data = av_upshifted,
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

    av_task_downshift <- sl3::sl3_Task$new(
      data = av_downshifted,
      covariates = covars,
      outcome = "y",
      outcome_type = "quasibinomial"
    )

    sl_fit <- sl$train(at_task_noshift)

    # fit new Super Learner to the natural (no shift) data and predict
    av_pred_star_qn <- bound_precision(sl_fit$predict(av_task_noshift))
    # predict with Super Learner from unshifted data on the shifted data
    av_pred_star_qn_upshifted <- bound_precision(sl_fit$predict(av_task_upshift))
    # predict with Super Learner from unshifted data on the shifted data
    av_pred_star_qn_upupshifted <- bound_precision(sl_fit$predict(av_task_upupshift))
    av_pred_star_qn_downshifted <- bound_precision(sl_fit$predict(av_task_downshift))

    # create output data frame and return result
    out <- data.table::as.data.table(cbind(
      av_pred_star_qn,
      av_pred_star_qn_upshifted,
      av_pred_star_qn_upupshifted,
      av_pred_star_qn_downshifted
    ))

    data.table::setnames(out, c("noshift", "upshift", "upupshift", "downshift"))

    results[[i]] <- out
  }

  return(results)
}
