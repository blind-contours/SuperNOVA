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
#' @export
#' @return A \code{data.table} with two columns, containing estimates of the
#'  outcome mechanism at the natural value of the exposure Q(A, W) and an
#'  upshift of the exposure Q(A + delta, W).

est_Q_w_shifted_joint_mediation <- function(exposures,
                                            mediator,
                                            deltas,
                                            mu_learner,
                                            covars,
                                            av,
                                            at,
                                            zn_estims) {
  future::plan(future::sequential, gc = TRUE)

  # scale the outcome for logit transform
  y_star_av <- scale_to_unit(vals = av$y)
  y_star_at <- scale_to_unit(vals = at$y)

  av$y <- y_star_av
  at$y <- y_star_at

  at <- as.data.frame(at)
  av <- as.data.frame(av)

  # need a data set with the exposure stochastically shifted DOWNWARDS A-delta
  # do this for both AV and AT as AT is used in mediation

  outcome_results_list <- list()

  for (i in 1:length(exposures)) {
    exposure <- exposures[[i]]
    delta <- deltas[[i]]
    zn_estim <- zn_estims[[i]]

    # need a data set with the exposure stochastically shifted DOWNWARDS A-delta
    av_z_downshifted <- av_downshifted <- av
    av_downshifted[exposure] <- av_downshifted[exposure] - delta

    av_z_downshifted[exposure] <- av_z_downshifted[exposure] - delta
    av_z_downshifted[mediator] <- zn_estim$q_av$downshift

    # need a data set with the exposure stochastically shifted UPWARDS A+delta
    av_z_upshifted <- av_upshifted <- av
    av_upshifted[exposure] <- av_upshifted[exposure] + delta

    av_z_upshifted[exposure] <- av_z_upshifted[exposure] + delta
    av_z_upshifted[mediator] <- zn_estim$q_av$upshift

    # need a data set with the exposure stochastically shifted UPWARDS A+2delta
    av_z_upupshifted <- av_upupshifted <- av
    av_upupshifted[exposure] <- av_upupshifted[exposure] + 2 * delta

    av_z_upupshifted[exposure] <- av_z_upupshifted[exposure] + 2 * delta
    av_z_upupshifted[mediator] <- zn_estim$q_av$upupshift
    # at

    # need a data set with the exposure stochastically shifted DOWNWARDS A-delta
    at_z_downshifted <- at_downshifted <- at
    at_downshifted[exposure] <- at_downshifted[exposure] - delta

    at_z_downshifted[exposure] <- at_z_downshifted[exposure] - delta
    at_z_downshifted[mediator] <- zn_estim$q_at$downshift

    # need a data set with the exposure stochastically shifted UPWARDS A+delta
    at_z_upshifted <- at_upshifted <- at
    at_upshifted[exposure] <- at_upshifted[exposure] + delta

    at_z_upshifted[exposure] <- at_z_upshifted[exposure] + delta
    at_z_upshifted[mediator] <- zn_estim$q_at$upshift

    # need a data set with the exposure stochastically shifted UPWARDS A+2delta
    at_z_upupshifted <- at_upupshifted <- at
    at_upupshifted[exposure] <- at_upupshifted[exposure] + 2 * delta

    at_z_upupshifted[exposure] <- at_z_upupshifted[exposure] + 2 * delta
    at_z_upupshifted[mediator] <- zn_estim$q_at$upupshift


    # Outcome mechanism
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

    av_task_a_upshift <- sl3::sl3_Task$new(
      data = av_upshifted,
      covariates = covars,
      outcome = "y",
      outcome_type = "quasibinomial"
    )

    av_task_a_upupshift <- sl3::sl3_Task$new(
      data = av_upupshifted,
      covariates = covars,
      outcome = "y",
      outcome_type = "quasibinomial"
    )

    av_task_a_downshift <- sl3::sl3_Task$new(
      data = av_downshifted,
      covariates = covars,
      outcome = "y",
      outcome_type = "quasibinomial"
    )

    av_task_a_z_upshift <- sl3::sl3_Task$new(
      data = av_z_upshifted,
      covariates = covars,
      outcome = "y",
      outcome_type = "quasibinomial"
    )

    av_task_a_z_upupshift <- sl3::sl3_Task$new(
      data = av_z_upupshifted,
      covariates = covars,
      outcome = "y",
      outcome_type = "quasibinomial"
    )

    av_task_a_z_downshift <- sl3::sl3_Task$new(
      data = av_z_downshifted,
      covariates = covars,
      outcome = "y",
      outcome_type = "quasibinomial"
    )

    sl_fit <- sl$train(at_task_noshift)

    # fit new Super Learner to the natural (no shift) data and predict
    at_q_pred <- bound_precision(sl_fit$predict(at_task_noshift))

    av_q_pred <- bound_precision(sl_fit$predict(av_task_noshift))
    av_pred_a_upshift <- bound_precision(sl_fit$predict(av_task_a_upshift))
    av_pred_a_downshift <- bound_precision(sl_fit$predict(av_task_a_downshift))
    av_pred_a_upupshift <- bound_precision(sl_fit$predict(av_task_a_upupshift))

    at_pred_a_z_upshift <- bound_precision(sl_fit$predict(av_task_a_z_upshift))
    at_pred_a_z_upupshift <- bound_precision(sl_fit$predict(av_task_a_z_upupshift))
    at_pred_a_z_downshift <- bound_precision(sl_fit$predict(av_task_a_z_downshift))

    # create output data frame and return result
    av_a_shifted <- data.table::as.data.table(cbind(
      av_q_pred,
      av_pred_a_upshift,
      av_pred_a_upupshift,
      av_pred_a_downshift
    ))

    av_a_z_shifted <- data.table::as.data.table(cbind(
      av_q_pred,
      at_pred_a_z_upshift,
      at_pred_a_z_upupshift,
      at_pred_a_z_downshift
    ))

    data.table::setnames(av_a_shifted, c("noshift", "upshift", "upupshift", "downshift"))
    data.table::setnames(av_a_z_shifted, c("noshift", "upshift", "upupshift", "downshift"))

    outcome_results_list[[i]] <- list("a_shifted" = av_a_shifted, "a_z_shifted" = av_a_z_shifted)
  }

  return(outcome_results_list)
}
