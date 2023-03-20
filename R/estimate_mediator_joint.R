#' Estimate the Mediator Mechanism for Joint Shifts
#'
#' @details Compute the outcome regression for the observed data, including
#'  with the shift imposed by the intervention. This returns the mediator
#'  regression for the observed data (at A) and under the counterfactual shift
#'  shift (at A + delta).
#'
#' @param exposure A \code{character} vector of exposures to be shifted.
#' @param mediator A \code{character} vector of mediator outcome
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
#' @export
#' @return A \code{data.table} with two columns, containing estimates of the
#'  outcome mechanism at the natural value of the exposure Q(A, W) and an
#'  upshift of the exposure Q(A + delta, W).

estimate_mediator_joint <- function(mediator,
                                    exposures,
                                    deltas,
                                    mu_learner,
                                    covars,
                                    av,
                                    at) {
  future::plan(future::sequential, gc = TRUE)

  av <- as.data.frame(av)
  at <- as.data.frame(at)

  med_outcome_list <- list()

  for (i in 1:length(exposures)) {
    exposure <- exposures[[i]]
    delta <- deltas[[i]]

    # need a data set with the exposure stochastically shifted DOWNWARDS A-delta
    av_downshifted <- av
    av_downshifted[exposure] <- av_downshifted[exposure] - delta

    # need a data set with the exposure stochastically shifted UPWARDS A+delta
    av_upshifted <- av
    av_upshifted[exposure] <- av_upshifted[exposure] + delta

    # need a data set with the exposure stochastically shifted UPWARDS A+2delta
    av_upupshifted <- av
    av_upupshifted[exposure] <- av_upupshifted[exposure] + 2 * delta

    # at

    # need a data set with the exposure stochastically shifted DOWNWARDS A-delta
    at_downshifted <- at
    at_downshifted[exposure] <- at_downshifted[exposure] - delta

    # need a data set with the exposure stochastically shifted UPWARDS A+delta
    at_upshifted <- at
    at_upshifted[exposure] <- at_upshifted[exposure] + delta

    # need a data set with the exposure stochastically shifted UPWARDS A+2delta
    at_upupshifted <- at
    at_upupshifted[exposure] <- at_upupshifted[exposure] + 2 * delta


    # Outcome mechanism
    sl <- Lrnr_sl$new(
      learners = mu_learner,
      metalearner = sl3::Lrnr_nnls$new()
    )

    at_task_noshift <- sl3::sl3_Task$new(
      data = at,
      covariates = covars,
      outcome = mediator,
      outcome_type = "continuous"
    )

    av_task_noshift <- sl3::sl3_Task$new(
      data = av,
      covariates = covars,
      outcome = mediator,
      outcome_type = "continuous"
    )

    at_task_upshift <- sl3::sl3_Task$new(
      data = at_upshifted,
      covariates = covars,
      outcome = mediator,
      outcome_type = "continuous"
    )

    av_task_upshift <- sl3::sl3_Task$new(
      data = av_upshifted,
      covariates = covars,
      outcome = mediator,
      outcome_type = "continuous"
    )

    at_task_upupshift <- sl3::sl3_Task$new(
      data = at_upupshifted,
      covariates = covars,
      outcome = mediator,
      outcome_type = "continuous"
    )

    av_task_upupshift <- sl3::sl3_Task$new(
      data = av_upupshifted,
      covariates = covars,
      outcome = mediator,
      outcome_type = "continuous"
    )

    at_task_downshift <- sl3::sl3_Task$new(
      data = at_downshifted,
      covariates = covars,
      outcome = mediator,
      outcome_type = "continuous"
    )

    av_task_downshift <- sl3::sl3_Task$new(
      data = av_downshifted,
      covariates = covars,
      outcome = mediator,
      outcome_type = "continuous"
    )

    sl_fit <- sl$train(at_task_noshift)

    # fit new Super Learner to the natural (no shift) data and predict
    at_q_pred <- sl_fit$predict(at_task_noshift)
    av_q_pred <- sl_fit$predict(av_task_noshift)

    # predict with Super Learner from unshifted data on the shifted data
    at_pred_upshifted <- sl_fit$predict(at_task_upshift)
    av_pred_upshifted <- sl_fit$predict(av_task_upshift)

    # predict with Super Learner from unshifted data on the shifted data
    at_pred_upupshifted <- sl_fit$predict(at_task_upupshift)
    av_pred_upupshifted <- sl_fit$predict(av_task_upupshift)

    # predict with Super Learner from unshifted data on the shifted data
    at_pred_downshifted <- sl_fit$predict(at_task_downshift)
    av_pred_downshifted <- sl_fit$predict(av_task_downshift)


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

    med_outcome_list[[i]] <- out_list
  }

  return(med_outcome_list)
}
