#' Estimate the Exposure Mechanism via Generalized Propensity Score for One
#' Mediator Variable and Shift an Exposure
#'
#' @details Compute the propensity score (mediator mechanism) for the observed
#'  data, including the shift. This gives the propensity score for the observed
#'  data for the mediator (at the observed A) the counterfactual shifted exposure levels (at
#'  {A - delta}, {A + delta}, and {A + 2 * delta}).
#'
#' @param exposure A \code{character} labeling the exposure variable.
#' @param mediator Z \code{character} labeling the mediator variable.
#' @param covars A \code{character} labeling the covariate variables
#' @param delta A \code{numeric} value identifying a shift in the observed
#'  value of the exposure under which observations are to be evaluated.
#' @param pi_learner Object containing a set of instantiated learners
#'  from \pkg{sl3}, to be used in fitting an ensemble model.
#' @param av A \code{dataframe} of validation data specific to the fold
#' @param at A \code{dataframe} of training data specific to the fold
#'
#' @importFrom data.table as.data.table setnames set copy
#' @importFrom stats predict
#' @importFrom haldensify haldensify
#' @importFrom assertthat assert_that
#' @import sl3
#'
#' @export
#'
#' @return A \code{data.table} with four columns, containing estimates of the
#'  generalized propensity score at a downshift (g(A - delta | W)), no shift
#'  (g(A | W)), an upshift (g(A + delta) | W), and an upshift of magnitude two
#'  (g(A + 2 delta) | W).

joint_stoch_shift_est_z_exp <- function(exposures,
                                        mediator,
                                        deltas,
                                        pi_learner,
                                        w_names,
                                        a_names,
                                        z_names,
                                        av,
                                        at) {
  future::plan(future::sequential, gc = TRUE)

  results <- list()
  delta_results <- list()
  Hn_result <- list()

  av <- as.data.frame(av)
  at <- as.data.frame(at)

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

    sl_task <- sl3::sl3_Task$new(
      data = at,
      outcome = mediator,
      covariates = c(w_names)
    )

    sl_task_noshift <- sl3::sl3_Task$new(
      data = av,
      outcome = mediator,
      covariates = c(w_names)
    )

    # sl3 task for data with exposure shifted DOWNWARDS A-delta
    sl_task_downshifted <- sl3::sl3_Task$new(
      data = av_downshifted,
      outcome = mediator,
      covariates = c(w_names)
    )

    # sl3 task for data with exposure shifted UPWARDS A+delta
    sl_task_upshifted <- sl3::sl3_Task$new(
      data = av_upshifted,
      outcome = mediator,
      covariates = c(w_names)
    )

    # sl3 task for data with exposure shifted UPWARDS A+2delta
    sl_task_upupshifted <- sl3::sl3_Task$new(
      data = av_upupshifted,
      outcome = mediator,
      covariates = c(w_names)
    )

    sl_fit <- pi_learner$train(sl_task)

    pred_g_exp_noshift <- sl_fit$predict(sl_task_noshift)
    pred_g_exp_downshifted <- sl_fit$predict(sl_task_downshifted)
    pred_g_exp_upshifted <- sl_fit$predict(sl_task_upshifted)
    pred_g_exp_upupshifted <- sl_fit$predict(sl_task_upupshifted)

    # create output data.tables
    out <- data.table::as.data.table(cbind.data.frame(
      pred_g_exp_downshifted,
      pred_g_exp_noshift,
      pred_g_exp_upshifted,
      pred_g_exp_upupshifted
    ))

    data.table::setnames(out, c("downshift", "noshift", "upshift", "upupshift"))

    Hn <- est_hn(gn_exp = out)

    results[[i]] <- out
    delta_results[[i]] <- delta
    Hn_result[[i]] <- Hn
  }


  return(list(
    "gn_results" = results,
    "delta_results" = delta_results,
    "Hn_results" = Hn_result,
    "model" = sl_fit
  ))
}
