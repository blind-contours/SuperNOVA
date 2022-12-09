#' Estimate the Exposure Mechanism via Generalized Propensity Score for Two
#' Exposure Variable
#'
#' @details Compute the propensity score (exposure mechanism) for the observed
#'  data, including the shift. This gives the propensity score for the observed
#'  data (at the observed A) the counterfactual shifted exposure levels (at
#'  {A - delta}, {A + delta}, and {A + 2 * delta}).
#'
#' @param exposure A \code{vector} of characters labeling the exposure
#' variables.
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
#'
#' @export
#'
#' @return A \code{data.table} with four columns, containing estimates of the
#'  generalized propensity score at a downshift (g(A - delta | W)), no shift
#'  (g(A | W)), an upshift (g(A + delta) | W), and an upshift of magnitude two
#'  (g(A + 2 delta) | W).

joint_stoch_shift_est_g_exp <- function(exposures,
                                        delta,
                                        pi_learner,
                                        covars,
                                        av,
                                        at) {
  future::plan(future::sequential, gc = TRUE)

  results <- list()

  for (i in 1:length(exposures)) {
    if (i == 3) {
      exposure <- exposures[[i]][2]
      covariates <- c(covars, exposures[[1]])
    } else {
      exposure <- exposures[[i]]
      covariates <- covars
    }

    # need a data set with the exposure stochastically shifted DOWNWARDS A-delta
    av_downshifted <- data.table::copy(av)

    data.table::set(av_downshifted, j = exposure, value = shift_additive(
      a = subset(av, select = exposure), delta = -delta
    ))

    # need a data set with the exposure stochastically shifted UPWARDS A+delta
    av_upshifted <- data.table::copy(av)
    data.table::set(av_upshifted, j = exposure, value = shift_additive(
      a = subset(av, select = exposure), delta = delta
    ))

    # need a data set with the exposure stochastically shifted UPWARDS A+2delta
    av_upupshifted <- data.table::copy(av)
    data.table::set(av_upupshifted, j = exposure, value = shift_additive(
      a = subset(av, select = exposure), delta = 2 * delta
    ))

    sl_task <- sl3::sl3_Task$new(
      data = at,
      outcome = exposure,
      covariates = covariates
    )

    sl_task_noshift <- sl3::sl3_Task$new(
      data = av,
      outcome = exposure,
      covariates = covariates
    )

    # sl3 task for data with exposure shifted DOWNWARDS A-delta
    sl_task_downshifted <- sl3::sl3_Task$new(
      data = av_downshifted,
      outcome = exposure,
      covariates = covariates
    )

    # sl3 task for data with exposure shifted UPWARDS A+delta
    sl_task_upshifted <- sl3::sl3_Task$new(
      data = av_upshifted,
      outcome = exposure,
      covariates = covariates
    )

    # sl3 task for data with exposure shifted UPWARDS A+2delta
    sl_task_upupshifted <- sl3::sl3_Task$new(
      data = av_upupshifted,
      outcome = exposure,
      covariates = covariates
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

    results[[i]] <- out
  }


  return(results)
}
