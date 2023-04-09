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
#'
#' @export
#'
#' @return A \code{data.table} with four columns, containing estimates of the
#'  generalized propensity score at a downshift (g(A - delta | W)), no shift
#'  (g(A | W)), an upshift (g(A + delta) | W), and an upshift of magnitude two
#'  (g(A + 2 delta) | W).

indiv_stoch_shift_e_estimator <- function(exposure,
                                          mediator,
                                          delta,
                                          pi_learner,
                                          w_names = w_names,
                                          a_names = a_names,
                                          z_names = z_names,
                                          av,
                                          at) {
  future::plan(future::sequential, gc = TRUE)

  covars <- c(w_names, z_names)

  # need a data set with the exposure stochastically shifted DOWNWARDS A-delta
  av_downshifted <- data.table::copy(av)

  data.table::set(av_downshifted,
    j = mediator,
    value = shift_additive(
      a = subset(av, select = mediator), delta = -delta
    )
  )

  at_downshifted <- data.table::copy(at)

  data.table::set(at_downshifted,
    j = mediator,
    value = shift_additive(
      a = subset(at, select = mediator), delta = -delta
    )
  )

  # need a data set with the exposure stochastically shifted UPWARDS A+delta
  av_upshifted <- data.table::copy(av)
  data.table::set(av_upshifted,
    j = mediator,
    value = shift_additive(
      a = subset(av, select = mediator), delta = delta
    )
  )

  at_upshifted <- data.table::copy(at)
  data.table::set(at_upshifted,
    j = mediator,
    value = shift_additive(
      a = subset(at, select = mediator), delta = delta
    )
  )

  # need a data set with the exposure stochastically shifted UPWARDS A+2delta
  av_upupshifted <- data.table::copy(av)
  data.table::set(av_upupshifted,
    j = mediator,
    value = shift_additive(
      a = subset(av, select = mediator), delta = 2 * delta
    )
  )

  at_upupshifted <- data.table::copy(at)
  data.table::set(at_upupshifted,
    j = mediator,
    value = shift_additive(
      a = subset(at, select = mediator), delta = 2 * delta
    )
  )

  sl_task <- sl3::sl3_Task$new(
    data = at,
    outcome = exposure,
    covariates = covars
  )

  sl_task_noshift_at <- sl3::sl3_Task$new(
    data = at,
    outcome = exposure,
    covariates = covars
  )

  sl_task_noshift_av <- sl3::sl3_Task$new(
    data = av,
    outcome = exposure,
    covariates = covars
  )

  # sl3 task for data with exposure shifted DOWNWARDS A-delta
  sl_task_downshifted_at <- sl3::sl3_Task$new(
    data = at_downshifted,
    outcome = exposure,
    covariates = covars
  )

  sl_task_downshifted_av <- sl3::sl3_Task$new(
    data = av_downshifted,
    outcome = exposure,
    covariates = covars
  )

  # sl3 task for data with exposure shifted UPWARDS A+delta
  sl_task_upshifted_at <- sl3::sl3_Task$new(
    data = at_upshifted,
    outcome = exposure,
    covariates = covars
  )

  sl_task_upshifted_av <- sl3::sl3_Task$new(
    data = av_upshifted,
    outcome = exposure,
    covariates = covars
  )

  # sl3 task for data with exposure shifted UPWARDS A+2delta
  sl_task_upupshifted_at <- sl3::sl3_Task$new(
    data = at_upupshifted,
    outcome = exposure,
    covariates = covars
  )

  sl_task_upupshifted_av <- sl3::sl3_Task$new(
    data = av_upupshifted,
    outcome = exposure,
    covariates = covars
  )

  sl_fit <- pi_learner$train(sl_task)

  # at predictions -----------
  pred_z_exp_noshift_at <- sl_fit$predict(sl_task_noshift_at)
  pred_z_exp_downshifted_at <- sl_fit$predict(sl_task_downshifted_at)
  pred_z_exp_upshifted_at <- sl_fit$predict(sl_task_upshifted_at)
  pred_z_exp_upupshifted_at <- sl_fit$predict(sl_task_upupshifted_at)

  # av predictions -----------

  pred_z_exp_noshift_av <- sl_fit$predict(sl_task_noshift_av)
  pred_z_exp_downshifted_av <- sl_fit$predict(sl_task_downshifted_av)
  pred_z_exp_upshifted_av <- sl_fit$predict(sl_task_upshifted_av)
  pred_z_exp_upupshifted_av <- sl_fit$predict(sl_task_upupshifted_av)


  # create output data.tables
  av_out <- cbind.data.frame(
    pred_z_exp_downshifted_av,
    pred_z_exp_noshift_av,
    pred_z_exp_upshifted_av,
    pred_z_exp_upupshifted_av
  )

  at_out <- cbind.data.frame(
    pred_z_exp_downshifted_at,
    pred_z_exp_noshift_at,
    pred_z_exp_upshifted_at,
    pred_z_exp_upupshifted_at
  )

  data.table::setnames(av_out, c("downshift", "noshift", "upshift", "upupshift"))
  data.table::setnames(at_out, c("downshift", "noshift", "upshift", "upupshift"))

  return(list("av" = av_out, "at" = at_out))
}
