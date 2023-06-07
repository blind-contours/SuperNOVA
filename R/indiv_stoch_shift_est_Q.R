#' Estimate the Outcome Mechanism
#'
#' @description
#' Compute the outcome regression for the observed data, taking into account
#' the shift imposed by the intervention. This function returns the outcome
#' regression for the observed data (at A) and under the counterfactual shift
#' (at A + delta).
#'
#' @param exposure A character vector representing the exposures to be shifted.
#' @param delta A numeric value indicating the magnitude of the shift to be
#'   computed for the exposure 'A'. This is passed to the internal
#'   'shift_additive' function and is currently limited to additive shifts.
#' @param mu_learner An object containing a set of instantiated learners from the
#'   'sl3' package, to be used in fitting an ensemble model.
#' @param covars A character vector of covariates to adjust for.
#' @param av A dataframe containing the validation data specific to the fold.
#' @param at A dataframe containing the training data specific to the fold.
#' @param lower_bound A numeric value representing the lower bound for the shifted exposure (optional).
#' @param upper_bound A numeric value representing the upper bound for the shifted exposure (optional).
#'
#' @importFrom stats glm as.formula predict
#' @importFrom data.table as.data.table setnames copy set
#' @importFrom stringr str_detect
#' @importFrom assertthat assert_that
#' @import sl3
#' @export
#'
#' @return A data.table with two columns, containing estimates of the
#'   outcome mechanism at the natural value of the exposure Q(A, W) and an
#'   upshift of the exposure Q(A + delta, W).
#'
#' @examples
#' \dontrun{
#' # Load required libraries
#' library(SuperNOVA)
#' library(sl3)
#'
#' # Create example data
#' set.seed(123)
#' n <- 100
#' W <- rnorm(n)
#' A <- rnorm(n, mean = W)
#' Y <- rnorm(n, mean = A + W)
#' data <- data.frame(W = W, A = A, Y = Y)
#'
#' # Split the data into training and validation sets
#' train_idx <- sample(seq_len(n), size = floor(0.8 * n), replace = FALSE)
#' at <- data[train_idx, ]
#' av <- data[-train_idx, ]
#'
#' # Define learners from the sl3 package
#' mu_learner <- Lrnr_sl("Lrnr_glm_fast")
#' covars <- c("W")
#'
#' # Run the indiv_stoch_shift_est_Q function
#' results <- indiv_stoch_shift_est_Q(
#'   exposure = "A",
#'   delta = 1,
#'   mu_learner = mu_learner,
#'   covars = covars,
#'   av = av,
#'   at = at
#' )
#' }
#'
indiv_stoch_shift_est_Q <- function(exposure,
                                    delta,
                                    mu_learner,
                                    covars,
                                    av,
                                    at,
                                    lower_bound = lower_bound,
                                    upper_bound = upper_bound,
                                    outcome_type = "continuous") {
  future::plan(future::sequential, gc = TRUE)

  # scale the outcome for logit transform
  if (outcome_type != "binary") {
    y_star_av <- scale_to_unit(vals = av$y)
    y_star_at <- scale_to_unit(vals = at$y)

    av$y <- y_star_av
    at$y <- y_star_at
  }


  # need a data set with the exposure stochastically shifted DOWNWARDS A-delta
  # do this for both AV and AT as AT is used in mediation

  av_downshifted <- data.table::copy(av)

  data.table::set(av_downshifted,
    j = exposure,
    value = shift_additive(
      a = av[[exposure]],
      delta = -delta,
      lower_bound = lower_bound,
      upper_bound = upper_bound
    )
  )

  at_downshifted <- data.table::copy(at)
  data.table::set(at_downshifted,
    j = exposure,
    value = shift_additive(
      a = at[[exposure]],
      delta = -delta,
      lower_bound = lower_bound,
      upper_bound = upper_bound
    )
  )

  # need a data set with the exposure stochastically shifted UPWARDS A+delta
  # do this for both AV and AT as AT is used in mediation

  av_upshifted <- data.table::copy(av)
  data.table::set(av_upshifted,
    j = exposure,
    value = shift_additive(
      a = av[[exposure]],
      delta = delta,
      lower_bound = lower_bound,
      upper_bound = upper_bound
    )
  )

  at_upshifted <- data.table::copy(at)
  data.table::set(at_upshifted,
    j = exposure,
    value = shift_additive(
      a = at[[exposure]],
      delta = delta,
      lower_bound = lower_bound,
      upper_bound = upper_bound
    )
  )

  # need a data set with the exposure stochastically shifted UPWARDS A+2delta
  # do this for both AV and AT as AT is used in mediation

  av_upupshifted <- data.table::copy(av)
  data.table::set(av_upupshifted,
    j = exposure,
    value = shift_additive(
      a = av[[exposure]],
      delta = 2 * delta,
      lower_bound = lower_bound,
      upper_bound = upper_bound
    )
  )

  at_upupshifted <- data.table::copy(at)
  data.table::set(at_upupshifted,
    j = exposure,
    value = shift_additive(
      a = at[[exposure]],
      delta = 2 * delta,
      lower_bound = lower_bound,
      upper_bound = upper_bound
    )
  )

  # Outcome mechanism
  sl <- sl3::Lrnr_sl$new(
    learners = mu_learner,
    metalearner = sl3::Lrnr_nnls$new()
  )

  at_task_noshift <- suppressMessages(sl3::sl3_Task$new(
    data = at,
    covariates = covars,
    outcome = "y",
    outcome_type = "quasibinomial"
  ))

  av_task_noshift <- suppressMessages(sl3::sl3_Task$new(
    data = av,
    covariates = covars,
    outcome = "y",
    outcome_type = "quasibinomial"
  ))

  at_task_upshift <- suppressMessages(sl3::sl3_Task$new(
    data = at_upshifted,
    covariates = covars,
    outcome = "y",
    outcome_type = "quasibinomial"
  ))

  av_task_upshift <- suppressMessages(sl3::sl3_Task$new(
    data = av_upshifted,
    covariates = covars,
    outcome = "y",
    outcome_type = "quasibinomial"
  ))

  at_task_upupshift <- suppressMessages(sl3::sl3_Task$new(
    data = at_upupshifted,
    covariates = covars,
    outcome = "y",
    outcome_type = "quasibinomial"
  ))

  av_task_upupshift <- suppressMessages(sl3::sl3_Task$new(
    data = av_upupshifted,
    covariates = covars,
    outcome = "y",
    outcome_type = "quasibinomial"
  ))

  at_task_downshift <- suppressMessages(sl3::sl3_Task$new(
    data = at_downshifted,
    covariates = covars,
    outcome = "y",
    outcome_type = "quasibinomial"
  ))

  av_task_downshift <- suppressMessages(sl3::sl3_Task$new(
    data = av_downshifted,
    covariates = covars,
    outcome = "y",
    outcome_type = "quasibinomial"
  ))

  sl_fit <- suppressWarnings(suppressMessages(sl$train(at_task_noshift)))

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
