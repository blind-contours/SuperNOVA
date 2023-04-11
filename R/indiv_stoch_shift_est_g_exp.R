#' Estimate the Exposure Mechanism via Generalized Propensity Score for One
#' Exposure Variable
#'
#' @details Compute the propensity score (exposure mechanism) for the observed
#'  data, including the shift. This gives the propensity score for the observed
#'  data (at the observed A) the counterfactual shifted exposure levels (at
#'  {A - delta}, {A + delta}, and {A + 2 * delta}).
#'
#' @param exposure A \code{character} labeling the exposure variable.
#' @param covars A \code{character} labeling the covariate variables
#' @param delta A \code{numeric} value identifying a shift in the observed
#'  value of the exposure under which observations are to be evaluated.
#' @param pi_learner Object containing a set of instantiated learners
#'  from \pkg{sl3}, to be used in fitting an ensemble model.
#' @param av A \code{dataframe} of validation data specific to the fold
#' @param at A \code{dataframe} of training data specific to the fold
#' @param adaptive_delta Whether or not to adaptively adjust delta based on
#' positivity (estimated from the clever covariate) meeting the hn_trunc_thresh
#' level
#' @param hn_trunc_thresh Level of the clever covariate in the adaptive delta
#' procedure
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

indiv_stoch_shift_est_g_exp <- function(exposure,
                                        delta,
                                        g_learner,
                                        covars,
                                        av,
                                        at,
                                        adaptive_delta,
                                        hn_trunc_thresh,
                                        exposure_quantized,
                                        lower_bound,
                                        upper_bound,
                                        outcome_type,
                                        density_type,
                                        n_bins,
                                        max_degree) {
  future::plan(future::sequential, gc = TRUE)

  # need a data set with the exposure stochastically shifted DOWNWARDS A-delta
  av_downshifted <- data.table::copy(av)

  data.table::set(av_downshifted,
    j = exposure,
    value = shift_additive(
      a = subset(av, select = exposure),
      delta = -delta,
      lower_bound = lower_bound,
      upper_bound = upper_bound
    )
  )

  at_downshifted <- data.table::copy(at)

  data.table::set(at_downshifted,
    j = exposure,
    value = shift_additive(
      a = subset(at, select = exposure),
      delta = -delta,
      lower_bound = lower_bound,
      upper_bound = upper_bound
    )
  )

  # need a data set with the exposure stochastically shifted UPWARDS A+delta
  av_upshifted <- data.table::copy(av)
  data.table::set(av_upshifted,
    j = exposure,
    value = shift_additive(
      a = subset(av, select = exposure),
      delta = delta,
      lower_bound = lower_bound,
      upper_bound = upper_bound
    )
  )

  at_upshifted <- data.table::copy(at)
  data.table::set(at_upshifted,
    j = exposure,
    value = shift_additive(
      a = subset(at, select = exposure),
      delta = delta,
      lower_bound = lower_bound,
      upper_bound = upper_bound
    )
  )

  # need a data set with the exposure stochastically shifted UPWARDS A+2delta
  av_upupshifted <- data.table::copy(av)
  data.table::set(av_upupshifted,
    j = exposure,
    value = shift_additive(
      a = subset(av, select = exposure),
      delta = 2 * delta,
      lower_bound = lower_bound,
      upper_bound = upper_bound
    )
  )

  at_upupshifted <- data.table::copy(at)
  data.table::set(at_upupshifted,
    j = exposure,
    value = shift_additive(
      a = subset(at, select = exposure),
      delta = 2 * delta,
      lower_bound = lower_bound,
      upper_bound = upper_bound
    )
  )

  if (density_type == "sl") {
    sl_task <- sl3::sl3_Task$new(
      data = at,
      outcome = exposure,
      covariates = covars,
      outcome_type = outcome_type
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

    g_model <- suppressWarnings(suppressMessages(g_learner$train(sl_task)))

    # at predictions -----------
    pred_g_exp_noshift_at <- g_model$predict(sl_task_noshift_at)
    pred_g_exp_downshifted_at <- g_model$predict(sl_task_downshifted_at)
    pred_g_exp_upshifted_at <- g_model$predict(sl_task_upshifted_at)
    pred_g_exp_upupshifted_at <- g_model$predict(sl_task_upupshifted_at)

    # av predictions -----------

    pred_g_exp_noshift_av <- g_model$predict(sl_task_noshift_av)
    pred_g_exp_downshifted_av <- g_model$predict(sl_task_downshifted_av)
    pred_g_exp_upshifted_av <- g_model$predict(sl_task_upshifted_av)
    pred_g_exp_upupshifted_av <- g_model$predict(sl_task_upupshifted_av)

    if (exposure_quantized == TRUE) {
      get_value_from_column <- function(a, row_quantile_predictions, lower_bound = lower_bound, upper_bound = upper_bound) {
        a <- ifelse(a < lower_bound, lower_bound, a)
        a <- ifelse(a > upper_bound, upper_bound, a)
        row_quantile_predictions[[a]]
      }
      pred_g_exp_noshift_at <- mapply(get_value_from_column, a = at$a, row_quantile_predictions = unlist(pred_g_exp_noshift_at, recursive = FALSE), lower_bound = lower_bound, upper_bound = upper_bound)
      pred_g_exp_downshifted_at <- mapply(get_value_from_column, a = at$a - delta, row_quantile_predictions = unlist(pred_g_exp_downshifted_at, recursive = FALSE), lower_bound = lower_bound, upper_bound = upper_bound)
      pred_g_exp_upshifted_at <- mapply(get_value_from_column, a = at$a + delta, row_quantile_predictions = unlist(pred_g_exp_upshifted_at, recursive = FALSE), lower_bound = lower_bound, upper_bound = upper_bound)
      pred_g_exp_upupshifted_at <- mapply(get_value_from_column, a = at$a + 2 * delta, row_quantile_predictions = unlist(pred_g_exp_upupshifted_at, recursive = FALSE), lower_bound = lower_bound, upper_bound = upper_bound)

      pred_g_exp_noshift_av <- mapply(get_value_from_column, a = av$a, row_quantile_predictions = unlist(pred_g_exp_noshift_av, recursive = FALSE), lower_bound = lower_bound, upper_bound = upper_bound)
      pred_g_exp_downshifted_av <- mapply(get_value_from_column, a = av$a - delta, row_quantile_predictions = unlist(pred_g_exp_downshifted_av, recursive = FALSE), lower_bound = lower_bound, upper_bound = upper_bound)
      pred_g_exp_upshifted_av <- mapply(get_value_from_column, a = av$a + delta, row_quantile_predictions = unlist(pred_g_exp_upshifted_av, recursive = FALSE), lower_bound = lower_bound, upper_bound = upper_bound)
      pred_g_exp_upupshifted_av <- mapply(get_value_from_column, a = av$a + 2 * delta, row_quantile_predictions = unlist(pred_g_exp_upupshifted_av, recursive = FALSE), lower_bound = lower_bound, upper_bound = upper_bound)
    }
  } else {
    at <- as.data.frame(at)
    at_upshifted <- as.data.frame(at_upshifted)
    at_upupshifted <- as.data.frame(at_upupshifted)
    at_downshifted <- as.data.frame(at_downshifted)

    av <- as.data.frame(av)
    av_upshifted <- as.data.frame(av_upshifted)
    av_upupshifted <- as.data.frame(av_upupshifted)
    av_downshifted <- as.data.frame(av_downshifted)


    g_model <- suppressMessages(haldensify(
      A = at[[exposure]], W = at[w_names], n_bins = n_bins, lambda_seq = exp(seq(-1, -10, length = 100)),
      # the following arguments are passed to hal9001::fit_hal()
      max_degree = max_degree
    ))

    # at predictions -----------
    pred_g_exp_noshift_at <- predict(g_model, new_A = at[[exposure]], new_W = at[w_names])
    pred_g_exp_downshifted_at <- predict(g_model, new_A = at_downshifted[[exposure]], new_W = at_downshifted[w_names])
    pred_g_exp_upshifted_at <- predict(g_model, new_A = at_upshifted[[exposure]], new_W = at_upshifted[w_names])
    pred_g_exp_upupshifted_at <- predict(g_model, new_A = at_upupshifted[[exposure]], new_W = at_upupshifted[w_names])

    # av predictions -----------

    pred_g_exp_noshift_av <- predict(g_model, new_A = av[[exposure]], new_W = av[w_names])
    pred_g_exp_downshifted_av <- predict(g_model, new_A = av_downshifted[[exposure]], new_W = av_downshifted[w_names])
    pred_g_exp_upshifted_av <- predict(g_model, new_A = av_upshifted[[exposure]], new_W = av_upshifted[w_names])
    pred_g_exp_upupshifted_av <- predict(g_model, new_A = av_upupshifted[[exposure]], new_W = av_upupshifted[w_names])
  }

  # create output data.tables
  av_out <- cbind.data.frame(
    pred_g_exp_downshifted_av,
    pred_g_exp_noshift_av,
    pred_g_exp_upshifted_av,
    pred_g_exp_upupshifted_av
  )

  at_out <- cbind.data.frame(
    pred_g_exp_downshifted_at,
    pred_g_exp_noshift_at,
    pred_g_exp_upshifted_at,
    pred_g_exp_upupshifted_at
  )

  data.table::setnames(av_out, c("downshift", "noshift", "upshift", "upupshift"))
  data.table::setnames(at_out, c("downshift", "noshift", "upshift", "upupshift"))

  Hn_av <- est_hn(gn_exp = av_out)
  Hn_at <- est_hn(gn_exp = at_out)

  delta_reduced <- delta

  if (adaptive_delta == TRUE) {
    max_hn <- max(Hn_at$shift)
    if (max_hn > hn_trunc_thresh) {
      pos_violation <- TRUE
    } else {
      pos_violation <- FALSE
    }

    while (pos_violation == TRUE) {
      if (max_hn > hn_trunc_thresh) {
        delta_reduced <- delta_reduced - (0.1 * delta_reduced)
        # need a data set with the exposure stochastically shifted DOWNWARDS A-delta
        av_downshifted <- data.table::copy(av)

        data.table::set(av_downshifted,
          j = exposure,
          value = shift_additive(
            a = subset(av, select = exposure),
            delta = -delta_reduced,
            lower_bound = lower_bound,
            upper_bound = upper_bound
          )
        )

        # need a data set with the exposure stochastically shifted DOWNWARDS A-delta
        av_downshifted <- data.table::copy(av)

        data.table::set(av_downshifted,
          j = exposure,
          value = shift_additive(
            a = subset(av, select = exposure),
            delta = -delta_reduced,
            lower_bound = lower_bound,
            upper_bound = upper_bound
          )
        )

        at_downshifted <- data.table::copy(at)

        data.table::set(at_downshifted,
          j = exposure,
          value = shift_additive(
            a = subset(at, select = exposure),
            delta = -delta_reduced,
            lower_bound = lower_bound,
            upper_bound = upper_bound
          )
        )

        # need a data set with the exposure stochastically shifted UPWARDS A+delta
        av_upshifted <- data.table::copy(av)
        data.table::set(av_upshifted,
          j = exposure,
          value = shift_additive(
            a = subset(av, select = exposure),
            delta = delta_reduced,
            lower_bound = lower_bound,
            upper_bound = upper_bound
          )
        )

        at_upshifted <- data.table::copy(at)
        data.table::set(at_upshifted,
          j = exposure,
          value = shift_additive(
            a = subset(at, select = exposure),
            delta = delta_reduced,
            lower_bound = lower_bound,
            upper_bound = upper_bound
          )
        )

        # need a data set with the exposure stochastically shifted UPWARDS A+2delta
        av_upupshifted <- data.table::copy(av)
        data.table::set(av_upupshifted,
          j = exposure,
          value = shift_additive(
            a = subset(av, select = exposure),
            delta = 2 * delta_reduced,
            lower_bound = lower_bound,
            upper_bound = upper_bound
          )
        )

        at_upupshifted <- data.table::copy(at)
        data.table::set(at_upupshifted,
          j = exposure,
          value = shift_additive(
            a = subset(at, select = exposure),
            delta = 2 * delta_reduced,
            lower_bound = lower_bound,
            upper_bound = upper_bound
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

        # at predictions -----------
        pred_g_exp_noshift_at <- sl_fit$predict(sl_task_noshift_at)
        pred_g_exp_downshifted_at <- sl_fit$predict(sl_task_downshifted_at)
        pred_g_exp_upshifted_at <- sl_fit$predict(sl_task_upshifted_at)
        pred_g_exp_upupshifted_at <- sl_fit$predict(sl_task_upupshifted_at)

        # av predictions -----------

        pred_g_exp_noshift_av <- sl_fit$predict(sl_task_noshift_av)
        pred_g_exp_downshifted_av <- sl_fit$predict(sl_task_downshifted_av)
        pred_g_exp_upshifted_av <- sl_fit$predict(sl_task_upshifted_av)
        pred_g_exp_upupshifted_av <- sl_fit$predict(sl_task_upupshifted_av)

        # create output data.tables
        av_out <- cbind.data.frame(
          pred_g_exp_downshifted_av,
          pred_g_exp_noshift_av,
          pred_g_exp_upshifted_av,
          pred_g_exp_upupshifted_av
        )

        # create output data.tables
        at_out <- cbind.data.frame(
          pred_g_exp_downshifted_at,
          pred_g_exp_noshift_at,
          pred_g_exp_upshifted_at,
          pred_g_exp_upupshifted_at
        )

        data.table::setnames(av_out, c("downshift", "noshift", "upshift", "upupshift"))
        data.table::setnames(at_out, c("downshift", "noshift", "upshift", "upupshift"))

        Hn_at <- est_hn(gn_exp = at_out)
        max_hn <- max(Hn_at$shift)
        Hn_av <- est_hn(gn_exp = av_out)
      } else {
        pos_violation <- FALSE
      }
    }
  }

  return(list(
    "av" = av_out,
    "at" = at_out,
    "delta" = delta_reduced,
    "Hn_at" = Hn_at,
    "Hn_av" = Hn_av,
    "model" = g_model
  ))
}
