#' Estimate the Exposure Mechanism via Generalized Propensity Score for Two
#' Exposure Variable
#'
#' @details Compute the propensity score (exposure mechanism) for the observed
#'  data, including the shift. This gives the propensity score for the observed
#'  data (at the observed A) the counterfactual shifted exposure levels (at
#'  {A - delta}, {A + delta}, and {A + 2 * delta}).
#'
#' @param exposures A \code{vector} of characters labeling the exposure
#' variables.
#' @param covars A \code{character} labeling the covariate variables
#' @param deltas A \code{numeric} value identifying a shift in the observed
#'  value of the exposure under which observations are to be evaluated.
#' @param pi_learner Object containing a set of instantiated learners
#'  from \pkg{sl3}, to be used in fitting an ensemble model.
#' @param av A \code{dataframe} of validation data specific to the fold
#' @param at A \code{dataframe} of training data specific to the fold
#' @param adaptive_delta Whether to adaptively change the delta based on positivity
#' determined from the clever covariate being below the hn_trunc_thresh level
#' @param hn_trunc_thresh Truncation level of the clever covariate used in the
#' adaptive delta method
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
                                        deltas,
                                        g_learner,
                                        covars,
                                        av,
                                        at,
                                        adaptive_delta,
                                        hn_trunc_thresh,
                                        lower_bound,
                                        upper_bound) {
  future::plan(future::sequential, gc = TRUE)

  # Function to create shifted data
  create_shifted_data <- function(data, exposure, delta, lower_bound, upper_bound) {
    shifted_data <- data.table::copy(data)
    data.table::set(shifted_data,
                    j = exposure,
                    value = shift_additive(
                      a = subset(data, select = exposure),
                      delta = delta,
                      lower_bound = lower_bound,
                      upper_bound = upper_bound
                    )
    )
    return(shifted_data)
  }

  results <- list()
  delta_results <- list()
  Hn_result <- list()

  for (i in 1:length(exposures)) {
    if (i == 3) {
      exposure <- exposures[[i]][2]
      covariates <- c(covars, exposures[[1]])
      delta <- delta_results[[2]]
      adaptive_delta <- FALSE
    } else {
      exposure <- exposures[[i]]
      delta <- deltas[[i]]
      covariates <- covars
    }

    # Creating shifted data
    av_downshifted <- create_shifted_data(av, exposure, -delta, lower_bound, upper_bound)
    at_downshifted <- create_shifted_data(at, exposure, -delta, lower_bound, upper_bound)

    av_upshifted <- create_shifted_data(av, exposure, delta, lower_bound, upper_bound)
    at_upshifted <- create_shifted_data(at, exposure, delta, lower_bound, upper_bound)

    av_upupshifted <- create_shifted_data(av, exposure, 2 * delta, lower_bound, upper_bound)
    at_upupshifted <- create_shifted_data(at, exposure, 2 * delta, lower_bound, upper_bound)


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

    delta_reduced <- delta

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

          av_downshifted <- create_shifted_data(av, exposure, -delta_reduced, lower_bound, upper_bound)
          at_downshifted <- create_shifted_data(at, exposure, -delta_reduced, lower_bound, upper_bound)
          av_upshifted <- create_shifted_data(av, exposure, delta_reduced, lower_bound, upper_bound)
          at_upshifted <- create_shifted_data(at, exposure, delta_reduced, lower_bound, upper_bound)
          av_upupshifted <- create_shifted_data(av, exposure, 2 * delta_reduced, lower_bound, upper_bound)
          at_upupshifted <- create_shifted_data(at, exposure, 2 * delta_reduced, lower_bound, upper_bound)

          if (density_type == "sl") {

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

          max_hn <- max(Hn_at$shift)
        } else {
          pos_violation <- FALSE
        }
      }
    }


    results[[i]] <- av_out
    delta_results[[i]] <- delta_reduced
    Hn_result[[i]] <- Hn_av
  }


  return(list(
    "gn_results" = results,
    "delta_results" = delta_results,
    "Hn_results" = Hn_result
  ))
}
