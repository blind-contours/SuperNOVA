#' Compute the Pooled Shift Parameter Estimate From the Fold Specific Results
#'
#' @details Estimate the value of the pooledcausal parameter alongside statistical
#'  inference for the parameter estimate based on the nuisance parameters from
#'  the fold specific results.
#'
#' @param indiv_shift_results List of individual shift results from the
#' parallelized CV estimation
#' @param estimator The type of estimator to be fit, either \code{"tmle"} for
#'  targeted maximum likelihood estimation or \code{"onestep"} for a one-step
#'  estimator.
#' @param fluc_mod_out An object giving values of the logistic tilting model
#'  for targeted minimum loss estimation. This type of object should be the
#'  output of the internal routines to perform this step of the TML estimation
#'  procedure, as given by \code{\link{fit_fluctuation}}.
#' @param fluctuation Type of fluctuation to use
#'
#' @importFrom stats var
#'
#' @return A \code{list} containing the parameter estimate, estimated variance
#'  based on the efficient influence function (EIF), the estimate of the EIF
#'  incorporating inverse probability of censoring weights, and the estimate of
#'  the EIF without the application of such weights.
calc_pooled_indiv_shifts <- function(indiv_shift_results,
                                     estimator = c("tmle", "onestep"),
                                     fluc_mod_out = NULL,
                                     fluctuation) {
  # set TMLE as default estimator type
  estimator <- match.arg(estimator)

  results_list <- list()

  names <- names(indiv_shift_results)
  names <- gsub("^.+?:", "", names)
  names <- stringr::str_trim(names)


  for (var_set in unique(names)) {
    var_set_results <- indiv_shift_results[stringr::str_detect(
      names(indiv_shift_results), var_set
    )]

    if (length(var_set_results) != 0) {
      test <- unlist(var_set_results, recursive = FALSE)

      Hn <- do.call(rbind, test[stringr::str_detect(
        names(test), "Hn"
      )])

      Qn_scaled <- do.call(rbind, test[stringr::str_detect(
        names(test), "Qn_scaled"
      )])

      data <- do.call(rbind, test[stringr::str_detect(
        names(test), "data"
      )])

      k_fold_results <- do.call(rbind, test[stringr::str_detect(
        names(test), "k_fold"
      )])

      deltas <- do.call(rbind, test[stringr::str_detect(
        names(test), "Delta"
      )])

      tmle_fit <- tmle_exposhift(
        data_internal = data,
        Qn_scaled = Qn_scaled,
        Hn = Hn,
        fluctuation = fluctuation,
        y = data$y,
        delta = mean(deltas)
      )

      indiv_shift_in_fold <- calc_final_ind_shift_param(
        tmle_fit = tmle_fit,
        exposure = var_set,
        fold_k = "Pooled TMLE"
      )

      indiv_shift_in_fold$Delta <- mean(deltas)

      results_df <- rbind(k_fold_results, indiv_shift_in_fold)

      results_list[[var_set]] <- results_df
    } else {

    }
  }


  return(results_list)
}
