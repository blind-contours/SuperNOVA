#' Compute the Pooled Shift Parameter Estimate From the Fold Specific Results
#' For the Effect Modification Parameter
#'
#' @details Estimate the value of the pooled causal parameter alongside statistical
#'  inference for the parameter estimate based on the nuisance parameters from
#'  the fold specific results.
#'
#' @param y A \code{numeric} vector of the observed outcomes.
#' @param fluc_mod_out An object giving values of the logistic tilting model
#'  for targeted minimum loss estimation. This type of object should be the
#'  output of the internal routines to perform this step of the TML estimation
#'  procedure, as given by \code{\link{fit_fluctuation}}.
#' @param em_shift_results Shift results for effect modification found across
#' the folds.
#' @param estimator Type of estimator to use (TMLE or one-step)
#' @param fluctuation type of fluctuation to use
#' @param w_names Names of covariates
#' @param a_names Names of exposure
#' @param em_learner Stack of decision trees to regress influence function
#' onto covariate space
#'
#' @importFrom stats var
#'
#' @return A \code{list} containing the parameter estimate, estimated variance
#'  based on the efficient influence function (EIF), the estimate of the EIF
#'  incorporating inverse probability of censoring weights, and the estimate of
#'  the EIF without the application of such weights.
calc_pooled_em_shifts <- function(y,
                                  em_shift_results,
                                  estimator = c("tmle", "onestep"),
                                  fluc_mod_out = NULL,
                                  w_names,
                                  a_names,
                                  fluctuation,
                                  em_learner) {
  # set TMLE as default estimator type
  estimator <- match.arg(estimator)

  names <- names(em_shift_results)
  names <- gsub("^.+?:", "", names)
  names <- stringr::str_trim(names)

  results_list <- list()

  for (var_set in unique(names)) {
    var_set_results <- em_shift_results[stringr::str_detect(
      names(em_shift_results), var_set
    )]

    if (length(var_set_results) != 0) {
      test <- unlist(var_set_results, recursive = FALSE)

      Hn_at <- do.call(rbind, test[stringr::str_detect(
        names(test), "Hn_at"
      )])

      Hn_av <- do.call(rbind, test[stringr::str_detect(
        names(test), "Hn_av"
      )])

      Qn_scaled_at <- do.call(rbind, test[stringr::str_detect(
        names(test), "Qn_scaled_at"
      )])

      Qn_scaled_av <- do.call(rbind, test[stringr::str_detect(
        names(test), "Qn_scaled_av"
      )])

      data <- do.call(rbind, test[stringr::str_detect(
        names(test), "av_data"
      )])

      k_fold_results <- do.call(rbind, test[stringr::str_detect(
        names(test), "k_fold"
      )])

      deltas <- do.call(rbind, test[stringr::str_detect(
        names(test), "delta"
      )])

      tmle_fit <- tmle_exposhift(
        data_internal = data,
        delta = mean(deltas),
        Qn_scaled = Qn_scaled_av,
        Hn = Hn_av,
        fluctuation = fluctuation,
        y = data$y
      )

      effect_m_name <- stringr::str_extract(var_set, paste(w_names, collapse = "|"))
      exposure <- stringr::str_extract(var_set, paste(a_names, collapse = "|"))

      effect_mod_in_fold <- calc_final_effect_mod_param(
        tmle_fit_av = tmle_fit,
        tmle_fit_at = tmle_fit,
        exposure = exposure,
        at = data,
        av = data,
        effect_m_name = effect_m_name,
        fold_k = "Pooled TMLE",
        em_learner = em_learner
      )

      effect_mod_in_fold$Delta <- mean(deltas)

      results_df <- rbind(k_fold_results, effect_mod_in_fold)

      results_list[[var_set]] <- results_df
    } else {

    }
  }


  return(results_list)
}
