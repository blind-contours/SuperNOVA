#' Compute the Pooled Joint Mediation Shift Parameter Estimate From the Fold
#' Specific Results
#'
#' @details Estimate the value of the pooled causal parameter alongside statistical
#'  inference for the parameter estimate based on the nuisance parameters from
#'  the fold specific results for the mediation parameter.
#'
#' @param med_shift_results List of data frames of the mediation results across
#' the parallelized CV procedure
#' @param estimator The type of estimator to be fit, either \code{"tmle"} for
#'  targeted maximum likelihood estimation or \code{"onestep"} for a one-step
#'  estimator.
#' @param fluc_mod_out An object giving values of the logistic tilting model
#'  for targeted minimum loss estimation. This type of object should be the
#'  output of the internal routines to perform this step of the TML estimation
#'  procedure, as given by \code{\link{fit_fluctuation}}.
#' @param a_names List of exposure names
#' @param w_names List of covariate names
#' @param z_names List of mediator names
#' @param y_name Outcome name
#' @param fluctuation Type of fluctuation to be used
#' @importFrom stats var
#'
#' @return A \code{list} containing the parameter estimate, estimated variance
#'  based on the efficient influence function (EIF), the estimate of the EIF
#'  incorporating inverse probability of censoring weights, and the estimate of
#'  the EIF without the application of such weights.
calc_pooled_joint_med_shifts <- function(joint_med_shift_results,
                                         estimator = c("tmle", "onestep"),
                                         fluc_mod_out = NULL,
                                         a_names,
                                         w_names,
                                         z_names,
                                         y_name,
                                         fluctuation) {
  # set TMLE as default estimator type
  estimator <- match.arg(estimator)

  names <- names(joint_med_shift_results)
  names <- gsub("^.+?:", "", names)
  names <- stringr::str_trim(names)

  results_list <- list()

  for (var_set in unique(names)) {
    var_set_results <- joint_med_shift_results[stringr::str_detect(
      names(joint_med_shift_results), var_set
    )]

    if (length(var_set_results) != 0) {
      test <- unlist(var_set_results, recursive = FALSE)

      # Get clever covariate for each shift for each fold ----

      Hn_a_shift <- test[stringr::str_detect(names(test), "Hn_a_shift")]
      Hn_a_shift <- do.call(rbind, Hn_a_shift)

      Hn_a_z_shift <- test[stringr::str_detect(names(test), "Hn_az_shift")]
      Hn_a_z_shift <- do.call(rbind, Hn_a_z_shift)

      # Get scaled Qn for each shift for each fold ----

      Qn_a_shift <- test[stringr::str_detect(names(test), "Qn_a_shift_scaled")]
      Qn_a_shift <- do.call(rbind, Qn_a_shift)

      Qn_az_shift <- test[stringr::str_detect(names(test), "Qn_a_z_shift_scaled")]
      Qn_az_shift <- do.call(rbind, Qn_az_shift)


      data <- do.call(rbind, test[stringr::str_detect(
        names(test), "data"
      )])

      k_fold_results <- do.call(rbind, test[stringr::str_detect(
        names(test), "k_fold"
      )])

      deltas <- do.call(rbind, test[stringr::str_detect(
        names(test), "delta"
      )])

      tmle_fit_a_shift <- tmle_exposhift(
        data_internal = data,
        delta = mean(deltas),
        Qn_scaled = Qn_a_shift,
        Hn = Hn_a_shift,
        fluctuation = fluctuation,
        y = data$y
      )

      tmle_fit_az_shift <- tmle_exposhift(
        data_internal = data,
        delta = mean(deltas),
        Qn_scaled = Qn_az_shift,
        Hn = Hn_a_z_shift,
        fluctuation = fluctuation,
        y = data$y
      )

      var_names <- extract_vars_from_basis(
        var_set, 1,
        a_names, w_names, z_names
      )

      exposure <- stringr::str_extract(
        var_set,
        paste(c(a_names), collapse = "|")
      )
      mediator <- stringr::str_extract(var_set, paste(c(z_names),
        collapse = "|"
      ))
      exposure <- exposure[!is.na(exposure)]
      mediator <- mediator[!is.na(mediator)]

      mediation_in_fold <- calc_mediation_param(
        tmle_fit_a_shift = tmle_fit_a_shift,
        tmle_fit_a_z_shift = tmle_fit_az_shift,
        exposure,
        mediator,
        y = data$y,
        fold_k = "Pooled TMLE",
        delta = mean(deltas)
      )

      results_df <- rbind(k_fold_results, mediation_in_fold)
      rownames(results_df) <- NULL

      results_list[[var_set]] <- results_df
    } else {

    }
  }


  return(results_list)
}
