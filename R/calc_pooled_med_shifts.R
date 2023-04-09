#' Compute the Pooled Mediation Shift Parameter Estimate From the Fold
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
#' @importFrom dplyr bind_rows
#'
#' @return A \code{list} containing the parameter estimate, estimated variance
#'  based on the efficient influence function (EIF), the estimate of the EIF
#'  incorporating inverse probability of censoring weights, and the estimate of
#'  the EIF without the application of such weights.
calc_pooled_med_shifts <- function(med_shift_results,
                                   estimator = c("tmle", "onestep"),
                                   fluc_mod_out = NULL,
                                   a_names,
                                   w_names,
                                   z_names,
                                   y_name,
                                   fluctuation) {
  # set TMLE as default estimator type
  estimator <- match.arg(estimator)

  names <- names(med_shift_results)
  names <- gsub("^.+?:", "", names)
  names <- stringr::str_trim(names)

  results_list <- list()

  for (var_set in unique(names)) {
    var_set_results <- med_shift_results[stringr::str_detect(
      names(med_shift_results), var_set
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
        names(test), "delta"
      )])

      deltas <- mean(unlist(deltas))

      tmle_fit <- tmle_exposhift(
        data_internal = data,
        Qn_scaled = Qn_scaled,
        Hn = Hn,
        fluctuation = fluctuation,
        y = data$y,
        delta = mean(deltas)
      )

      total_effect <- tmle_fit$psi - mean(data$y)
      total_effect_eif <- tmle_fit$eif
      var_total_effect <- var(total_effect_eif) / length(total_effect_eif)

      total_effects <- list(
        "Parameter" = "Total-Pooled-TMLE",
        "Psi" = total_effect,
        "Variance" = var_total_effect,
        "SE" = sqrt(var_total_effect),
        "Lower CI" = calc_CIs(total_effect, sqrt(var_total_effect))[[1]],
        "Upper CI" = calc_CIs(total_effect, sqrt(var_total_effect))[[2]],
        "P-Value" = 2 * stats::pnorm(abs(total_effect / sqrt(var_total_effect)), lower.tail = F)
      )

      # NDE using pooled pseudo regression

      eif_comp_sum_w_pseudo <- as.vector(unlist(test[stringr::str_detect(names(test), "eif_comp_sum_w_pseudo")]))
      psi_comp_sum_w_pseudo <- mean(eif_comp_sum_w_pseudo)
      eif_comp_sum_w_pseudo <- eif_comp_sum_w_pseudo - psi_comp_sum_w_pseudo

      nde_effect_pseudo <- psi_comp_sum_w_pseudo - mean(data$y)

      eif_nde_no_shift <- test[stringr::str_detect(names(test), "eif_no_shift")]
      eif_nde_no_shift <- as.vector(unlist(eif_nde_no_shift))
      eif_nde_pseudo <- eif_comp_sum_w_pseudo - eif_nde_no_shift
      var_nde_pseudo <- var(eif_nde_pseudo) / length(eif_nde_pseudo)
      se_nde_pseudo <- sqrt(var_nde_pseudo)
      CI_nde_pseudo <- calc_CIs(nde_effect_pseudo, se_nde_pseudo)
      p_value_nde_pseudo <- 2 * stats::pnorm(abs(nde_effect_pseudo / se_nde_pseudo), lower.tail = F)

      nde_results_pseudo <- list(
        "Parameter" = "NDE-Pseudo-Reg",
        "Psi" = nde_effect_pseudo,
        "Variance" = var_nde_pseudo,
        "SE" = se_nde_pseudo,
        "Lower CI" = CI_nde_pseudo[[1]],
        "Upper CI" = CI_nde_pseudo[[2]],
        "P-Value" = p_value_nde_pseudo
      )

      # NDE using pooled integration

      eif_comp_sum_w_int <- as.vector(unlist(test[stringr::str_detect(names(test), "eif_comp_sum_w_double_int")]))
      psi_comp_sum_w_int <- mean(eif_comp_sum_w_int)
      eif_comp_sum_w_int <- eif_comp_sum_w_int - psi_comp_sum_w_int

      nde_effect_int <- psi_comp_sum_w_int - mean(data$y)

      eif_nde_no_shift <- test[stringr::str_detect(names(test), "eif_no_shift")]
      eif_nde_no_shift <- as.vector(unlist(eif_nde_no_shift))
      eif_nde_int <- eif_comp_sum_w_int - eif_nde_no_shift
      var_nde_int <- var(eif_nde_int) / length(eif_nde_int)
      se_nde_int <- sqrt(var_nde_int)
      CI_nde_int <- calc_CIs(nde_effect_int, se_nde_int)
      p_value_nde_int <- 2 * stats::pnorm(abs(nde_effect_int / se_nde_int), lower.tail = F)

      nde_results_int <- list(
        "Parameter" = "NDE-Integrated",
        "Psi" = nde_effect_int,
        "Variance" = var_nde_int,
        "SE" = se_nde_int,
        "Lower CI" = CI_nde_int[[1]],
        "Upper CI" = CI_nde_int[[2]],
        "P-Value" = p_value_nde_int
      )

      # NIE using pooled pseudo-regression

      nie_effect_pseudo <- total_effect - nde_effect_pseudo
      eif_nie_pseudo <- total_effect_eif - eif_nde_pseudo
      var_nie_pseudo <- var(eif_nie_pseudo) / length(eif_nie_pseudo)
      se_nie_pseudo <- sqrt(var_nie_pseudo)
      CI_nie_pseudo <- calc_CIs(nie_effect_pseudo, se_nie_pseudo)
      p_value_nie_pseudo <- 2 * stats::pnorm(abs(nie_effect_pseudo / se_nie_pseudo), lower.tail = F)

      nie_results_pseudo <- list(
        "Parameter" = "NIE-Pseudo-Reg",
        "Psi" = nie_effect_pseudo,
        "Variance" = var_nie_pseudo,
        "SE" = se_nie_pseudo,
        "Lower CI" = CI_nie_pseudo[[1]],
        "Upper CI" = CI_nie_pseudo[[2]],
        "P-Value" = p_value_nie_pseudo
      )


      # NIE using pooled integration

      nie_effect_int <- total_effect - nde_effect_int
      eif_nie_int <- total_effect_eif - eif_nde_int
      var_nie_int <- var(eif_nie_int) / length(eif_nie_int)
      se_nie_int <- sqrt(var_nie_int)
      CI_nie_int <- calc_CIs(nie_effect_int, se_nie_int)
      p_value_nie_int <- 2 * stats::pnorm(abs(nie_effect_int / se_nie_int), lower.tail = F)

      nie_results_int <- list(
        "Parameter" = "NIE-Integrated",
        "Psi" = nie_effect_int,
        "Variance" = var_nie_int,
        "SE" = se_nie_int,
        "Lower CI" = CI_nie_int[[1]],
        "Upper CI" = CI_nie_int[[2]],
        "P-Value" = p_value_nie_int
      )


      mediation_pooled <- bind_rows(
        nde_results_pseudo, nde_results_int,
        nie_results_pseudo, nie_results_int, total_effects
      )
      mediation_pooled$id <- "Pooled"

      rownames(mediation_pooled) <- NULL

      k_fold_results <- bind_rows(test[stringr::str_detect(names(test), "k_fold")], .id = "id")
      k_fold_results$id <- gsub(":.*$", "", k_fold_results$id)


      results_df <- bind_rows(k_fold_results, mediation_pooled)
      rownames(results_df) <- NULL

      results_list[[var_set]] <- results_df
    } else {

    }
  }


  return(results_list)
}
