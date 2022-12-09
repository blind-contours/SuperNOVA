#' @title Calculates the Individual Shift Parameter
#' @description Simply subtract the psi estimates of a shift compared to no
#' shift. Use the delta method to do the same thing for the eif to derive
#' variance estimates for this new parameter.
#' @export
calc_final_ind_shift_param <- function(tmle_fit, exposure, fold_k) {
  condition <- exposure
  psi_param <- tmle_fit$psi - tmle_fit$noshift_psi
  variance_est <- var(tmle_fit$eif - tmle_fit$noshift_eif) /
    length(tmle_fit$eif)
  se_est <- sqrt(variance_est)
  CI <- calc_CIs(psi_param, se_est)

  Lower_CI <- CI[1]
  Upper_CI <- CI[2]
  p.value <- 2 * stats::pnorm(abs(psi_param / se_est), lower.tail = F)

  n <- length(tmle_fit$eif)

  results <- data.table::data.table(
    condition, psi_param, variance_est, se_est,
    Lower_CI, Upper_CI, p.value, fold_k,
    "Indiv Shift", exposure, n
  )

  names(results) <- c(
    "Condition", "Psi", "Variance", "SE", "Lower CI",
    "Upper CI", "P-value", "Fold", "Type", "Variables", "N"
  )

  return(results)
}

#' @title Calculates the Effect Modification Shift Parameter
#' @description Use a decision tree estimator to regress the eif of an exposure
#' on to a covariate - both of these have been data adaptively determined.
#' @importFrom partykit glmtree
#' @export

calc_final_effect_mod_param <- function(tmle_fit_av,
                                        tmle_fit_at,
                                        exposure,
                                        at,
                                        av,
                                        effect_m_name,
                                        fold_k,
                                        em_learner) {
  pseudo_outcome_at <- tmle_fit_at$eif - tmle_fit_at$noshift_eif
  pseudo_outcome_av <- tmle_fit_av$eif - tmle_fit_av$noshift_eif

  names(pseudo_outcome_at) <- "pseudo_outcome_at"
  names(pseudo_outcome_av) <- "pseudo_outcome_av"

  em_model_data_at <- cbind.data.frame(
    pseudo_outcome_at, subset(at,
      select = effect_m_name
    )
  )

  em_model_data_av <- cbind.data.frame(
    pseudo_outcome_av, subset(av,
      select = effect_m_name
    )
  )

  task <- sl3::make_sl3_Task(
    data = em_model_data_at,
    covariates = effect_m_name,
    outcome = "pseudo_outcome_at",
    outcome_type = "continuous"
  )

  discrete_sl_metalrn <- sl3::Lrnr_cv_selector$new(sl3::loss_squared_error)

  discrete_sl <- sl3::Lrnr_sl$new(
    learners = em_learner,
    metalearner = discrete_sl_metalrn,
  )

  sl_fit <- suppressWarnings(discrete_sl$train(task))

  selected_learner <- sl_fit$learner_fits[[
    which(sl_fit$coefficients == 1)
  ]]


  rules <- list_rules_party(selected_learner$fit_object)

  if (all(rules == "") == TRUE) {
    # create a rule for the median -----
    medians <- apply(at[, ..effect_m_name], 2, median)
    median_rule_list <- list()
    for (i in seq(length(medians))) {
      median_rule <- paste(names(medians)[i], "<=", medians[i])
      median_rule_list[[i]] <- median_rule
    }

    rules <- median_rule_list
  }

  results_list <- list()

  for (i in seq(rules)) {
    rule <- rules[i]

    em_split_data <- em_model_data_av %>%
      dplyr::mutate(ind = ifelse(eval(parse(text = rule)), 1, 0))

    inverse_prop_positive <- ifelse(em_split_data$ind == 1,
      1 / (table(em_split_data$ind)[[2]] /
        length(em_split_data$ind)), 0
    )

    inverse_prop_negative <- ifelse(em_split_data$ind == 0,
      1 / (table(em_split_data$ind)[[1]]
      / length(em_split_data$ind)), 0
    )

    inverse_prop_eif_pos <- inverse_prop_positive * tmle_fit_av$eif
    inverse_prop_eif_neg <- inverse_prop_negative * tmle_fit_av$eif

    psi_em_one <- mean(tmle_fit_av$qn_shift_star[em_split_data$ind == 1])
    psi_em_zero <- mean(tmle_fit_av$qn_shift_star[em_split_data$ind == 0])

    psi_one_var <- var(tmle_fit_av$eif[em_split_data$ind == 1]) /
      length(tmle_fit_av$eif)

    psi_zero_var <- var(tmle_fit_av$eif[em_split_data$ind == 0]) /
      length(tmle_fit_av$eif)

    em_one_ci <- calc_CIs(psi_em_one, sqrt(psi_one_var))
    em_zero_ci <- calc_CIs(psi_em_zero, sqrt(psi_zero_var))

    level_1_p_val <- calc_pvals(psi_em_one, psi_one_var)
    level_0_p_val <- calc_pvals(psi_em_zero, psi_zero_var)

    psi_params <- c(psi_em_one, psi_em_zero)
    variance_ests <- c(psi_one_var, psi_zero_var)
    se_ests <- c(sqrt(psi_one_var), sqrt(psi_zero_var))

    Lower_CIs <- c(em_one_ci[1], em_zero_ci[1])
    Upper_CIs <- c(em_one_ci[2], em_zero_ci[2])

    P_val_ests <- c(level_1_p_val, level_0_p_val)

    condition <- c(
      paste("Level 1 Shift Diff in ", rule),
      paste("Level 0 Shift Diff in", rule)
    )

    n <- length(tmle_fit_av$eif)

    results <- data.table::data.table(
      condition, psi_params, variance_ests,
      se_ests, Lower_CIs, Upper_CIs,
      P_val_ests, fold_k, "Effect Mod",
      paste(exposure, effect_m_name, sep = ""),
      n
    )
    names(results) <- c(
      "Condition", "Psi", "Variance", "SE", "Lower CI",
      "Upper CI", "P-value", "Fold", "Type", "Variables", "N"
    )

    results_list[[i]] <- results

    results_df <- do.call(rbind, results_list)
  }

  results_df <- results_df[!duplicated(results_df$Psi), ]

  return(results_df)
}

#' @title Calculates the Joint Shift Parameter
#' @description Estimates the shift parameter for a joint shift
#' @export

calc_final_joint_shift_param <- function(joint_shift_fold_results,
                                         exposures,
                                         fold_k) {
  conditions <- exposures
  conditions[[3]] <- paste(conditions[[1]], conditions[[2]], sep = "&")
  conditions[[4]] <- "Psi"

  results <- lapply(joint_shift_fold_results, calc_joint_results)
  results_table <- do.call(rbind, results)

  intxn_results <- calc_intxn_results(results_table, joint_shift_fold_results)

  joint_intxn_results <- rbind(
    as.data.frame(results_table),
    t(as.data.frame(unlist(intxn_results)))
  )

  joint_intxn_results <- as.data.frame(cbind(
    unlist(conditions),
    joint_intxn_results,
    rep(fold_k, 4),
    "Interaction",
    conditions[[3]],
    length(joint_shift_fold_results[[3]]$eif)
  ))

  names(joint_intxn_results) <- c(
    "Condition", "Psi", "Variance",
    "SE", "Lower CI", "Upper CI",
    "P-value", "Fold", "Type",
    "Variables", "N"
  )
  rownames(joint_intxn_results) <- NULL
  return(joint_intxn_results)
}
