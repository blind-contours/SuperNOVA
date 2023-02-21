#' @title Calculates the Individual Shift Parameter
#' @description Simply subtract the psi estimates of a shift compared to no
#' shift. Use the delta method to do the same thing for the eif to derive
#' variance estimates for this new parameter.
#' @param tmle_fit TMLE results for the individual shift
#' @param exposure The exposure identified
#' @param fold_k Fold exposure is being shifted in
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
#' @param tmle_fit_av TMLE results for the validation fold
#' @param tmle_fit_at TMLE results for the training fold
#' @param exposure The exposure identified
#' @param at Training data
#' @param av Validation data
#' @param effect_m_name The effect modifier
#' @param fold_k The fold the effect modification was found
#' @param em_learner The Super Learner of decision trees
#' @importFrom partykit glmtree
#' @importFrom stats median
#' @importFrom dplyr mutate

#' @export

calc_final_effect_mod_param <- function(tmle_fit_av,
                                        tmle_fit_at,
                                        exposure,
                                        at,
                                        av,
                                        effect_m_name,
                                        fold_k,
                                        em_learner) {
  pseudo_outcome_at <- tmle_fit_at$qn_shift_star - tmle_fit_at$qn_noshift_star
  pseudo_outcome_av <- tmle_fit_av$qn_shift_star - tmle_fit_av$qn_noshift_star

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

    if (all(apply(at[, ..effect_m_name], 2, function(x) {
      all(x %in% 0:1)
    })) == TRUE) {
      medians <- apply(at[, ..effect_m_name], 2, median)

      median_rule_list <- list()
      for (i in seq(length(medians))) {
        median_rule <- paste(names(medians)[i], "==", medians[i])
        median_rule_list[[i]] <- median_rule
      }
      rules <- median_rule_list
    } else {
      medians <- apply(at[, ..effect_m_name], 2, median)
      median_rule_list <- list()
      for (i in seq(length(medians))) {
        median_rule <- paste(names(medians)[i], "<=", medians[i])
        median_rule_list[[i]] <- median_rule
      }

      rules <- median_rule_list
    }
  }

  results_list <- list()

  for (i in seq(rules)) {
    rule <- rules[i]

    em_split_data <- em_model_data_av %>%
      dplyr::mutate(ind = ifelse(eval(parse(text = rule)), 1, 0))

    # TODO: if the rule does not partition the validation sample and leads to
    # all one value then pass to next rule - think about better implementation.
    if (dim(table(em_split_data$ind)) == 1) {
      next()
    }

    inverse_prop_positive <- ifelse(em_split_data$ind == 1,
      1 / (table(em_split_data$ind)[[2]] /
        length(em_split_data$ind)), 0
    )

    inverse_prop_negative <- ifelse(em_split_data$ind == 0,
      1 / (table(em_split_data$ind)[[1]]
      / length(em_split_data$ind)), 0
    )

    inverse_prop_eif_pos <- inverse_prop_positive * (tmle_fit_av$eif - tmle_fit_av$noshift_eif)
    inverse_prop_eif_neg <- inverse_prop_negative * (tmle_fit_av$eif - tmle_fit_av$noshift_eif)

    diff <- tmle_fit_av$qn_shift_star - tmle_fit_av$qn_noshift_star

    psi_em_one <- mean(diff[em_split_data$ind == 1])
    psi_em_zero <- mean(diff[em_split_data$ind == 0])

    psi_one_var <- var(inverse_prop_eif_pos[em_split_data$ind == 1]) /
      table(em_split_data$ind)[[2]]

    psi_zero_var <- var(inverse_prop_eif_neg[em_split_data$ind == 0]) /
      table(em_split_data$ind)[[1]]

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
      paste("Level 1 Shift Diff in", rule),
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
  }

  results_df <- do.call(rbind, results_list)

  results_df <- results_df[!duplicated(results_df$Psi), ]

  return(results_df)
}

#' @title Calculates the Joint Shift Parameter
#' @description Estimates the shift parameter for a joint shift
#' @param joint_shift_fold_results Results of the joint shift
#' @param exposures Exposures shifted
#' @param fold_k Fold the joint shift is identified
#' @param deltas_updated The new delta, could be updated if Hn has positivity
#' violations
#' @export

calc_final_joint_shift_param <- function(joint_shift_fold_results,
                                         exposures,
                                         fold_k,
                                         deltas_updated) {
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
    length(joint_shift_fold_results[[3]]$eif),
    deltas_updated[[1]],
    deltas_updated[[2]]
  ))

  names(joint_intxn_results) <- c(
    "Condition", "Psi", "Variance",
    "SE", "Lower CI", "Upper CI",
    "P-value", "Fold", "Type",
    "Variables", "N", paste("Delta", conditions[[1]]),
    paste("Delta", conditions[[2]])
  )
  rownames(joint_intxn_results) <- NULL
  return(joint_intxn_results)
}

#' @title Calculates the Mediation Shift Parameters
#' @description Estimates the shift parameter for a natural direct and indirect
#' effect
#' @param tmle_fit_a_shift TMLE results for a shift in A alone
#' @param tmle_fit_a_z_shift TMLE results for a shift in A and Z
#' @param exposure The exposure found
#' @param mediator The mediator found
#' @param fold_k Fold the exposure and mediator were found
#' @param y Outcome
#' @param delta The shift amount


#' @export

calc_mediation_param <- function(tmle_fit_a_shift,
                                 tmle_fit_a_z_shift,
                                 exposure,
                                 mediator,
                                 fold_k,
                                 y,
                                 delta) {
  condition <- paste(exposure, mediator, sep = "&")

  nde_est <- tmle_fit_a_shift$psi - mean(y)
  nie_est <- tmle_fit_a_z_shift$psi - tmle_fit_a_shift$psi

  no_shift_eif <- tmle_fit_a_z_shift$noshift_eif

  nde_eif <- tmle_fit_a_shift$eif - no_shift_eif
  nie_eif <- tmle_fit_a_z_shift$eif - tmle_fit_a_shift$eif

  var_nde <- var(nde_eif) /
    length(nde_eif)
  se_nde <- sqrt(var_nde)
  ci_nde <- calc_CIs(nde_est, se_nde)
  lower_ci_nde <- ci_nde[1]
  upper_ci_nde <- ci_nde[2]
  p_value_nde <- 2 * stats::pnorm(abs(nde_est / se_nde), lower.tail = F)

  var_nie <- var(nie_eif) /
    length(nie_eif)
  se_nie <- sqrt(var_nie)
  ci_nie <- calc_CIs(nie_est, se_nie)
  lower_ci_nie <- ci_nie[1]
  upper_ci_nie <- ci_nie[2]
  p_value_nie <- 2 * stats::pnorm(abs(nie_est / se_nie), lower.tail = F)

  n <- length(nie_eif)

  nde_results <- data.table::data.table(
    "NDE", nde_est, var_nde, se_nde,
    lower_ci_nde, upper_ci_nde, p_value_nde, fold_k,
    "Mediation", n
  )

  nie_results <- data.table::data.table(
    "NIE", nie_est, var_nie, se_nie,
    lower_ci_nie, upper_ci_nie, p_value_nie, fold_k,
    "Mediation", n
  )

  names(nde_results) <- c(
    "Parameter", "Estimate", "Variance", "SE", "Lower CI",
    "Upper CI", "P-value", "Fold", "Type", "N"
  )


  names(nie_results) <- c(
    "Parameter", "Estimate", "Variance", "SE", "Lower CI",
    "Upper CI", "P-value", "Fold", "Type", "N"
  )

  nde_results$Delta <- delta
  nie_results$Delta <- delta

  mediation_results_table <- rbind(nde_results, nie_results)

  return(mediation_results_table)
}
