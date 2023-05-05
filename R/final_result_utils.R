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

#' Calculates the Effect Modification Shift Parameter
#'
#' @description This function uses a decision tree estimator to regress the difference in expected outcomes
#'  of an exposure onto a covariate. Both the exposure and the covariate have been determined through data-adaptive methods.
#'  If no rules are found, a median split is applied.
#'
#' @param tmle_fit_av TMLE results for the validation fold.
#' @param tmle_fit_at TMLE results for the training fold.
#' @param exposure The identified exposure variable as a \code{character} string.
#' @param at Training dataset.
#' @param av Validation dataset.
#' @param effect_m_name The name of the effect modifier variable as a \code{character} string.
#' @param fold_k The fold in which the effect modification was found.
#' @importFrom partykit glmtree
#' @importFrom stats median
#' @importFrom dplyr mutate
#' @export
#'

calc_final_effect_mod_param <- function(tmle_fit_av,
                                        tmle_fit_at,
                                        exposure,
                                        at,
                                        av,
                                        effect_m_name,
                                        fold_k) {
  # Calculate pseudo_outcome for training and validation datasets
  pseudo_outcome_at <- tmle_fit_at$qn_shift_star - tmle_fit_at$qn_noshift_star
  pseudo_outcome_av <- tmle_fit_av$qn_shift_star - tmle_fit_av$qn_noshift_star

  # Rename columns
  names(pseudo_outcome_at) <- "pseudo_outcome_at"
  names(pseudo_outcome_av) <- "pseudo_outcome_av"

  # Prepare data for the tree model
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

  # Create the formula for the tree model
  response_var <-  "pseudo_outcome_at"
  formula_string <- paste(response_var, "~", effect_m_name)
  formula <- as.formula(formula_string)

  # Fit the tree model
  tree_model <- glmtree(formula, data = em_model_data_at)

  # Extract the rules from the tree model
  rules <- list_rules_party(tree_model)

  # If no rules are found, use median split
  if (all(rules == "") == TRUE) {
    # Determine if the effect modifier is binary (0 or 1)
    is_binary <- all(apply(at[, ..effect_m_name], 2, function(x) {
      all(x %in% 0:1)
    }))

    # Calculate medians for the effect modifier
    medians <- apply(at[, ..effect_m_name], 2, median)

    # Create median rules
    median_rule_list <- lapply(seq_along(medians), function(i) {
      if (is_binary) {
        paste(names(medians)[i], "==", medians[i])
      } else {
        paste(names(medians)[i], "<=", medians[i])
      }
    })

    rules <- median_rule_list
  }

  # Initialize results_list
  results_list <- list()

  # Loop through the rules
  for (i in seq_along(rules)) {
    rule <- rules[[i]]

    # Apply the rule to split the data
    em_split_data <- em_model_data_av %>%
      dplyr::mutate(ind = ifelse(eval(parse(text = rule)), 1, 0))

    # If the rule does not partition the validation sample and leads to all one value,
    # skip to the next rule
    if (length(unique(em_split_data$ind)) == 1) {
      next()
    }

    # Calculate the inverse propensity weights
    inverse_prop_positive <- ifelse(em_split_data$ind == 1,
                                    1 / (table(em_split_data$ind)[[2]] /
                                           length(em_split_data$ind)), 0
    )

    inverse_prop_negative <- ifelse(em_split_data$ind == 0,
                                    1 / (table(em_split_data$ind)[[1]]
                                         / length(em_split_data$ind)), 0
    )

    # Calculate the weighted EIFs
    inverse_prop_eif_pos <- inverse_prop_positive * (tmle_fit_av$eif - tmle_fit_av$noshift_eif)
    inverse_prop_eif_neg <- inverse_prop_negative * (tmle_fit_av$eif - tmle_fit_av$noshift_eif)

    # Calculate the shift difference for both groups
    diff <- tmle_fit_av$qn_shift_star - tmle_fit_av$qn_noshift_star
    psi_em_one <- mean(diff[em_split_data$ind == 1])
    psi_em_zero <- mean(diff[em_split_data$ind == 0])

    # Calculate variance and standard error
    psi_one_var <- var(inverse_prop_eif_pos[em_split_data$ind == 1]) /
      table(em_split_data$ind)[[2]]

    psi_zero_var <- var(inverse_prop_eif_neg[em_split_data$ind == 0]) /
      table(em_split_data$ind)[[1]]

    # Calculate confidence intervals
    em_one_ci <- calc_CIs(psi_em_one, sqrt(psi_one_var))
    em_zero_ci <- calc_CIs(psi_em_zero, sqrt(psi_zero_var))

    # Calculate p-values
    level_1_p_val <- calc_pvals(psi_em_one, psi_one_var)
    level_0_p_val <- calc_pvals(psi_em_zero, psi_zero_var)

    # Store results in a data frame
    results <- data.table::data.table(
      Condition = c(paste("Level 1 Shift Diff in", rule),
                    paste("Level 0 Shift Diff in", rule)),
      Psi = c(psi_em_one, psi_em_zero),
      Variance = c(psi_one_var, psi_zero_var),
      SE = c(sqrt(psi_one_var), sqrt(psi_zero_var)),
      Lower_CI = c(em_one_ci[1], em_zero_ci[1]),
      Upper_CI = c(em_one_ci[2], em_zero_ci[2]),
      P_value = c(level_1_p_val, level_0_p_val),
      Fold = fold_k,
      Type = "Effect Mod",
      Variables = paste(exposure, effect_m_name, sep = ""),
      N = length(tmle_fit_av$eif)
    )

    # Add the results to the results_list
    results_list[[i]] <- results
  }

  # Combine all results into a single data frame
  results_df <- do.call(rbind, results_list)

  # Remove duplicated results
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
