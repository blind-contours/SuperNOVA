#' Calculate Final Results
#'
#' @details Set of functions that calculate the final target parameter
#' depending on single variable shift, effect modification shift,
#' or a joint stochastic intervention.
#'

#'
#' @importFrom assertthat assert_that
#'
#'


calc_final_ind_shift_param <- function(tmle_fit, exposure, fold_k) {
  condition <- exposure
  psi_param <- tmle_fit$psi - tmle_fit$noshift_psi
  variance_est <- var(tmle_fit$eif - tmle_fit$noshift_eif) / length(tmle_fit$eif)
  se_est <- sqrt(variance_est)
  CI <- calc_CIs(psi_param, se_est)

  Lower_CI <- CI[1]
  Upper_CI <- CI[2]
  p.value <- 2 * stats::pnorm(abs(psi_param / se_est), lower.tail = F)

  n <- length(tmle_fit$eif)

  results <- data.table::data.table(condition, psi_param, variance_est, se_est, Lower_CI, Upper_CI, p.value, fold_k, "Indiv Shift", exposure, n)

  names(results) <- c("Condition", "Psi", "Variance", "SE", "Lower CI", "Upper CI", "P-value", "Fold", "Type", "Variables", "N")

  return(results)
}

calc_final_effect_mod_param <- function(tmle_fit, exposure, effect_modifier, effect_m_name, fold_k) {

  inverse_prop_positive <- ifelse(effect_modifier == 1,  1/(table(effect_modifier)[[2]] / length(effect_modifier)), 0)
  inverse_prop_negative <- ifelse(effect_modifier == 0,  1/(table(effect_modifier)[[1]] / length(effect_modifier)), 0)

  inverse_prop_eif_pos <- inverse_prop_positive * tmle_fit$eif
  inverse_prop_eif_neg <- inverse_prop_negative * tmle_fit$eif

  psi_em_one <- mean(tmle_fit$Qn_shift_star[effect_modifier == 1])
  psi_em_zero <- mean(tmle_fit$Qn_shift_star[effect_modifier == 0])

  psi_diff <- psi_em_one - psi_em_zero

  psi_diff_eif <- inverse_prop_eif_pos - inverse_prop_eif_neg
  psi_diff_var <- var(psi_diff_eif) / length(psi_diff_eif)

  psi_one_var <- var(tmle_fit$eif[effect_modifier == 1]) / length(tmle_fit$eif)
  psi_zero_var <- var(tmle_fit$eif[effect_modifier == 0]) / length(tmle_fit$eif)

  psi_diff_se <- sqrt(psi_diff_var)

  psi_diff_CI <- calc_CIs(psi_diff, psi_diff_se)
  CIs <- mapply(calc_CIs, c(psi_em_one, psi_em_zero, psi_diff), c(sqrt(psi_one_var), sqrt(psi_zero_var), psi_diff_se))

  psi_p_val <- calc_pvals(psi_diff, psi_diff_se)
  level_1_p_val <-calc_pvals(psi_em_one, psi_one_var)
  level_0_p_val <- calc_pvals(psi_em_zero, psi_zero_var)

  psi_params <- c(psi_em_one, psi_em_zero, psi_diff)
  variance_ests <- c(psi_one_var, psi_zero_var, psi_diff_var)
  se_ests <- c(sqrt(psi_one_var), sqrt(psi_zero_var), psi_diff_se)

  Lower_CIs <- c(CIs[,1][1], CIs[,2][1], CIs[,3][1])
  Upper_CIs <- c(CIs[,1][2], CIs[,2][2], CIs[,3][2])

  P_val_ests <- c(level_1_p_val, level_0_p_val, psi_p_val)
  condition <- c(paste("Level 1 Shift Diff in ", effect_m_name), paste("Level 0 Shift Diff in", effect_m_name), "Effect Mod")

  n <- length(tmle_fit$eif)

  results <- data.table::data.table(condition, psi_params, variance_ests, se_ests, Lower_CIs, Upper_CIs, P_val_ests, fold_k, "Effect Mod", paste(exposure, effect_m_name, sep = ""), n)
  names(results) <- c("Condition", "Psi", "Variance", "SE", "Lower CI", "Upper CI", "P-value", "Fold", "Type", "Variables", "N")

  return(results)
}

calc_final_joint_shift_param <- function(joint_shift_fold_results, exposures, fold_k) {
  conditions <- exposures
  conditions[[3]] <- paste(conditions[[1]], conditions[[2]], sep = "&")
  conditions[[4]] <- "Psi"

  results <- lapply(joint_shift_fold_results, calc_joint_results)
  results_table <- do.call(rbind, results)

  intxn_results <- calc_intxn_results(results_table, joint_shift_fold_results)

  joint_intxn_results <- rbind(results_table,intxn_results)

  joint_intxn_results <- as.data.frame(cbind(conditions,joint_intxn_results, rep(fold_k, 4), "Interaction", conditions[[3]], length(joint_shift_fold_results[[3]]$eif)))

  names(joint_intxn_results) <- c("Condition","Psi", "Variance", "SE", "Lower CI", "Upper CI", "P-value", "Fold", "Type", "Variables", "N")
  rownames(joint_intxn_results) <- NULL
  return(joint_intxn_results)
}
