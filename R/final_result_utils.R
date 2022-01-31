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


calc_final_ind_shift_param <- function(indiv_shift_fold_results, exposure, fold_k) {
  condition <- exposure
  psi_param <- indiv_shift_fold_results[[1]]$psi - indiv_shift_fold_results[[2]]$psi
  variance_est <- var(indiv_shift_fold_results[[1]]$eif - indiv_shift_fold_results[[2]]$eif) / length(indiv_shift_fold_results[[1]]$eif)
  se_est <- sqrt(variance_est)
  CI <- c(
    round(psi_param + stats::qnorm(0.05 / 2, lower.tail = T) * se_est, 4),
    round(psi_param + stats::qnorm(0.05 / 2, lower.tail = F) * se_est, 4)
  )

  Lower_CI <- CI[1]
  Upper_CI <- CI[2]
  p.value <- 2 * stats::pnorm(abs(psi_param / se_est), lower.tail = F)

  n <- length(indiv_shift_fold_results[[1]]$eif)

  results <- data.table::data.table(condition, psi_param, variance_est, se_est, Lower_CI, Upper_CI, p.value, fold_k, "Indiv Shift", exposure, n)
  names(results) <- c("Condition", "Psi", "Variance", "SE", "Lower CI", "Upper CI", "P-value", "Fold", "Type", "Variables", "N")

  return(results)
}

calc_final_effect_mod_param <- function(effect_mod_fold_results, exposure, effect_m_name, fold_k) {
  level_1_psi <- effect_mod_fold_results[[1]]$`1`$psi - effect_mod_fold_results[[2]]$`1`$psi
  level_0_psi <- effect_mod_fold_results[[1]]$`0`$psi - effect_mod_fold_results[[2]]$`0`$psi

  psi_diff <- level_1_psi - level_0_psi

  level_1_var <- effect_mod_fold_results[[1]]$`1`$var - effect_mod_fold_results[[2]]$`1`$var
  level_0_var <- effect_mod_fold_results[[1]]$`0`$var - effect_mod_fold_results[[2]]$`0`$var

  psi_se <- sqrt(sum(level_1_var, level_0_var))

  psi_CI <- c(
    round(psi_diff + stats::qnorm(0.05 / 2, lower.tail = T) * psi_se, 4),
    round(psi_diff + stats::qnorm(0.05 / 2, lower.tail = F) * psi_se, 4)
  )

  level_1_CI <- c(
    round(level_1_psi + stats::qnorm(0.05 / 2, lower.tail = T) * sqrt(level_1_var), 4),
    round(level_1_psi + stats::qnorm(0.05 / 2, lower.tail = F) * sqrt(level_1_var), 4)
  )

  level_0_CI <- c(
    round(level_0_psi + stats::qnorm(0.05 / 2, lower.tail = T) * sqrt(level_0_var), 4),
    round(level_0_psi + stats::qnorm(0.05 / 2, lower.tail = F) * sqrt(level_0_var), 4)
  )

  psi_p_val <- 2 * stats::pnorm(abs(psi_diff / psi_se), lower.tail = F)
  level_1_p_val <- 2 * stats::pnorm(abs(level_1_psi / sqrt(level_1_var)), lower.tail = F)
  level_0_p_val <- 2 * stats::pnorm(abs(level_0_psi / sqrt(level_0_var)), lower.tail = F)

  psi_params <- c(level_1_psi, level_0_psi, psi_diff)
  variance_ests <- c(level_1_var, level_0_var, sum(level_1_var, level_0_var))
  se_ests <- c(sqrt(level_1_var), sqrt(level_0_var), psi_se)

  Lower_CIs <- c(level_1_CI[1], level_0_CI[1], psi_CI[1])
  Upper_CIs <- c(level_1_CI[2], level_0_CI[2], psi_CI[2])

  P_val_ests <- c(level_1_p_val, level_0_p_val, psi_p_val)
  condition <- c(paste("Level 1 Shift Diff in ", effect_m_name), paste("Level 0 Shift Diff in", effect_m_name), "Effect Mod")

  n <- sum(length(effect_mod_fold_results[[1]]$`1`$eif), length(effect_mod_fold_results[[1]]$`0`$eif))

  results <- data.table::data.table(condition, psi_params, variance_ests, se_ests, Lower_CIs, Upper_CIs, P_val_ests, fold_k, "Effect Mod", paste(exposure, effect_m_name, sep = ""), n)
  names(results) <- c("Condition", "Psi", "Variance", "SE", "Lower CI", "Upper CI", "P-value", "Fold", "Type", "Variables", "N")

  return(results)
}

calc_final_joint_shift_param <- function(joint_shift_fold_results, exposures, fold_k) {
  conditions <- exposures
  conditions[[3]] <- paste(conditions[[1]], conditions[[2]], sep = "")
  conditions[[4]] <- "No Shift"
  conditions[[5]] <- "Psi"
  psi_target <- joint_shift_fold_results[[3]]$psi - joint_shift_fold_results[[1]]$psi - joint_shift_fold_results[[2]]$psi + joint_shift_fold_results[[4]]$psi
  psi_params <- c(joint_shift_fold_results[[1]]$psi, joint_shift_fold_results[[2]]$psi, joint_shift_fold_results[[3]]$psi, joint_shift_fold_results[[4]]$psi, psi_target)

  variance_ests <- c(joint_shift_fold_results[[1]]$var, joint_shift_fold_results[[2]]$var, joint_shift_fold_results[[3]]$var, joint_shift_fold_results[[4]]$var)
  psi_variance <- var(joint_shift_fold_results[[3]]$eif - joint_shift_fold_results[[1]]$eif - joint_shift_fold_results[[2]]$eif + joint_shift_fold_results[[4]]$eif) / length(joint_shift_fold_results[[3]]$eif)
  variance_ests <- c(variance_ests, psi_variance)

  psi_se <- sqrt(psi_variance)
  se_ests <- c(joint_shift_fold_results[[1]]$se, joint_shift_fold_results[[2]]$se, joint_shift_fold_results[[3]]$se, joint_shift_fold_results[[4]]$se, psi_se)

  psi_CI <- c(
    round(psi_target + stats::qnorm(0.05 / 2, lower.tail = T) * psi_se, 4),
    round(psi_target + stats::qnorm(0.05 / 2, lower.tail = F) * psi_se, 4)
  )

  Lower_CIs <- c(
    joint_shift_fold_results[[1]]$CI1,
    joint_shift_fold_results[[2]]$CI1,
    joint_shift_fold_results[[3]]$CI1,
    joint_shift_fold_results[[4]]$CI1,
    psi_CI[1]
  )

  Upper_CIs <- c(
    joint_shift_fold_results[[1]]$CI2,
    joint_shift_fold_results[[2]]$CI2,
    joint_shift_fold_results[[3]]$CI2,
    joint_shift_fold_results[[4]]$CI2,
    psi_CI[2]
  )


  psi_p_val <- 2 * stats::pnorm(abs(psi_target / psi_se), lower.tail = F)
  P_val_ests <- c(joint_shift_fold_results[[1]]$p_value, joint_shift_fold_results[[2]]$p_value, joint_shift_fold_results[[3]]$p_value, joint_shift_fold_results[[4]]$p_value, psi_p_val)
  results <- data.table::data.table(conditions, psi_params, variance_ests, se_ests, Lower_CIs, Upper_CIs, P_val_ests, rep(fold_k, 5), "Interaction", conditions[[3]], length(joint_shift_fold_results[[3]]$eif))
  names(results) <- c("Condition", "Psi", "Variance", "SE", "Lower CI", "Upper CI", "P-value", "Fold", "Type", "Variables", "N")

  return(results)
}
