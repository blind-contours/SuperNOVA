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
  conditions[[3]] <- paste(conditions[[1]], conditions[[2]], sep = "&")
  conditions[[4]] <- "Psi"

  joint_psi <- joint_shift_fold_results[[3]]$psi - joint_shift_fold_results[[4]]$psi
  a_psi <- joint_shift_fold_results[[1]]$psi - joint_shift_fold_results[[4]]$psi
  b_psi <- joint_shift_fold_results[[2]]$psi - joint_shift_fold_results[[4]]$psi

  intxn_psi <- joint_psi - a_psi - b_psi

  psi_estimates <- c(a_psi, b_psi, joint_psi, intxn_psi)

  joint_variance <- var(joint_shift_fold_results[[3]]$eif - joint_shift_fold_results[[4]]$eif) / length(joint_shift_fold_results[[4]]$eif)
  a_psi_variance <- var(joint_shift_fold_results[[1]]$eif - joint_shift_fold_results[[4]]$eif) / length(joint_shift_fold_results[[4]]$eif)
  b_psi_variance <- var(joint_shift_fold_results[[2]]$eif - joint_shift_fold_results[[4]]$eif) / length(joint_shift_fold_results[[4]]$eif)

  # sum(joint_variance, a_psi_variance, b_psi_variance)

  psi_variance <- var((joint_shift_fold_results[[3]]$eif - joint_shift_fold_results[[4]]$eif) -
    (joint_shift_fold_results[[2]]$eif - joint_shift_fold_results[[4]]$eif) -
    (joint_shift_fold_results[[1]]$eif - joint_shift_fold_results[[4]]$eif)) / length(joint_shift_fold_results[[4]]$eif)

  variance_ests <- c(a_psi_variance, b_psi_variance, joint_variance, psi_variance)

  se_ests <- sapply(variance_ests, sqrt)

  calc_CI <- function(psi, psi_se) {
    psi_CI <- c(
      round(psi + stats::qnorm(0.05 / 2, lower.tail = T) * psi_se, 4),
      round(psi + stats::qnorm(0.05 / 2, lower.tail = F) * psi_se, 4)
    )
  }

  calc_p_value <- function(psi, psi_se) {
    2 * stats::pnorm(abs(psi / psi_se), lower.tail = F)
  }

  CIs <- mapply(calc_CI, psi_estimates, se_ests)
  p_vals <- mapply(calc_p_value, psi_estimates, se_ests)


  Lower_CIs <- CIs[1, ]
  Upper_CIs <- CIs[2, ]

  results <- data.table::data.table(conditions, psi_estimates, variance_ests, se_ests, Lower_CIs, Upper_CIs, p_vals, rep(fold_k, 4), "Interaction", conditions[[3]], length(joint_shift_fold_results[[3]]$eif))
  names(results) <- c("Condition", "Psi", "Variance", "SE", "Lower CI", "Upper CI", "P-value", "Fold", "Type", "Variables", "N")

  return(results)
}
