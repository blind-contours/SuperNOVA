fit_estimators <- function(w,
                           a,
                           y,
                           seed,
                           true_effects,
                           m14_effect_truth,
                           m14_intxn_truth,
                           true_em_effects,
                           deltas,
                           cv_folds,
                           var_sets){

  sim_results <- SuperNOVA(
    w = w,
    a = a,
    y = y,
    delta = deltas,
    estimator = "tmle",
    fluctuation = "standard",
    n_folds = cv_folds,
    family = "continuous",
    quantile_thresh = 0,
    verbose = TRUE,
    parallel = TRUE,
    seed = seed,
    var_sets = var_sets
  )


  indiv_shift_results <- sim_results$`Indiv Shift Results`
  em_results <- sim_results$`Effect Mod Results`
  joint_shift_results <- sim_results$`Joint Shift Results`

  indiv_biases <- vector()
  indiv_coverage <- vector()
  indiv_estimates <- vector()

  for (i in seq(true_effects)) {
    true_effect <- true_effects[i]
    exposure_results <- indiv_shift_results[[i]]
    exposure_results_pooled <- exposure_results[exposure_results$Fold == "Pooled TMLE", ]
    pooled_exposure_est <- exposure_results_pooled$Psi
    indiv_biases[[i]] <- true_effect - pooled_exposure_est
    indiv_coverage[[i]] <- ifelse(
      (exposure_results_pooled$`Lower CI` <= true_effect &
         true_effect <= exposure_results_pooled$`Upper CI`), 1, 0
    )

    indiv_estimates[[i]] <- pooled_exposure_est
  }


  indiv_results <- list(
    "indiv_est" = mean(indiv_estimates),
    "indiv_bias" = mean(indiv_biases),
    "indiv_cov" = mean(indiv_coverage)
  )


  M3W3_results <- em_results$M3W3
  M3W3_pooled <- M3W3_results[M3W3_results$Fold == "Pooled TMLE", ]

  level_1_em_results <- M3W3_pooled[str_detect(M3W3_pooled$Condition, "Level 1"),]
  level_0_em_results <- M3W3_pooled[str_detect(M3W3_pooled$Condition, "Level 0"),]

  level_1_truth <- true_em_effects[str_detect(names(true_em_effects), "Level 1")]
  level_0_truth <- true_em_effects[str_detect(names(true_em_effects), "Level 0")]

  est_level_1_bias <- mean(level_1_truth[[1]] - level_1_em_results$Psi)
  est_level_0_bias <- mean(level_0_truth[[1]] - level_0_em_results$Psi)

  est_biases <- mean(est_level_1_bias, est_level_0_bias)

  level_1_cov <- ifelse(
    (level_1_em_results$`Lower CI` <= level_1_truth[[1]] &
       level_1_truth[[1]] <= level_1_em_results$`Upper CI`), 1, 0
  )

  level_0_cov <- ifelse(
    (level_0_em_results$`Lower CI` <= level_0_truth[[1]] &
       level_0_truth[[1]] <= level_0_em_results$`Upper CI`), 1, 0
  )

  effect_mod_results <- list(
    "em_bias" = est_biases,
    "em_est" = mean(M3W3_pooled$Psi),
    "em_cov" = mean(level_1_cov, level_0_cov)
  )

  joint_intxn_results <- joint_shift_results$M1M4
  pooled_tmle_joint_intxn <- joint_intxn_results[joint_intxn_results$Fold == "Pooled TMLE", ]
  pooled_tmle_joint <- pooled_tmle_joint_intxn[pooled_tmle_joint_intxn$Condition == "M1&M4", ]
  pooled_tmle_intxn <- pooled_tmle_joint_intxn[pooled_tmle_joint_intxn$Condition == "Psi", ]

  intxn_bias <- pooled_tmle_intxn$Psi - m14_intxn_truth
  joint_bias <- pooled_tmle_joint$Psi - m14_effect_truth

  intxn_cov <- ifelse(
    (pooled_tmle_intxn$`Lower CI` <= m14_intxn_truth &
       m14_intxn_truth <= pooled_tmle_intxn$`Upper CI`), 1, 0
  )

  joint_cov <- ifelse(
    (pooled_tmle_joint$`Lower CI` <= m14_effect_truth &
       m14_effect_truth <= pooled_tmle_joint$`Upper CI`), 1, 0
  )

  joint_results <- list(
    "joint_bias" = joint_bias,
    "joint_est" = pooled_tmle_joint$Psi,
    "joint_cov" = joint_cov
  )

  intxn_results <- list(
    "intxn_bias" = intxn_bias,
    "intxn_est" = pooled_tmle_intxn$Psi,
    "intxn_cov" = intxn_cov
  )

  # bundle estimators in list
  estimates <- c(indiv_results, effect_mod_results, joint_results, intxn_results)

  return(estimates)
}


fit_estimators_mediation <- function(w,
                                     a,
                                     z,
                                     y,
                                     seed,
                                     nde_effects,
                                     nie_effects,
                                     ate_effects,
                                     deltas,
                                     cv_folds,
                                     var_sets){

  sim_results <- SuperNOVA(
    w = w,
    a = a,
    z = z,
    y = y,
    deltas = deltas,
    estimator = "tmle",
    fluctuation = "standard",
    n_folds = cv_folds,
    family = "continuous",
    quantile_thresh = 0,
    verbose = TRUE,
    parallel = TRUE,
    seed = seed,
    var_sets = var_sets
  )


  indiv_shift_results <- sim_results$`Indiv Shift Results`
  em_results <- sim_results$`Effect Mod Results`
  joint_shift_results <- sim_results$`Joint Shift Results`

  indiv_biases <- vector()
  indiv_coverage <- vector()
  indiv_estimates <- vector()

  for (i in seq(true_effects)) {
    true_effect <- true_effects[i]
    exposure_results <- indiv_shift_results[[i]]
    exposure_results_pooled <- exposure_results[exposure_results$Fold == "Pooled TMLE", ]
    pooled_exposure_est <- exposure_results_pooled$Psi
    indiv_biases[[i]] <- true_effect - pooled_exposure_est
    indiv_coverage[[i]] <- ifelse(
      (exposure_results_pooled$`Lower CI` <= true_effect &
         true_effect <= exposure_results_pooled$`Upper CI`), 1, 0
    )

    indiv_estimates[[i]] <- pooled_exposure_est
  }


  indiv_results <- list(
    "indiv_est" = mean(indiv_estimates),
    "indiv_bias" = mean(indiv_biases),
    "indiv_cov" = mean(indiv_coverage)
  )


  M3W3_results <- em_results$M3W3
  M3W3_pooled <- M3W3_results[M3W3_results$Fold == "Pooled TMLE", ]

  level_1_em_results <- M3W3_pooled[str_detect(M3W3_pooled$Condition, "Level 1"),]
  level_0_em_results <- M3W3_pooled[str_detect(M3W3_pooled$Condition, "Level 0"),]

  level_1_truth <- true_em_effects[str_detect(names(true_em_effects), "Level 1")]
  level_0_truth <- true_em_effects[str_detect(names(true_em_effects), "Level 0")]

  est_level_1_bias <- mean(level_1_truth[[1]] - level_1_em_results$Psi)
  est_level_0_bias <- mean(level_0_truth[[1]] - level_0_em_results$Psi)

  est_biases <- mean(est_level_1_bias, est_level_0_bias)

  level_1_cov <- ifelse(
    (level_1_em_results$`Lower CI` <= level_1_truth[[1]] &
       level_1_truth[[1]] <= level_1_em_results$`Upper CI`), 1, 0
  )

  level_0_cov <- ifelse(
    (level_0_em_results$`Lower CI` <= level_0_truth[[1]] &
       level_0_truth[[1]] <= level_0_em_results$`Upper CI`), 1, 0
  )

  effect_mod_results <- list(
    "em_bias" = est_biases,
    "em_est" = mean(M3W3_pooled$Psi),
    "em_cov" = mean(level_1_cov, level_0_cov)
  )

  joint_intxn_results <- joint_shift_results$M1M4
  pooled_tmle_joint_intxn <- joint_intxn_results[joint_intxn_results$Fold == "Pooled TMLE", ]
  pooled_tmle_joint <- pooled_tmle_joint_intxn[pooled_tmle_joint_intxn$Condition == "M1&M4", ]
  pooled_tmle_intxn <- pooled_tmle_joint_intxn[pooled_tmle_joint_intxn$Condition == "Psi", ]

  intxn_bias <- pooled_tmle_intxn$Psi - m14_intxn_truth
  joint_bias <- pooled_tmle_joint$Psi - m14_effect_truth

  intxn_cov <- ifelse(
    (pooled_tmle_intxn$`Lower CI` <= m14_intxn_truth &
       m14_intxn_truth <= pooled_tmle_intxn$`Upper CI`), 1, 0
  )

  joint_cov <- ifelse(
    (pooled_tmle_joint$`Lower CI` <= m14_effect_truth &
       m14_effect_truth <= pooled_tmle_joint$`Upper CI`), 1, 0
  )

  joint_results <- list(
    "joint_bias" = joint_bias,
    "joint_est" = pooled_tmle_joint$Psi,
    "joint_cov" = joint_cov
  )

  intxn_results <- list(
    "intxn_bias" = intxn_bias,
    "intxn_est" = pooled_tmle_intxn$Psi,
    "intxn_cov" = intxn_cov
  )

  # bundle estimators in list
  estimates <- c(indiv_results, effect_mod_results, joint_results, intxn_results)

  return(estimates)
}

