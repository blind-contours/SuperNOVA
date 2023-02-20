est_nde_nie <- function(data, m_learners, g_learners, e_learners,
                        psi_Z_learners) {
  node_list <- list(
    W = c("W_1", "W_2", "W_3"),
    A = "A",
    Z = c("Z_1", "Z_2", "Z_3"),
    Y = "Y"
  )
  learner_list <- list(
    Y = m_learners,
    A = g_learners
  )
  tmle_spec_NIE <- tmle_NIE(
    e_learners = e_learners,
    psi_Z_learners = psi_Z_learners,
    max_iter = 100
  )
  NIE_est <- tmle3(tmle_spec_NIE, data, node_list, learner_list)

  tmle_spec_NDE <- tmle_NDE(
    e_learners = e_learners,
    psi_Z_learners = psi_Z_learners,
    max_iter = 100
  )
  NDE_est <- tmle3(tmle_spec_NDE, data, node_list, learner_list)

  return(list(nie = NIE_est, nde = NDE_est))
}

###############################################################################
# fit estimators
###############################################################################
fit_estimators <- function(data,
                           covars,
                           exposures,
                           outcome,
                           seed,
                           true_effects,
                           m14_effect_truth,
                           m14_intxn_truth,
                           true_em_effects,
                           deltas,
                           shift_var,
                           cv_folds,
                           var_sets){
  ## setup learners

  w <- data[, covars]
  a <- data[, exposures]
  y <- data[, outcome]

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

  est_biases <- mean(true_em_effects[1] -  M3W3_pooled$Psi[2], true_em_effects[2] -  M3W3_pooled$Psi[1])

  level_1_cov <- ifelse(
    (M3W3_pooled$`Lower CI`[1] <= true_em_effects[2] &
       true_em_effects[2] <= M3W3_pooled$`Upper CI`[1]), 1, 0
  )

  level_2_cov <- ifelse(
    (M3W3_pooled$`Lower CI`[2] <= true_em_effects[1] &
       true_em_effects[1] <= M3W3_pooled$`Upper CI`[2]), 1, 0
  )

  effect_mod_results <- list(
    "em_bias" = est_biases,
    "em_est" = mean(M1W3_pooled$Psi),
    "em_cov" = mean(level_1_cov, level_2_cov)
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
