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
                           effect_truth,
                           deltas,
                           shift_var,
                           cv_folds){
  ## setup learners

  w <- data[, covars]
  a <- data[, exposures]
  y <- data[, outcome]

  sim_results <- SuperNOVA(
    w = w,
    a = a,
    y = y,
    var_sets = shift_var,
    delta = deltas,
    estimator = "tmle",
    fluctuation = "standard",
    n_folds = cv_folds,
    family = "continuous",
    quantile_thresh = 0,
    verbose = TRUE,
    parallel = TRUE,
    seed = seed
  )


  indiv_shift_results <- sim_results$`Indiv Shift Results`
  em_results <- sim_results$`Effect Mod Results`
  joint_shift_results <- sim_results$`Joint Shift Results`
  med_shift_results <- sim_results$`Mediation Shift Results`

  if (is.null(indiv_shift_results[[shift_var]]) == FALSE) {
    exposure_results <- indiv_shift_results[[shift_var]]
    exposure_results_pooled <- exposure_results[exposure_results$Fold == "Pooled TMLE", ]
    pooled_exposure_est <- exposure_results_pooled$Psi
    pooled_bias <- effect - pooled_exposure_est
    exposure_cov <- ifelse(
      (exposure_results_pooled$`Lower CI` <= effect &
        effect <= exposure_results_pooled$`Upper CI`), 1, 0
    )
  } else {
    pooled_exposure_est <- NULL
    pooled_bias <- NULL
    exposure_cov <- NULL
  }

  indiv_results <- list(
    "pooled_exposure_est" = pooled_exposure_est,
    "pooled_bias" = pooled_bias,
    "exposure_cov" = exposure_cov
  )


  if (is.null(em_results$a_2w_2) == FALSE) {
    a_2w_2_results <- em_results$a_2w_2
    a_2w_2_pooled <- a_2w_2_results[a_2w_2_results$Fold == "Pooled TMLE", ]

    a_2w_2_pooled_0 <- a_2w_2_pooled[1, ]
    a_2w_2_pooled_1 <- a_2w_2_pooled[2, ]


    a_2w_2_pooled_0_est <- a_2w_2_pooled_0$Psi
    a_2w_2_pooled_1_est <- a_2w_2_pooled_1$Psi

    a_2w_2_0_bias <- a_2w_2_pooled_0_est - a2_effect_in_w2_0
    a_2w_2_1_bias <- a_2w_2_pooled_1_est - a2_effect_in_w2_1

    a_2w_2_0_cov <- ifelse(
      (a_2w_2_pooled_0$`Lower CI` <= a2_effect_in_w2_0 &
        a2_effect_in_w2_0 <= a_2w_2_pooled_0$`Upper CI`), 1, 0
    )

    a_2w_2_1_cov <- ifelse(
      (a_2w_2_pooled_1$`Lower CI` <= a2_effect_in_w2_1 &
        a2_effect_in_w2_1 <= a_2w_2_pooled_1$`Upper CI`), 1, 0
    )

    em_results <- list(
      "a_2w_2_pooled_0_est" = a_2w_2_pooled_0_est,
      "a_2w_2_pooled_1_est" = a_2w_2_pooled_1_est,
      "a_2w_2_0_bias" = a_2w_2_0_bias,
      "a_2w_2_1_bias" = a_2w_2_1_bias,
      "a_2w_2_0_cov" = a_2w_2_0_cov,
      "a_2w_2_1_cov" = a_2w_2_1_cov
    )

  } else {
    em_results <-  NULL
  }

  if (is.null(names(med_shift_results)) == FALSE) {
    if (is.null(med_shift_results$a_1z) == FALSE) {
      a_1z_results <- med_shift_results$a_1z
      a_1z_pooled <- a_1z_results[a_1z_results$Fold == "Pooled TMLE", ]

      a_1z_pooled_NDE <- a_1z_pooled[1, ]
      a_1z_pooled_NIE <- a_1z_pooled[2, ]

      a_1z_nde <- a_1z_pooled_NDE$Estimate
      a_1z_nie <- a_1z_pooled_NIE$Estimate

      a_1z_nde_bias <- a_1z_nde - nde_truth
      a_1z_nie_bias <- a_1z_nie - nie_truth

      a_1z_nde_cov <- ifelse(
        (a_1z_pooled_NDE$`Lower CI` <= nde_truth &
          nde_truth <= a_1z_pooled_NDE$`Upper CI`), 1, 0
      )

      a_1z_nie_cov <- ifelse(
        (a_1z_pooled_NIE$`Lower CI` <= nie_truth &
          nie_truth <= a_1z_pooled_NIE$`Upper CI`), 1, 0
      )
    } else {
      a_1z_nde <- NULL
      a_1z_nie <- NULL

      a_1z_nde_bias <- NULL
      a_1z_nie_bias <- NULL

      a_1z_nde_cov <- NULL
      a_1z_nie_cov <- NULL
    }

    med_results <- list(
      "a_1z_nde" = a_1z_nde,
      "a_1z_nie" = a_1z_nie,
      "a_1z_nde_bias" = a_1z_nde_bias,
      "a_1z_nie_bias" = a_1z_nie_bias,
      "a_1z_nde_cov" = a_1z_nde_cov,
      "a_1z_nie_cov" = a_1z_nie_cov
    )

  }else{
    med_results <- NULL
  }

  # bundle estimators in list
  estimates <- c(indiv_results, em_results, med_results)

  return(estimates)
}
