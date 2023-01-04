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
                           mediators,
                           outcome,
                           seed,
                           nde_truth,
                           nie_truth,
                           ate_a_1,
                           ate_a_2,
                           a2_effect_in_w2_1,
                           a2_effect_in_w2_0,
                           deltas,
                           cv_folds) {
  ## setup learners

  w <- data[,covars]
  a <- data[,exposures]
  z <- data[,mediators]
  y <- data[,outcome]

  sim_results <- SuperNOVA(w = w,
                           a = a,
                           z = z,
                           y = y,
                           var_sets = c("a_1"),
                           delta = deltas,
                           estimator = "tmle",
                           fluctuation = "standard",
                           n_folds = cv_folds,
                           family = "continuous",
                           quantile_thresh = 0,
                           verbose = TRUE,
                           parallel = TRUE,
                           seed = seed)


  indiv_shift_results <- sim_results$`Indiv Shift Results`
  em_results <- sim_results$`Effect Mod Results`
  joint_shift_results <- sim_results$`Joint Shift Results`
  med_shift_results <- sim_results$`Mediation Shift Results`

  if (is.null(indiv_shift_results$a_1) == FALSE) {
    a1_results <- indiv_shift_results$a_1
    a1_pooled <- a1_results[a1_results$Fold == "Pooled TMLE",]
    a1_est <- a1_pooled$Psi
    a1_bias <- a1_est - ate_a_1
    a1_cov = ifelse(
      (a1_pooled$`Lower CI` <= ate_a_1 &
         ate_a_1 <= a1_pooled$`Upper CI`), 1,0
    )
  }else{
  a1_est <- NULL
  a1_bias <- NULL
  a1_cov <- NULL
  }

  if (is.null(indiv_shift_results$a_2) == FALSE) {
    a2_results <- indiv_shift_results$a_2
    a2_pooled <- a2_results[a2_results$Fold == "Pooled TMLE",]
    a2_est <- a2_pooled$Psi
    a2_bias <- a2_est - ate_a_2
    a2_cov = ifelse(
      (a2_pooled$`Lower CI` <= ate_a_2 &
         ate_a_2 <= a2_pooled$`Upper CI`), 1,0
    )
  }else{
    a2_est <- NULL
    a2_bias <- NULL
    a2_cov <- NULL
  }

  if (is.null(indiv_shift_results$a_3) == FALSE) {
    a3_results <- indiv_shift_results$a_3
    a3_pooled <- a3_results[a3_results$Fold == "Pooled TMLE",]
    a3_est <- a3_pooled$Psi
    a3_bias <- a3_est - ate_a_3
    a3_cov = ifelse(
      (a3_pooled$`Lower CI` <= ate_a_3 &
         ate_a_3 <= a3_pooled$`Upper CI`), 1,0
    )
  }else{
    a3_est <- NULL
    a3_bias <- NULL
    a3_cov <- NULL
  }

  indiv_results <- list("a1_est" = a1_est,
                        "a1_bias" = a1_bias,
                        "a1_cov" = a1_cov,
                        "a2_est" = a2_est,
                        "a2_bias" = a2_bias,
                        "a2_cov" = a2_cov,
                        "a3_est" = a3_est,
                        "a3_bias" = a3_bias,
                        "a3_cov" = a3_cov)


  if (is.null(em_results$a_2w_2) == FALSE) {
    a_2w_2_results <- em_results$a_2w_2
    a_2w_2_pooled <- a_2w_2_results[a_2w_2_results$Fold == "Pooled TMLE",]

    a_2w_2_pooled_0 <- a_2w_2_pooled[1,]
    a_2w_2_pooled_1 <- a_2w_2_pooled[2,]


    a_2w_2_pooled_0_est <- a_2w_2_pooled_0$Psi
    a_2w_2_pooled_1_est <- a_2w_2_pooled_1$Psi

    a_2w_2_0_bias <- a_2w_2_pooled_0_est - a2_effect_in_w2_0
    a_2w_2_1_bias <- a_2w_2_pooled_1_est - a2_effect_in_w2_1

    a_2w_2_0_cov = ifelse(
      (a_2w_2_pooled_0$`Lower CI` <= a2_effect_in_w2_0 &
         a2_effect_in_w2_0 <= a_2w_2_pooled_0$`Upper CI`), 1,0
    )

    a_2w_2_1_cov = ifelse(
      (a_2w_2_pooled_1$`Lower CI` <= a2_effect_in_w2_1 &
         a2_effect_in_w2_1 <= a_2w_2_pooled_1$`Upper CI`), 1,0
    )
  }else{
    a_2w_2_pooled_0_est <- NULL
    a_2w_2_pooled_1_est <- NULL

    a_2w_2_0_bias <- NULL
    a_2w_2_1_bias <- NULL

    a_2w_2_0_cov <- NULL
    a_2w_2_1_cov <- NULL
  }

  em_results <- list("a_2w_2_pooled_0_est" = a_2w_2_pooled_0_est,
                        "a_2w_2_pooled_1_est" = a_2w_2_pooled_1_est,
                        "a_2w_2_0_bias" = a_2w_2_0_bias,
                        "a_2w_2_1_bias" = a_2w_2_1_bias,
                        "a_2w_2_0_cov" = a_2w_2_0_cov,
                        "a_2w_2_1_cov" = a_2w_2_1_cov
                        )

  if (is.null(names(med_shift_results)) == FALSE) {
    if (is.null(med_shift_results$a_1z) == FALSE) {
      a_1z_results <- med_shift_results$a_1z
      a_1z_pooled <- a_1z_results[a_1z_results$Fold == "Pooled TMLE",]

      a_1z_pooled_NDE <- a_1z_pooled[1,]
      a_1z_pooled_NIE <- a_1z_pooled[2,]

      a_1z_nde <- a_1z_pooled_NDE$Estimate
      a_1z_nie <- a_1z_pooled_NIE$Estimate

      a_1z_nde_bias <- a_1z_nde - nde_truth
      a_1z_nie_bias <- a_1z_nie - nie_truth

      a_1z_nde_cov = ifelse(
        (a_1z_pooled_NDE$`Lower CI` <= nde_truth &
           nde_truth <= a_1z_pooled_NDE$`Upper CI`), 1,0
      )

      a_1z_nie_cov = ifelse(
        (a_1z_pooled_NIE$`Lower CI` <= nie_truth &
           nie_truth <= a_1z_pooled_NIE$`Upper CI`), 1,0
      )

    }else{
      a_1z_nde <- NULL
      a_1z_nie <- NULL

      a_1z_nde_bias <- NULL
      a_1z_nie_bias <- NULL

      a_1z_nde_cov <- NULL
      a_1z_nie_cov <- NULL
    }
  }


  med_results <- list("a_1z_nde" = a_1z_nde,
                       "a_1z_nie" = a_1z_nie,
                       "a_1z_nde_bias" = a_1z_nde_bias,
                       "a_1z_nie_bias" = a_1z_nie_bias,
                       "a_1z_nde_cov" = a_1z_nde_cov,
                       "a_1z_nie_cov" = a_1z_nie_cov
    )



  # bundle estimators in list
  estimates <- c(indiv_results, em_results, med_results)

  return(estimates)
}
