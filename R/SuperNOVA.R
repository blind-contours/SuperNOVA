#' @title Data-Adaptive Efficient Estimation of Interactions and Effect Modification using Stochastic Shift Interventions
#'
#' @details Treat variables identified in the same basis functions in ensemble b-spline models as
#' a data-adaptive parameter. Given variable sets identified, depending on if basis functions contain
#' variables for `A` or `AV` construct targeted minimum loss estimations of
#'  the counterfactual mean differences under various target parameters for individual variables, effect modifying variables,
#'  or interacting variables. Ensemble machine learning may be used to construct the initial
#'  estimates of nuisance functions using \pkg{sl3}.
#' @param W A \code{matrix}, \code{data.frame}, or similar containing a set of
#'  baseline covariates.
#' @param V \code{matrix}, \code{data.frame}, or similar containing a set of
#'  baseline covariates considered to be possible effect modifiers of exposure
#' @param A \code{matrix}, \code{data.frame}, or similar containing an individual or set of
#'  exposures
#' @param Y \code{numeric} vector of the observed outcomes.
#' @param delta A \code{numeric} value indicating the shift in the exposures to
#'  be used in defining the target parameter. This is defined with respect to
#'  the scale of the exposures (A).
#' @param LOD_vals vector of values indicating the limit of detection for each exposure variable
#' @param estimator The type of estimator to be fit, either \code{"tmle"} for
#'  targeted maximum likelihood or \code{"onestep"} for a one-step estimator.
#' @param fluctuation The method to be used in the submodel fluctuation step
#'  (targeting step) to compute the TML estimator. The choices are "standard"
#'  and "weighted" for where to place the auxiliary covariate in the logistic
#'  tilting regression.
#' @param max_iter A \code{numeric} integer giving the maximum number of steps
#'  to be taken in iterating to a solution of the efficient influence function.
#' @param sl_density_lrnr Learners used to fit Super Learner ensembles to densities via \pkg{sl3}
#' @param Q1_stack Learners used to fit Super Learner ensembles to the outcome model via \pkg{sl3}
#' @param n_folds Number of folds to use in cross-validation
#' @param family Outcome type family
#' @param quantile_thresh Threshold based on quantiles of the f-statistic used to identify "important" basis functions in the data-adaptive procedure
#' @param verbose Whether to run verbosely
#' @param parallel TRUE/FALSE parallelize across cores
#'
#' @return
#' @export
#' @importFrom MASS mvrnorm
#' @importFrom bindata rmvbin
#' @importFrom foreach %dopar%
#' @importFrom magrittr %>%
#' @importFrom stats as.formula glm p.adjust plogis predict qlogis qnorm qunif rnorm runif
#' @importFrom rlang :=
#' @importFrom data.table rbindlist

#' @examples
SuperNOVA <- function(W,
                      V,
                      A,
                      Y,
                      delta,
                      LOD_vals = NULL,
                      estimator = "tmle",
                      fluctuation = "standard",
                      max_iter = 10,
                      Density_stack,
                      Exposures_stack,
                      Covariate_stack,
                      Outcome_stack,
                      n_folds,
                      family,
                      quantile_thresh,
                      verbose = FALSE,
                      parallel = TRUE) {

  # check arguments and set up some objects for programmatic convenience
  call <- match.call(expand.dots = TRUE)
  estimator <- match.arg(estimator)
  fluctuation <- match.arg(fluctuation)
  # eif_reg_type <- match.arg(eif_reg_type)

  # coerce W to matrix and, if no names in W, assign them generically
  if (!is.data.frame(W)) W <- as.data.frame(W)
  W_names <- colnames(W)
  if (is.null(W_names)) {
    W_names <- paste0("W", seq_len(ncol(W)))
    colnames(W) <- W_names
  }

  if (!is.null(V)) {
    if (!is.data.frame(V)) V <- as.data.frame(V)
    V_names <- colnames(V)
    if (is.null(V_names)) {
      V_names <- paste0("V", seq_len(ncol(V)))
      colnames(V) <- V_names
    }
  } else {
    V_names <- NULL
  }

  # coerce W to matrix and, if no names in W, assign them generically
  A <- data.frame(A)
  A_names <- colnames(A)

  exposure_names <- c(V_names, A_names)

  if (is.null(A_names)) {
    A_names <- paste0("A", seq_len(ncol(A)))
    colnames(A) <- A_names
  }
  for (i in 1:dim(A)[2]) {
    check_A <- A[i]
    A_name <- A_names[i]
    if (length(unique(check_A)) < 20) {
      A[, A_name] <- A[, A_name] + runif(dim(A)[1], min = -0.05, max = 0.05)
    }
  }

  # coerce W to matrix and, if no names in W, assign them generically
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  colnames(Y) <- "Y"
  Y_names <- colnames(Y)

  # # dissociate fit type from other arguments to simplify passing to do.call
  # LOD_fit_type <- unlist(LOD_fit_args[names(LOD_fit_args) == "fit_type"],
  #   use.names = FALSE
  # )

  if (parallel == TRUE) {
    num_cores <- parallel::detectCores()
    doParallel::registerDoParallel(num_cores)
    future::plan(future::multisession)
  }


  V_in <- NULL

  if (!is.null(V)) {
    data_internal <- data.table::data.table(W, V, A, Y)
  } else {
    data_internal <- data.table::data.table(W, A, Y)
  }

  if (family == "binomial") {
    ## create the CV folds
    data_internal$folds <- create_cv_folds(n_folds, data_internal$Y)

    outcome_type <- "binomial"
  } else {
    # data_internal$Y <- scale_to_unit(data_internal$Y)
    # ## create the CV folds
    data_internal$folds <- create_cv_folds(n_folds, data_internal$Y)
    outcome_type <- "continous"
  }


  fold_basis_results <- foreach::foreach(fold_k = unique(data_internal$folds), .combine = "c") %dopar% {
    At <- data_internal[data_internal$folds != fold_k, ]
    Av <- data_internal[data_internal$folds == fold_k, ]

    covars <- W_names
    exposures <- c(A_names, V_names)

    basis_results <- fit_basis_estimators(
      At = At,
      covars = covars,
      exposures = exposures,
      outcome = Y_names,
      family = family,
      quantile_thresh = quantile_thresh,
      Exposures_stack = Exposures_stack,
      Covariate_stack = Covariate_stack,
      fold = fold_k,
      max_iter,
      verbose
    )

    basis_used <- basis_results$basis
    list(basis_used)
  }

  common_variables <- Reduce(intersect, fold_basis_results)

  fold_basis_results <- foreach::foreach(fold_k = unique(data_internal$folds), .combine = "rbind") %dopar% {

  # for (fold_k in unique(data_internal$folds)) {


    x <- c("Condition", "Psi", "Variance", "SE", "Lower CI", "Upper CI", "P-value", "Fold", "Type", "Variables", "N")
    within_fold_indiv_stoch_shift_df <- data.frame(matrix(ncol = 11, nrow = 0))
    within_fold_eff_mod_stoch_shift_df <- data.frame(matrix(ncol = 11, nrow = 0))
    within_fold_intxn_stoch_shift_df <- data.frame(matrix(ncol = 11, nrow = 0))

    colnames(within_fold_indiv_stoch_shift_df) <- x
    colnames(within_fold_eff_mod_stoch_shift_df) <- x
    colnames(within_fold_intxn_stoch_shift_df) <- x

    At <- data_internal[data_internal$folds != fold_k, ]
    Av <- data_internal[data_internal$folds == fold_k, ]

    for (i in 1:length(common_variables)) {

      hits <- extract_vars_from_basis(common_variables, i, A_names, V_names)

      matches <- hits$matches
      target <- hits$target

      if (length(matches) == 1 & any(grepl(paste(c(A_names, W_names, V_names), collapse = "|"), matches))) {
        exposure <- target
        covars <- c(W_names, V_names)
        ind_gn_exp_estim <- indiv_stoch_shift_est_g_exp(target, delta, Density_stack, covars, Av, At)

        # no_shift_gn_indiv <- data.table::data.table(cbind(ind_gn_exp_estim$noshift, ind_gn_exp_estim$noshift, ind_gn_exp_estim$noshift, ind_gn_exp_estim$noshift))
        # colnames(no_shift_gn_indiv) <- c("downshift", "noshift", "upshift", "upupshift")

        covars <- c(A_names, W_names, V_names)
        ind_Qn_estim <- indiv_stoch_shift_est_Q(exposure, delta, stack = Outcome_stack, covars, Av, At)

        # ind_no_shift_Q <- data.table::data.table(cbind(ind_Qn_estim$noshift, ind_Qn_estim$noshift, ind_Qn_estim$noshift, ind_Qn_estim$noshift))
        # colnames(ind_no_shift_Q) <- c("downshift", "noshift", "upshift", "upupshift")

        # ind_Hn_estims <- list()
        # ind_Qn_estims <- list()

        ind_Hn_estim <- est_Hn(gn_exp = ind_gn_exp_estim)
        # ind_no_shift_Hn_estims <- est_Hn(no_shift_gn_indiv)

        # ind_Hn_estims[[1]] <- ind_Hn_estim
        # ind_Qn_estims[[1]] <- ind_Qn_estim
        #
        # ind_Hn_estims[[2]] <- ind_no_shift_Hn_estims
        # ind_Qn_estims[[2]] <- ind_no_shift_Q

        # indiv_results_list <- list()

        # for (i in 1:length(ind_Qn_estims)) {
        #   Hn_estim <- ind_Hn_estims[[i]]
        #   Qn_estim <- ind_Qn_estims[[i]]

          tmle_fit <- tmle_exposhift(
            data_internal = Av,
            V = V_in,
            delta = delta,
            Qn_estim = ind_Qn_estim,
            Hn_estim = ind_Hn_estim,
            fluctuation = fluctuation,
            max_iter = max_iter,
            eif_reg_type = eif_reg_type,
            samp_fit_args = LOD_fit_args,
            ipcw_efficiency = ipcw_efficiency
          )


        tmle_fit$call <- call

        indiv_shift_in_fold <- calc_final_ind_shift_param(tmle_fit, exposure, fold_k)
        within_fold_indiv_stoch_shift_df <- rbind(within_fold_indiv_stoch_shift_df, indiv_shift_in_fold)

      } else if (length(matches) == 2 & any(grepl(paste(c(A_names), collapse = "|"), matches)) & any(grepl(paste(c(V_names), collapse = "|"), matches)) & !is.null(V_names)) {

        exposure <- stringr::str_extract(matches, paste(c(A_names), collapse = "|"))
        effect_m_name <- stringr::str_extract(matches, paste(c(V_names), collapse = "|"))
        exposure <- exposure[!is.na(exposure)]
        effect_m_name <- effect_m_name[!is.na(effect_m_name)]

        covars <- c(W_names)
        gn_exp_estim <- indiv_stoch_shift_est_g_exp(exposure, delta, Density_stack, covars, Av, At)
        Hn_estim <- est_Hn(gn_exp = gn_exp_estim)


        covars <- c(A_names, W_names)
        Qn_estim <- indiv_stoch_shift_est_Q(exposure, delta, stack = Outcome_stack, covars, Av, At)

        tmle_fit <- tmle_exposhift(
          data_internal = Av,
          V = V_in,
          delta = delta,
          Qn_estim = Qn_estim,
          Hn_estim = Hn_estim,
          fluctuation = fluctuation,
          max_iter = max_iter
        )

        if(length(unique(Av[[effect_m_name]])) > 2){
          V <- ifelse(Av[[effect_m_name]] >= quantile(Av[[effect_m_name]], 0.5)[[1]],1, 0)
        }else{
          V <- Av[[effect_m_name]]
        }

        effect_mod_in_fold <- calc_final_effect_mod_param(tmle_fit, exposure, effect_modifier = V, effect_m_name, fold_k)
        within_fold_eff_mod_stoch_shift_df <- rbind(within_fold_eff_mod_stoch_shift_df, effect_mod_in_fold)

      } else if (length(matches) == 2 & all(grepl(paste(c(A_names), collapse = "|"), matches))) {
        exposures <- as.list(matches)
        exposures[[3]] <- matches

        covars <- c(W_names, V_names)

        joint_gn_exp_estims <- joint_stoch_shift_est_g_exp(exposures, delta, Density_stack, covars, Av, At)

        joint_gn_exp_estims[[3]] <- joint_gn_exp_estims[[1]] * joint_gn_exp_estims[[3]]

        # joint_no_shift_joint <- data.table::data.table(cbind(joint_gn_exp_estims[[1]]$noshift, joint_gn_exp_estims[[1]]$noshift, joint_gn_exp_estims[[1]]$noshift, joint_gn_exp_estims[[1]]$noshift))
        # colnames(joint_no_shift_joint) <- c("downshift", "noshift", "upshift", "upupshift")

        joint_Hn_estims <- lapply(joint_gn_exp_estims, est_Hn)
        # joint_no_shift_Hn_estims <- est_Hn(joint_no_shift_joint)

        covars <- c(A_names, W_names, V_names)
        joint_Qn_estims <- joint_stoch_shift_est_Q(exposures, delta, stack = Outcome_stack, covars, Av, At)

        # joint_no_shift_Q <- data.table::data.table(cbind(joint_Qn_estims[[1]]$noshift, joint_Qn_estims[[1]]$noshift, joint_Qn_estims[[1]]$noshift, joint_Qn_estims[[1]]$noshift))
        # colnames(joint_no_shift_Q) <- c("downshift", "noshift", "upshift", "upupshift")

        # joint_Hn_estims[[4]] <- joint_no_shift_Hn_estims
        # joint_Qn_estims[[4]] <- joint_no_shift_Q

        intxn_results_list <- list()

        for (i in 1:length(joint_Qn_estims)) {
          Hn_estim <- joint_Hn_estims[[i]]
          Qn_estim <- joint_Qn_estims[[i]]

          tmle_fit <- tmle_exposhift(
            data_internal = Av,
            V = V_in,
            delta = delta,
            Qn_estim = Qn_estim,
            Hn_estim = Hn_estim,
            fluctuation = fluctuation,
            max_iter = max_iter,
            eif_reg_type = eif_reg_type,
            samp_fit_args = LOD_fit_args,
            ipcw_efficiency = ipcw_efficiency
          )

          intxn_results_list[[i]] <- tmle_fit
        }

        intxn_in_fold <- calc_final_joint_shift_param(intxn_results_list, matches, fold_k)
        within_fold_intxn_stoch_shift_df <- rbind(within_fold_intxn_stoch_shift_df, intxn_in_fold)
      }
    }

    results <- rbind(within_fold_indiv_stoch_shift_df, within_fold_eff_mod_stoch_shift_df, within_fold_intxn_stoch_shift_df)

    results
  }


  groups <-
    fold_basis_results %>% dplyr::group_by(Type)

  results_list <- dplyr::group_split(groups)

  for (i in seq(results_list)) {
    assign(unique(results_list[[i]]$Type), results_list[[i]])
  }

  for (target_param in c("Indiv Shift", "Effect Mod", "Interaction")) {
    if (any(grepl(target_param, ls())) == FALSE) {
      assign(target_param, NA)
    }
  }

  results_list <- list(
    "Effect Mod Results" = `Effect Mod`,
    "Indiv Shift Results" = `Indiv Shift`,
    "Joint Shift Results" = `Interaction`
  )

  return(results_list)
}
