#' @title Data-adaptive efficient estimation of interactions and effect
#' modification using stochastic shift interventions
#'
#' @details Treat variables identified in the same basis functions in ensemble
#' b-spline models as a data-adaptive parameter. Given variable sets identified,
#'  depending on if basis functions contain variables for `A` or `AV` construct
#'  targeted minimum loss estimations of the counterfactual mean differences
#'  under various target parameters for individual variables, effect modifying
#'  variables, or interacting variables. Ensemble machine learning may be used
#'  to construct the initial estimates of nuisance functions using \pkg{sl3}.
#' @param w A \code{matrix}, \code{data.frame}, or similar containing a set of
#'  baseline covariates.
#' @param a \code{matrix}, \code{data.frame}, or similar containing an
#' individual or set of
#'  exposures
#' @param z \code{matrix}, \code{data.frame}, or similar containing an
#' individual or set of mediators
#' @param y \code{numeric} vector of the observed outcomes.
#' @param deltas A \code{numeric} value indicating the shift in the exposures to
#'  be used in defining the target parameter. This is defined with respect to
#'  the scale of the exposures (A).
#' @param var_sets If using SuperNOVA deterministically, this parameter takes
#' in a list of the form var_sets <- c("A_1", "A_1-Z_2") etc. where the
#' analyst passes in variable sets for exposures, exosure-mediator, or exposure-
#' covariate.
#' @param estimator The type of estimator to be fit, either \code{"tmle"} for
#'  targeted maximum likelihood or \code{"onestep"} for a one-step estimator.
#' @param fluctuation The method to be used in the submodel fluctuation step
#'  (targeting step) to compute the TML estimator. The choices are "standard"
#'  and "weighted" for where to place the auxiliary covariate in the logistic
#'  tilting regression.
#' @param pi_learner Learners used to fit Super Learner ensembles to densities
#' via \pkg{sl3}
#' @param mu_learner Learners used to fit Super Learner ensembles to the outcome
#' model via \pkg{sl3}
#' @param g_learner Learners used to fit Super Learner ensembles to the g mech
#' g(A|W) - a probability estimator not a density estimator used in mediation
#' via \pkg{sl3}
#' @param e_learner Learners used to fit Super Learner ensembles to the e mech
#' g(A|Z,W) - a probability estimator not a density estimator used in mediation
#' via \pkg{sl3}
#' @param zeta_learner Learners used to fit Super Learner ensembles to the outcome
#' model via \pkg{sl3}
#' @param em_learner Super learner from \pkg{sl3} of decision trees used to
#' apply to the EIF of the shift applied to an exposure using covariates W.
#' @param n_folds Number of folds to use in cross-validation
#' @param family Outcome type family
#' @param quantile_thresh Threshold based on quantiles of the f-statistic used
#' to identify "important" basis functions in the data-adaptive procedure
#' @param verbose Whether to run verbosely
#' @param parallel TRUE/FALSE parallelize across cores
#' @param seed \code{numeric} seed value to be passed to all functions
#' @param hn_trunc_thresh Truncation level for the clever covariate
#' @param parallel_type If parallel is TRUE, type of parallelization to use. Default
#' is "multi_session", other values are multicore and sequential.
#' @param num_cores Number of CPU cores to use in parallelization
#' @param adaptive_delta If TRUE, this reduces the user specified delta until
#' the Hn calculated for a shift does not have any observation that is greater
#' than hn_trunc_thresh.
#'
#' @export
#' @importFrom MASS mvrnorm
#' @importFrom foreach %dopar%
#' @importFrom magrittr %>%
#' @importFrom stats as.formula glm p.adjust plogis predict qlogis qnorm qunif rnorm runif
#' @importFrom rlang :=
#' @importFrom stringr str_count
#' @import furrr
#' @import purrr
#' @importFrom data.table rbindlist
#' @return S3 object of class \code{SuperNOVA} containing the results of the
#'  procedure to compute a TML or one-step estimate of the counterfactual mean
#'  under a modified treatment policy that shifts a continuous-valued exposure
#'  by a scalar amount \code{delta}. These exposure are data-adaptively
#'  identified using the CV-TMLE procedure.
#' )

SuperNOVA <- function(w,
                      a,
                      z = NULL,
                      y,
                      deltas,
                      estimator = "tmle",
                      fluctuation = "standard",
                      var_sets = NULL,
                      pi_learner = NULL,
                      mu_learner = NULL,
                      g_learner = NULL,
                      e_learner = NULL,
                      em_learner = NULL,
                      zeta_learner = NULL,
                      n_folds = 2,
                      outcome_type = "continuous",
                      quantile_thresh = 0.5,
                      verbose = FALSE,
                      parallel = TRUE,
                      parallel_type = "multi_session",
                      num_cores = 2,
                      seed = seed,
                      hn_trunc_thresh = 10,
                      adaptive_delta = FALSE,
                      n_mc_sample = 1000,
                      exposure_quantized = FALSE,
                      mediator_quantized = FALSE,
                      density_type = "hal",
                      n_bins = 10,
                      max_degree = 1,
                      integration_method = "MC") {
  # check arguments and set up some objects for programmatic convenience
  call <- match.call(expand.dots = TRUE)
  estimator <- match.arg(estimator)
  fluctuation <- match.arg(fluctuation)

  # coerce W to matrix and, if no names in W, assign them generically
  if (!is.data.frame(w)) w <- as.data.frame(w)
  w_names <- colnames(w)
  if (is.null(w_names)) {
    w_names <- paste0("w", seq_len(ncol(w)))
    colnames(w) <- w_names
  }

  # coerce Z to matrix and, if no names in Z, assign them generically
  if (!is.data.frame(z)) z <- as.data.frame(z)
  z_names <- colnames(z)
  if (is.null(z_names)) {
    z_names <- paste0("z", seq_len(ncol(z)))
    colnames(z) <- z_names
  }

  # coerce W to matrix and, if no names in W, assign them generically
  a <- data.frame(a)
  a_names <- colnames(a)

  if (is.null(a_names)) {
    a_names <- paste0("a", seq_len(ncol(a)))
    colnames(a) <- a_names
  }
  # for (i in 1:dim(a)[2]) {
  #   check_a <- a[i]
  #   a_name <- a_names[i]
  #   if (length(unique(check_a)) < 20) {
  #     a[, a_name] <- a[, a_name] + runif(dim(a)[1], min = -0.05, max = 0.05)
  #   }
  # }

  if (is.null(pi_learner)) {
    sls <- create_sls()
    pi_learner <- sls$pi_learner
  }

  if (is.null(mu_learner)) {
    sls <- create_sls()
    mu_learner <- sls$mu_learner
  }

  if (is.null(zeta_learner)) {
    sls <- create_sls()
    zeta_learner <- sls$zeta_learner
  }

  if (is.null(g_learner)) {
    sls <- create_sls()
    g_learner <- sls$g_learner
  }

  if (is.null(e_learner)) {
    sls <- create_sls()
    e_learner <- sls$e_learner
  }

  if (is.null(em_learner)) {
    sls <- create_sls()
    em_learner <- sls$em_learner
  }

  if (exposure_quantized == TRUE) {
    sls <- create_sls()
    quant_learner <- sls$quant_learner
  }

  if (parallel == TRUE) {
    if (parallel_type == "multi_session") {
      future::plan(future::multisession,
        workers = num_cores,
        gc = TRUE
      )
    } else {
      future::plan(future::multicore,
        workers = num_cores,
        gc = TRUE
      )
    }
  } else {
    future::plan(future::sequential,
      gc = TRUE
    )
  }

  data_internal <- data.table::data.table(w, a, z, y)
  `%notin%` <- Negate(`%in%`)


  if (outcome_type == "binomial") {
    ## create the CV folds
    data_internal$folds <- create_cv_folds(n_folds, data_internal$y)

    outcome_type <- "binomial"
  } else {
    data_internal$folds <- create_cv_folds(n_folds, data_internal$y)
    outcome_type <- "continous"
  }

  if (is.null(var_sets)) {
    fold_basis_results <- furrr::future_map(unique(data_internal$folds),
      function(fold_k) {
        at <- data_internal[data_internal$folds != fold_k, ]
        av <- data_internal[data_internal$folds == fold_k, ]


        basis_results <- fit_basis_estimators(
          at = at,
          a_names = a_names,
          z_names = z_names,
          w_names = w_names,
          outcome = "y",
          family = family,
          quantile_thresh = quantile_thresh,
          zeta_learner = zeta_learner,
          fold = fold_k,
          seed = seed
        )

        basis_used <- basis_results$basis
      },
      .options = furrr::furrr_options(seed = seed, packages = "SuperNOVA")
    )

    basis_prop_in_fold <- calc_basis_freq(fold_basis_results, n_folds)

    if (length(basis_prop_in_fold) == 0) {
      print("No Exposures were Predictive")
      return(NULL)
    }
  } else {
    fold_basis_results <- rep(list(var_sets), n_folds)
    basis_prop_in_fold <- calc_basis_freq(fold_basis_results, n_folds)
  }

  fold_results_indiv <- list()
  fold_results_em <- list()
  fold_results_intxn <- list()
  fold_results_mediation <- list()
  joint_fold_results_mediation <- list()

  fold_SuperNOVA_results <- furrr::future_map(
    unique(data_internal$folds), function(fold_k) {
      at <- data_internal[data_internal$folds != fold_k, ]
      av <- data_internal[data_internal$folds == fold_k, ]

      basis_variables <- fold_basis_results[[fold_k]]

      basis_variables <- basis_variables[
        basis_variables %notin% w_names
      ]

      for (i in 1:length(basis_variables)) {
        hits <- extract_vars_from_basis(
          basis_variables, i,
          a_names, w_names, z_names
        )

        matches <- hits$matches
        target <- hits$target

        # Check if basis is one exposure ---------------------------

        if (length(matches) == 1 & any(grepl(paste(c(a_names),
          collapse = "|"
        ), matches))) {
          exposure <- target

          delta <- deltas[[exposure]]

          lower_bound <- min(min(av[[exposure]]), min(at[[exposure]]))
          upper_bound <- max(max(av[[exposure]]), max(at[[exposure]]))

          if (exposure_quantized == TRUE) {
            g_type <- "categorical"
          } else {
            g_type <- "continuous"
          }

          ind_gn_exp_estim <- indiv_stoch_shift_est_g_exp(
            exposure = target,
            delta = delta,
            g_learner = pi_learner,
            covars = w_names,
            av = av,
            at = at,
            adaptive_delta = adaptive_delta,
            hn_trunc_thresh = hn_trunc_thresh,
            exposure_quantized = exposure_quantized,
            lower_bound = lower_bound,
            upper_bound = upper_bound,
            outcome_type = g_type
          )

          delta <- ind_gn_exp_estim$delta

          covars <- c(a_names, w_names)

          ind_qn_estim <- indiv_stoch_shift_est_Q(
            exposure = exposure,
            delta = delta,
            mu_learner = mu_learner,
            covars = covars,
            av = av,
            at = at,
            lower_bound = lower_bound,
            upper_bound = upper_bound
          )

          Hn <- ind_gn_exp_estim$Hn_av

          tmle_fit <- tmle_exposhift(
            data_internal = av,
            delta = delta,
            Qn_scaled = ind_qn_estim$q_av,
            Hn = Hn,
            fluctuation = fluctuation,
            y = av$y
          )

          tmle_fit$call <- call

          indiv_shift_in_fold <- calc_final_ind_shift_param(
            tmle_fit,
            exposure,
            fold_k
          )

          indiv_shift_in_fold$Delta <- delta

          fold_results_indiv[[
            paste("Fold", fold_k, ":", target)
          ]] <- list(
            "data" = av,
            "Qn_scaled" = ind_qn_estim$q_av,
            "Hn" = Hn,
            "k_fold_result" = indiv_shift_in_fold,
            "Delta" = delta
          )

          # Check if basis is one exposure and baseline covariates ---------
        } else if (sum(stringr::str_count(matches, paste(c(a_names), collapse = "|"))) == 1 &
          sum(stringr::str_count(matches, paste(c(w_names), collapse = "|"))) >= 1) {
          exposure <- stringr::str_extract(
            matches,
            paste(c(a_names), collapse = "|")
          )
          effect_m_name <- stringr::str_extract(matches, paste(c(w_names),
            collapse = "|"
          ))
          exposure <- exposure[!is.na(exposure)]
          effect_m_name <- effect_m_name[!is.na(effect_m_name)]

          delta <- deltas[[exposure]]
          covars <- c(w_names)

          gn_exp_estim <- indiv_stoch_shift_est_g_exp(
            exposure = exposure,
            delta = delta,
            g_learner = pi_learner,
            covars = covars,
            av = av,
            at = at,
            adaptive_delta = adaptive_delta,
            hn_trunc_thresh = hn_trunc_thresh
          )

          delta <- gn_exp_estim$delta

          covars <- c(a_names, w_names, z_names)

          qn_estim <- indiv_stoch_shift_est_Q(
            exposure = exposure,
            delta = delta,
            mu_learner = mu_learner,
            covars = covars,
            av = av,
            at = at
          )

          Hn_av <- gn_exp_estim$Hn_av
          Hn_at <- gn_exp_estim$Hn_at

          tmle_fit_av <- tmle_exposhift(
            data_internal = av,
            delta = delta,
            Qn_scaled = qn_estim$q_av,
            Hn = Hn_av,
            fluctuation = fluctuation,
            y = av$y
          )

          tmle_fit_at <- tmle_exposhift(
            data_internal = at,
            delta = delta,
            Qn_scaled = qn_estim$q_at,
            Hn = Hn_at,
            fluctuation = fluctuation,
            y = at$y
          )

          effect_mod_in_fold <- calc_final_effect_mod_param(
            tmle_fit_av = tmle_fit_av,
            tmle_fit_at = tmle_fit_at,
            exposure = exposure,
            at = at,
            av = av,
            effect_m_name = effect_m_name,
            fold_k = fold_k,
            em_learner = em_learner
          )

          effect_mod_in_fold$Delta <- delta

          fold_results_em[[
            paste("Fold", fold_k, ":", target)
          ]] <- list(
            "av_data" = av,
            "at_data" = at,
            "Qn_scaled_av" = qn_estim$q_av,
            "Qn_scaled_at" = qn_estim$q_at,
            "Hn_av" = Hn_av,
            "Hn_at" = Hn_at,
            "k_fold_result" = effect_mod_in_fold,
            "delta" = delta
          )


          # Check if basis is two exposures ---------------------------
        } else if (length(matches) == 2 & all(grepl(
          paste(c(a_names),
            collapse = "|"
          ),
          matches
        ))) {
          exposures <- as.list(matches)
          delta <- deltas[unlist(exposures)]
          exposures[[3]] <- matches

          covars <- c(w_names)

          joint_gn_exp_estims <- joint_stoch_shift_est_g_exp(
            exposures,
            deltas,
            g_learner = pi_learner,
            covars,
            av,
            at,
            adaptive_delta = adaptive_delta,
            hn_trunc_thresh = hn_trunc_thresh,
            upper_bound = upper_bound,
            lower_bound = lower_bound
          )

          joint_gn_exp_estims$gn_results[[3]] <- mapply(
            `*`,
            joint_gn_exp_estims$gn_results[[1]],
            joint_gn_exp_estims$gn_results[[3]]
          )


          deltas_updated <- joint_gn_exp_estims$delta_results
          deltas_updated[[3]] <- c(deltas_updated[[1]], deltas_updated[[2]])

          covars <- c(a_names, w_names)

          joint_qn_estims <- joint_stoch_shift_est_Q(
            exposures,
            deltas = deltas_updated,
            mu_learner = mu_learner,
            covars,
            av,
            at,
            upper_bound,
            lower_bound
          )

          intxn_results_list <- list()
          qn_estim_scaled_list <- list()
          joint_hn_estims <- list()

          for (i in 1:length(joint_qn_estims)) {
            hn_estim <- joint_gn_exp_estims$Hn_results[[i]]

            qn_estim <- joint_qn_estims[[i]]
            delta <- deltas_updated[[i]]

            tmle_fit <- tmle_exposhift(
              data_internal = av,
              delta = delta,
              Qn_scaled = qn_estim,
              Hn = hn_estim,
              fluctuation = fluctuation,
              y = av$y
            )

            intxn_results_list[[i]] <- tmle_fit
            joint_hn_estims[[i]] <- hn_estim
          }

          intxn_in_fold <- calc_final_joint_shift_param(
            joint_shift_fold_results = intxn_results_list,
            matches,
            fold_k,
            deltas_updated
          )

          fold_results_intxn[[
            paste("Fold", fold_k, ":", target)
          ]] <- list(
            "data" = av,
            "Qn_scaled" = joint_qn_estims,
            "Hn" = joint_hn_estims,
            "k_fold_result" = intxn_in_fold,
            "deltas" = deltas_updated
          )
        } else if (sum(stringr::str_count(matches, paste(c(a_names), collapse = "|"))) == 1 &
          sum(stringr::str_count(matches, paste(c(z_names), collapse = "|"))) == 1) {
          ## get the exposure and mediator variables
          exposure <- stringr::str_extract(
            matches,
            paste(c(a_names), collapse = "|")
          )
          mediator <- stringr::str_extract(matches, paste(c(z_names),
            collapse = "|"
          ))
          exposure <- exposure[!is.na(exposure)]
          mediator <- mediator[!is.na(mediator)]

          ## get delta from the list

          delta <- deltas[[exposure]]

          lower_bound <- min(min(av[[exposure]]), min(at[[exposure]]))
          upper_bound <- max(max(av[[exposure]]), max(at[[exposure]]))


          if (exposure_quantized == TRUE) {
            g_learner <- quant_learner
            outcome_type <- "categorical"
          } else {
            g_learner <- pi_learner
            outcome_type <- "continuous"
          }

          if (mediator_quantized == TRUE) {
            med_learner <- quant_learner
            outcome_type <- "categorical"
          } else {
            med_learner <- pi_learner
            outcome_type <- "continuous"
          }

          ## get g(A|W) under shifts and no shift
          gn_exp_estim <- indiv_stoch_shift_est_g_exp(
            exposure = exposure,
            delta = delta,
            g_learner = g_learner,
            covars = w_names,
            av = av,
            at = at,
            adaptive_delta = adaptive_delta,
            hn_trunc_thresh = hn_trunc_thresh,
            exposure_quantized = exposure_quantized,
            lower_bound = lower_bound,
            upper_bound = upper_bound,
            outcome_type = outcome_type,
            density_type = density_type,
            n_bins = n_bins,
            max_degree = max_degree
          )

          ## extract model for future integration
          g_model <- gn_exp_estim$model

          ## if delta adapt is true get the new delta, if not true
          ## delta_updated == delta

          delta_updated <- gn_exp_estim$delta

          ## get g(A|Z,W) under shifts and no shift

          gn_exp_estim_z <- indiv_stoch_shift_est_g_exp(
            exposure = exposure,
            delta = delta,
            g_learner = g_learner,
            covars = c(w_names, mediator),
            av = av,
            at = at,
            adaptive_delta = adaptive_delta,
            hn_trunc_thresh = hn_trunc_thresh,
            exposure_quantized = exposure_quantized,
            lower_bound = lower_bound,
            upper_bound = upper_bound,
            outcome_type = outcome_type,
            density_type = density_type,
            n_bins = n_bins,
            max_degree = max_degree
          )

          e_model <- gn_exp_estim_z$model

          ## get r(Z|W) under shifts and no shift

          if (mediator_quantized == TRUE) {
            zn_exp_estim <- indiv_stoch_shift_est_g_exp(
              exposure = mediator,
              delta = delta,
              g_learner = g_learner,
              covars = w_names,
              av = av,
              at = at,
              adaptive_delta = adaptive_delta,
              hn_trunc_thresh = hn_trunc_thresh,
              exposure_quantized = mediator_quantized,
              lower_bound = lower_bound,
              upper_bound = upper_bound,
              outcome_type = outcome_type,
              density_type = density_type,
              n_bins = n_bins,
              max_degree = max_degree
            )
          }else{

          zn_exp_estim <- joint_stoch_shift_est_z_exp(
            exposures = exposure,
            mediator = mediator,
            deltas = delta_updated,
            pi_learner = pi_learner,
            w_names = w_names,
            a_names = a_names,
            z_names = z_names,
            av = av,
            at = at
          )
          }

          r_model <- zn_exp_estim$model

          ## calc ratios for training and validation used for EIF

          g_shift_e_ratio_av <- (gn_exp_estim$av$upshift / gn_exp_estim_z$av$noshift)
          g_shift_e_ratio_at <- (gn_exp_estim$at$upshift / gn_exp_estim_z$at$noshift)

          ## estimate Q(Y|A,Z,W) under shift and no shift

          qn_estim <- est_Q_w_shifted_mediation(
            exposure = exposure,
            mediator = mediator,
            delta = delta_updated,
            mu_learner = mu_learner,
            covars = c(w_names, exposure, mediator),
            av = av,
            at = at,
            upper_bound = upper_bound,
            lower_bound = lower_bound
          )

          q_model <- qn_estim$model

          qn_estim_no_med <- est_Q_w_shifted_mediation(
            exposure = exposure,
            mediator = mediator,
            delta = delta_updated,
            mu_learner = mu_learner,
            covars = c(w_names, exposure),
            av = av,
            at = at,
            upper_bound = upper_bound,
            lower_bound = lower_bound
          )

          # g_shift_e_ratio <- ifelse(g_shift_e_ratio > hn_trunc_thresh, hn_trunc_thresh, g_shift_e_ratio)

          ## calculate dy part of EIF
          d_y <- g_shift_e_ratio_av * (av$y - qn_estim$av_predictions$noshift)

          if (exposure_quantized == TRUE) {
            d_z_w <- integrate_m_g_quant(
              av = av,
              at = at,
              covars = c(w_names, exposure, mediator),
              w_names = w_names,
              q_model = q_model,
              g_model = g_model,
              exposure = exposure,
              g_delta = delta_updated,
              m_delta = 0,
              n_bins = 4
            )
          } else {
            d_z_w <- integrate_q_g(
              av = av,
              at = at,
              covars = c(w_names, exposure, mediator),
              w_names = w_names,
              q_model = q_model,
              g_model = g_model,
              exposure = exposure,
              g_delta = delta_updated,
              m_delta = 0,
              n_samples = n_mc_sample,
              density_type = density_type,
              lower_bound = lower_bound,
              upper_bound = upper_bound,
              integration_method = integration_method
            )
          }

          if (exposure_quantized == TRUE) {
            d_a_int <- integrate_psi_g_discrete(
              av = av,
              at = at,
              covars = c(w_names, exposure, mediator),
              w_names = w_names,
              q_model = q_model,
              r_model = r_model,
              g_model = g_model,
              exposure = exposure,
              mediator = mediator,
              delta = delta_updated,
              n_samples = n_mc_sample,
              n_bins = 4,
              method = integration_method,
              mediator_quantized = mediator_quantized
            )
          } else {
            d_a_int <- integrate_psi_g(
              av = av,
              at = at,
              covars = c(w_names, exposure, mediator),
              w_names = w_names,
              q_model = q_model,
              r_model = r_model,
              g_model = g_model,
              exposure = exposure,
              mediator = mediator,
              delta = delta_updated,
              n_samples = n_mc_sample,
              n_iterations = 1,
              density_type = density_type,
              integration_method
            )
          }

          ## calculate the g/e no shift ratios needed for pseudo regression

          g_e_shift_ratio_av <- (gn_exp_estim$av$noshift / gn_exp_estim_z$av$noshift)
          g_e_shift_ratio_at <- (gn_exp_estim$at$noshift / gn_exp_estim_z$at$noshift)

          g_e_shift_ratio_av <- ifelse(g_e_shift_ratio_av >= hn_trunc_thresh, hn_trunc_thresh, g_e_shift_ratio_av)
          g_e_shift_ratio_at <- ifelse(g_e_shift_ratio_at >= hn_trunc_thresh, hn_trunc_thresh, g_e_shift_ratio_at)

          ##  pseudo regression to test against d_a
          pseudo_regression_at <- g_e_shift_ratio_at * qn_estim$at_predictions$upshift
          pseudo_regression_av <- g_e_shift_ratio_av * qn_estim$av_predictions$upshift


          at$pseudo_outcome_scaled <- scale_to_unit(pseudo_regression_at)
          av$pseudo_outcome_scaled <- scale_to_unit(pseudo_regression_av)

          at$pseudo_outcome <- pseudo_regression_at
          av$pseudo_outcome <- pseudo_regression_av

          sl_pseudo_task_at <- sl3::sl3_Task$new(
            data = at,
            outcome = "pseudo_outcome_scaled",
            covariates = c(exposure, w_names),
            outcome_type = "quasibinomial"
          )

          sl_pseudo_task_av <- sl3::sl3_Task$new(
            data = av,
            outcome = "pseudo_outcome_scaled",
            covariates = c(exposure, w_names),
            outcome_type = "quasibinomial"
          )

          sl <- sl3::Lrnr_sl$new(
            learners = mu_learner,
            metalearner = sl3::Lrnr_cv_selector$new(sl3::loss_loglik_binomial)
          )

          pseudo_model <- sl$train(sl_pseudo_task_at)
          psi_aw_av <- pseudo_model$predict(sl_pseudo_task_av)
          psi_aw_at <- pseudo_model$predict(sl_pseudo_task_at)

          psi_aw_av <- scale_to_original(psi_aw_av, max_orig = max(av$pseudo_outcome), min_orig = min(av$pseudo_outcome))

          if (exposure_quantized == TRUE) {
            d_a_pseudo <- integrate_psi_aw_g_quant(
              at = at,
              av = av,
              covars = c(exposure, w_names),
              w_names,
              pseudo_model,
              g_model,
              exposure,
              delta = 0,
              psi_aw = psi_aw_av,
              n_bins = 4
            )
          } else {
            d_a_pseudo <- integrate_psi_aw_g(
              at = at,
              av = av,
              covars = c(exposure, w_names),
              w_names,
              pseudo_model,
              g_model,
              exposure,
              psi_aw = psi_aw_av,
              n_samples = n_mc_sample,
              density_type = density_type,
              integration_method = integration_method
            )
          }

          calc_eif_comp_sum <- function(d_y, d_z_w, d_a) {
            return(d_y + d_z_w + d_a)
          }

          calc_psi_shift <- function(eif_comp_sum) {
            return(mean(eif_comp_sum))
          }

          calc_a_shift_fix_z_eif <- function(eif_comp_sum, psi_shift) {
            return(eif_comp_sum - psi_shift)
          }

          calc_psi_nde <- function(psi_shift, av_y) {
            return(psi_shift - mean(av_y))
          }

          calc_eif_nde <- function(a_shift_fix_z_eif, eif_no_shift) {
            return(a_shift_fix_z_eif - eif_no_shift)
          }

          calc_variance_est <- function(eif) {
            return(var(eif) / length(eif))
          }

          calc_se_est <- function(variance_est) {
            return(sqrt(variance_est))
          }

          calc_p_value <- function(psi_shift, se_est) {
            return(2 * stats::pnorm(abs(psi_shift / se_est), lower.tail = F))
          }


          # Sum the nuisance components with each strategy
          eif_comp_sum_w_pseudo <- calc_eif_comp_sum(d_y, d_z_w, d_a_pseudo)
          eif_comp_sum_w_double_int <- calc_eif_comp_sum(d_y, d_z_w, d_a_int$d_a)

          # Subtract off the average to get the EIF for each strategy
          psi_a_shift_fix_z_pseudo <- calc_psi_shift(eif_comp_sum_w_pseudo)
          psi_a_shift_fix_z_double_int <- calc_psi_shift(eif_comp_sum_w_double_int)

          a_shift_fix_z_eif_w_pseudo <- calc_a_shift_fix_z_eif(eif_comp_sum_w_pseudo, psi_a_shift_fix_z_pseudo)
          a_shift_fix_z_eif_w_double_int <- calc_a_shift_fix_z_eif(eif_comp_sum_w_double_int, psi_a_shift_fix_z_double_int)

          eif_no_shift <- av$y - qn_estim$av_predictions$noshift

          psi_nde_pseudo <- calc_psi_nde(psi_a_shift_fix_z_pseudo, av$y)
          psi_nde_double_int <- calc_psi_nde(psi_a_shift_fix_z_double_int, av$y)

          eif_nde_pseudo <- calc_eif_nde(a_shift_fix_z_eif_w_pseudo, eif_no_shift)
          eif_nde_double_int <- calc_eif_nde(a_shift_fix_z_eif_w_double_int, eif_no_shift)

          # Calculate variance, standard error, and confidence intervals
          variance_est_nde_pseudo <- calc_variance_est(eif_nde_pseudo)
          variance_est_nde_double_int <- calc_variance_est(eif_nde_double_int)

          se_est_nde_pseudo <- calc_se_est(variance_est_nde_pseudo)
          se_est_nde_double_int <- calc_se_est(variance_est_nde_double_int)

          CI_nde_pseudo <- calc_CIs(psi_nde_pseudo, se_est_nde_pseudo)
          CI_nde_double_int <- calc_CIs(psi_nde_double_int, se_est_nde_double_int)

          # Calculate p-values
          p_value_nde_pseudo <- calc_p_value(psi_nde_pseudo, se_est_nde_pseudo)
          p_value_nde_double_int <- calc_p_value(psi_nde_double_int, se_est_nde_double_int)

          ind_qn_estim <- indiv_stoch_shift_est_Q(
            exposure = exposure,
            delta = delta_updated,
            mu_learner = mu_learner,
            covars = c(w_names, exposure),
            av = av,
            at = at,
            upper_bound = upper_bound,
            lower_bound = lower_bound
          )

          Hn <- gn_exp_estim$Hn_av

          hn_truncated <- Hn

          hn_truncated$noshift <- ifelse(Hn$noshift > hn_trunc_thresh, hn_trunc_thresh, Hn$noshift)
          hn_truncated$shift <- ifelse(Hn$shift > hn_trunc_thresh, hn_trunc_thresh, Hn$shift)

          tmle_fit <- tmle_exposhift(
            data_internal = av,
            delta = delta_updated,
            Qn_scaled = ind_qn_estim$q_av,
            Qn_unscaled = qn_estim_no_med$av_predictions,
            Hn = Hn,
            fluctuation = fluctuation,
            y = av$y,
            estimator = "onestep"
          )

          tmle_fit$call <- call

          psi_nie_pseudo <- tmle_fit$psi - psi_a_shift_fix_z_pseudo
          psi_nie_double_int <- tmle_fit$psi - psi_a_shift_fix_z_double_int

          eif_nie_pseudo <- (tmle_fit$eif - a_shift_fix_z_eif_w_pseudo)
          eif_nie_double_int <- (tmle_fit$eif - a_shift_fix_z_eif_w_double_int)

          variance_est_nie_pseudo <- var(eif_nie_pseudo) / length(eif_nie_pseudo)
          variance_est_nie_double_int <- var(eif_nie_double_int) / length(eif_nie_double_int)

          se_est_nie_pseudo <- sqrt(variance_est_nie_pseudo)
          se_est_nie_double_int <- sqrt(variance_est_nie_double_int)

          CI_nie_pseudo <- calc_CIs(psi_nie_pseudo, se_est_nie_pseudo)
          CI_nie_double_int <- calc_CIs(psi_nie_double_int, se_est_nie_double_int)

          lower_CI_nie_pseudo <- CI_nie_pseudo[[1]]
          upper_CI_nie_pseudo <- CI_nie_pseudo[[2]]

          lower_CI_nie_double_int <- CI_nie_double_int[[1]]
          upper_CI_nie_double_int <- CI_nie_double_int[[2]]

          p_value_nie_pseudo <- calc_p_value(psi_nie_pseudo, se_est_nie_pseudo)
          p_value_nie_double_int <- calc_p_value(psi_nie_double_int, se_est_nie_double_int)


          nde_results_pseudo <- list(
            "Parameter" = "NDE-Pseudo-Reg",
            "Psi" = psi_nde_pseudo,
            "Variance" = variance_est_nde_pseudo,
            "SE" = se_est_nde_pseudo,
            "Lower CI" = CI_nde_pseudo[1],
            "Upper CI" = CI_nde_pseudo[2],
            "P-Value" = p_value_nde_pseudo
          )

          nde_results_double_int <- list(
            "Parameter" = "NDE-Double-Int",
            "Psi" = psi_nde_double_int,
            "Variance" = variance_est_nde_double_int,
            "SE" = se_est_nde_double_int,
            "Lower CI" = CI_nde_double_int[1],
            "Upper CI" = CI_nde_double_int[2],
            "P-Value" = p_value_nde_double_int
          )

          nie_results_pseudo <- list(
            "Parameter" = "NIE-Pseudo-Reg",
            "Psi" = psi_nie_pseudo,
            "Variance" = variance_est_nie_pseudo,
            "SE" = se_est_nie_pseudo,
            "Lower CI" = CI_nie_pseudo[1],
            "Upper CI" =  CI_nie_pseudo[2],
            "P-Value" = p_value_nie_pseudo
          )

          nie_results_double_int <- list(
            "Parameter" = "NIE-Double-Int",
            "Psi" = psi_nie_double_int,
            "Variance" = variance_est_nie_double_int,
            "SE" = se_est_nie_double_int,
            "Lower CI" = CI_nie_double_int[1],
            "Upper CI" = CI_nie_double_int[2],
            "P-Value" = p_value_nie_double_int
          )

          total_effect <- calc_final_ind_shift_param(
            tmle_fit,
            exposure,
            fold_k
          )

          total_effect <- total_effect[, 1:7]
          colnames(total_effect) <- c(
            "Parameter", "Psi", "Variance", "SE",
            "Lower CI", "Upper CI", "P-Value"
          )

          total_effect$Parameter <- "Total Effect"

          mediation_in_fold <- rbind(
            nde_results_pseudo, nde_results_double_int,
            nie_results_pseudo, nie_results_double_int,
            total_effect
          )

          rownames(mediation_in_fold) <- NULL

          eif_comp_sum_w_pseudo <- d_y + d_z_w + d_a_pseudo
          eif_comp_sum_w_double_int <- d_y + d_z_w + d_a_int$d_a


          fold_results_mediation[[
            paste("Fold", fold_k, ":", paste(exposure, mediator, sep = ""))
          ]] <- list(
            "data" = av,
            "delta" = delta_updated,
            "eif_comp_sum_w_pseudo" = eif_comp_sum_w_pseudo,
            "eif_comp_sum_w_double_int" = eif_comp_sum_w_double_int,
            "Qn_scaled" = ind_qn_estim$q_av,
            "Qn_unscaled" = qn_estim_no_med$av_predictions,
            "Hn" = Hn,
            "eif_no_shift" = eif_no_shift,
            "k_fold_result" = mediation_in_fold,
            "dy" = d_y,
            "dzw" = d_z_w,
            "da_pseudo" = d_a_pseudo,
            "da_int" = d_a_int$d_a
          )
        } else if (sum(stringr::str_count(matches, paste(c(a_names), collapse = "|"))) == 2 &
          sum(stringr::str_count(matches, paste(c(z_names), collapse = "|"))) == 1) {
          exposures <- stringr::str_extract(
            matches,
            paste(c(a_names), collapse = "|")
          )
          exposures <- exposures[!is.na(exposures)]

          exposures <- as.list(exposures)
          deltas_med <- deltas[unlist(exposures)]
          exposures[[3]] <- unlist(exposures)

          mediator <- stringr::str_extract(matches, paste(c(z_names),
            collapse = "|"
          ))

          mediator <- mediator[!is.na(mediator)]

          covars <- c(w_names)

          joint_gn_exp_estims <- joint_stoch_shift_est_g_exp(
            exposures = exposures,
            deltas = deltas_med,
            pi_learner = pi_learner,
            covars = covars,
            av = av,
            at = at,
            adaptive_delta = adaptive_delta,
            hn_trunc_thresh = hn_trunc_thresh
          )

          deltas_updated <- joint_gn_exp_estims$delta_results
          deltas_updated[[3]] <- c(deltas_updated[[1]], deltas_updated[[2]])

          joint_gn_exp_estims$gn_results[[3]] <- as.data.frame(mapply(
            `*`,
            joint_gn_exp_estims$gn_results[[1]],
            joint_gn_exp_estims$gn_results[[3]]
          ))

          covars <- c(w_names, a_names)

          zn_exp_estim <- joint_stoch_shift_est_z_exp(
            exposures = exposures,
            mediator = mediator,
            deltas = deltas_updated,
            pi_learner = pi_learner,
            w_names = w_names,
            a_names = a_names,
            z_names = z_names,
            av = av,
            at = at
          )

          # for the direct effects we just need the density of each A because
          # z is not shifted
          nde_gn_estim_a1 <- est_hn(gn_exp = joint_gn_exp_estims$gn_results[[1]])
          nde_gn_estim_a2 <- est_hn(gn_exp = joint_gn_exp_estims$gn_results[[2]])
          nde_gn_estim_a1a2 <- est_hn(gn_exp = joint_gn_exp_estims$gn_results[[3]])

          # truncate the NDE estimates

          nde_gn_estim_a1_trunc <- as.data.frame(apply(nde_gn_estim_a1, 2, function(x) {
            ifelse(x >= hn_trunc_thresh, hn_trunc_thresh, x)
          }))
          nde_gn_estim_a2_trunc <- as.data.frame(apply(nde_gn_estim_a2, 2, function(x) {
            ifelse(x >= hn_trunc_thresh, hn_trunc_thresh, x)
          }))
          nde_gn_estim_a1a2_trunc <- as.data.frame(apply(nde_gn_estim_a1a2, 2, function(x) {
            ifelse(x >= hn_trunc_thresh, hn_trunc_thresh, x)
          }))


          # for the indirect effects we need the density of each A and the density
          # of Z given a shift in the respective A.

          nie_gn_estim_a1 <- est_hn(gn_exp = as.data.frame(joint_gn_exp_estims$gn_results[[1]]) * as.data.frame(zn_exp_estim$gn_results[[1]]))
          nie_gn_estim_a2 <- est_hn(gn_exp = as.data.frame(joint_gn_exp_estims$gn_results[[2]]) * as.data.frame(zn_exp_estim$gn_results[[2]]))
          nie_gn_estim_a1a2 <- est_hn(gn_exp = as.data.frame(joint_gn_exp_estims$gn_results[[3]]) * as.data.frame(zn_exp_estim$gn_results[[3]]))

          # truncate the NIE estimates

          nie_gn_estim_a1_trunc <- as.data.frame(apply(nie_gn_estim_a1, 2, function(x) {
            ifelse(x >= hn_trunc_thresh, hn_trunc_thresh, x)
          }))
          nie_gn_estim_a2_trunc <- as.data.frame(apply(nie_gn_estim_a2, 2, function(x) {
            ifelse(x >= hn_trunc_thresh, hn_trunc_thresh, x)
          }))
          nie_gn_estim_a1a2_trunc <- as.data.frame(apply(nie_gn_estim_a1a2, 2, function(x) {
            ifelse(x >= hn_trunc_thresh, hn_trunc_thresh, x)
          }))

          nde_hns <- list(nde_gn_estim_a1_trunc, nde_gn_estim_a2_trunc, nde_gn_estim_a1a2_trunc)
          nie_hns <- list(nie_gn_estim_a1_trunc, nie_gn_estim_a2_trunc, nie_gn_estim_a1a2_trunc)

          covars <- c(a_names, w_names)

          zn_estim <- estimate_mediator_joint(
            mediator = mediator,
            exposure = exposures,
            deltas = deltas_updated,
            mu_learner = mu_learner,
            covars = covars,
            av = av,
            at = at
          )

          covars <- c(w_names, a_names, z_names)

          qn_estim <- est_Q_w_shifted_joint_mediation(
            exposures = exposures,
            mediator = mediator,
            delta = deltas_updated,
            mu_learner = mu_learner,
            covars = covars,
            av = av,
            at = at,
            zn_estim = zn_estim
          )

          for (i in seq(exposures)) {
            exposure <- exposures[i]
            delta <- deltas_updated[i]

            tmle_fit_a_shift <- tmle_exposhift(
              data_internal = av,
              delta = delta,
              Qn_scaled = qn_estim[[i]]$a_shifted,
              Hn = nde_hns[[i]],
              fluctuation = fluctuation,
              y = av$y
            )

            tmle_fit_a_z_shift <- tmle_exposhift(
              data_internal = av,
              delta = delta,
              Qn_scaled = qn_estim[[i]]$a_z_shifted,
              Hn = nie_hns[[i]],
              fluctuation = fluctuation,
              y = av$y
            )

            joint_mediation_in_fold <- calc_mediation_param(
              tmle_fit_a_shift = tmle_fit_a_shift,
              tmle_fit_a_z_shift = tmle_fit_a_z_shift,
              exposure,
              mediator,
              y = av$y,
              fold_k = fold_k,
              delta = delta
            )

            joint_fold_results_mediation[[
              paste("Fold", fold_k, ":", paste(paste(unlist(exposure), collapse = ""), mediator, sep = ""))
            ]] <- list(
              "data" = av,
              "Qn_a_shift_scaled" = qn_estim[[i]]$a_shifted,
              "Qn_a_z_shift_scaled" = qn_estim[[i]]$a_z_shifted,
              "Hn_a_shift" = nde_hns[[i]],
              "Hn_az_shift" = nie_hns[[i]],
              "k_fold_result" = joint_mediation_in_fold,
              "delta" = unlist(delta)
            )
          }
        }
      }

      results_list <- list(
        fold_results_indiv,
        fold_results_em,
        fold_results_intxn,
        fold_results_mediation,
        joint_fold_results_mediation
      )

      names(results_list) <- c(
        "indiv_shift",
        "em_shift",
        "intxn_shift",
        "med_shift",
        "joint_med_shift"
      )

      results_list
    },
    .options = furrr::furrr_options(seed = seed, packages = "SuperNOVA")
  )

  indiv_shift_results <- purrr::map(fold_SuperNOVA_results, c("indiv_shift"))
  em_shift_results <- purrr::map(fold_SuperNOVA_results, c("em_shift"))
  intxn_shift_results <- purrr::map(fold_SuperNOVA_results, c("intxn_shift"))
  med_shift_results <- purrr::map(fold_SuperNOVA_results, c("med_shift"))
  joint_med_shift_results <- purrr::map(fold_SuperNOVA_results, c("joint_med_shift"))

  indiv_shift_results <- unlist(indiv_shift_results, recursive = FALSE)
  em_shift_results <- unlist(em_shift_results, recursive = FALSE)
  intxn_shift_results <- unlist(intxn_shift_results, recursive = FALSE)
  med_shift_results <- unlist(med_shift_results, recursive = FALSE)
  joint_med_shift_results <- unlist(joint_med_shift_results, recursive = FALSE)


  if (!is.null(indiv_shift_results)) {
    pooled_indiv_shift_results <- calc_pooled_indiv_shifts(
      indiv_shift_results = indiv_shift_results,
      estimator = estimator,
      fluctuation = fluctuation
    )
  } else {
    indiv_shift_results <- NULL
  }

  if (!is.null(em_shift_results)) {
    pooled_em_shift_results <- calc_pooled_em_shifts(
      y,
      em_shift_results = em_shift_results,
      estimator = estimator,
      w_names = w_names,
      a_names = a_names,
      fluctuation = fluctuation,
      em_learner = em_learner
    )
  } else {
    pooled_em_shift_results <- NULL
  }

  if (!is.null(intxn_shift_results)) {
    pooled_intxn_shift_results <- calc_pooled_intxn_shifts(
      intxn_shift_results = intxn_shift_results,
      estimator = estimator,
      a_names = a_names,
      w_names = w_names,
      z_names = z_names,
      fluctuation = fluctuation
    )
  } else {
    pooled_intxn_shift_results <- NULL
  }

  if (!is.null(med_shift_results)) {
    pooled_med_shift_results <- calc_pooled_med_shifts(
      med_shift_results = med_shift_results,
      estimator = estimator,
      fluctuation = fluctuation,
      a_names = a_names,
      w_names = w_names,
      z_names = z_names
    )
  } else {
    pooled_med_shift_results <- NULL
  }


  results_list <- list(
    "Basis Fold Proportions" = basis_prop_in_fold,
    "Effect Mod Results" = pooled_em_shift_results,
    "Indiv Shift Results" = pooled_indiv_shift_results,
    "Joint Shift Results" = pooled_intxn_shift_results,
    "Mediation Shift Results" = pooled_med_shift_results
  )

  return(results_list)
}
