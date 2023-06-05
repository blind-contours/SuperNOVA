#' @title Data-adaptive estimation of interactions, effect modification, and
#' mediation using stochastic shift intervention target parameters. In many mixed exposure settings,
#' interactions in the mixture, effect modifiers in the covariates that modify the
#' impact of an exposure and mediating pathways from exposure to outcome are generally unknown. SuperNOVA finds these variable sets
#' on one part of the data and estimates counterfactual outcome changes given shifts to exposure on an estimation part of the data.
#' Using cross-validation and targeted learning, estimators are created that utlize machine learning that are unbiased and have the
#' minimum variance.
#'
#' @description The SuperNOVA function provides an efficient approach to estimate
#' interactions, effect modification, and mediation using targeted minimum loss
#' estimators for counterfactual mean differences under various target parameters.
#' The procedure employs data-adaptive ensemble b-spline models and stochastic interventions,
#' leveraging the \pkg{sl3} package for ensemble machine learning. The data is split into V folds, in each fold
#' the training data is used to find variable sets using flexible basis function estimators. Given the different variable sets,
#' stochastic intervention target parameters are applied with cross-validated targeted learning.
#'
#' @param w A \code{matrix}, \code{data.frame}, or similar containing a set of
#' baseline covariates. These variables are measured before exposures.
#' @param a \code{matrix}, \code{data.frame}, or similar containing individual or
#' multiple exposures.
#' @param z \code{matrix}, \code{data.frame}, or similar containing individual or
#' multiple mediators (optional).
#' @param y \code{numeric} vector of observed outcomes.
#' @param deltas A \code{numeric} value indicating the shift in exposures to
#' define the target parameter, with respect to the scale of the exposures (A). If adaptive_delta
#' is true, these values will be reduced.
#' @param var_sets A list specifying variable sets for deterministic SuperNOVA usage.
#' Example: var_sets <- c("A_1", "A_1-Z_2") where the analyst provides variable sets
#' for exposures, exposure-mediator, or exposure-covariate relationships.
#' @param estimator The type of estimator to fit: \code{"tmle"} for targeted
#' maximum likelihood estimation, or \code{"onestep"} for a one-step estimator.
#' @param fluctuation Method used in the targeting step for TML estimation: "standard" or "weighted".
#' This determines where to place the auxiliary covariate in the logistic tilting regression.
#' @param pi_learner Learners for fitting Super Learner ensembles to densities via \pkg{sl3}.
#' @param mu_learner Learners for fitting Super Learner ensembles to the outcome model via \pkg{sl3}.
#' @param g_learner Learners for fitting Super Learner ensembles to the g-mechanism
#' g(A|W) (a probability estimator, not a density estimator) for mediation via \pkg{sl3}.
#' @param e_learner Learners for fitting Super Learner ensembles to the e-mechanism
#' g(A|Z,W) (a probability estimator, not a density estimator) for mediation via \pkg{sl3}.
#' @param zeta_learner Learners for fitting Super Learner ensembles to the outcome model via \pkg{sl3}..
#' @param n_folds Number of folds to use in cross-validation, default is 2.
#' @param outcome_type Data type of the outcome, default is "continuous".
#' @param quantile_thresh Threshold based on quantiles of the F-statistic, used to
#' identify "important" basis functions in the data-adaptive procedure.
#' @param verbose Whether to run verbosely (default: FALSE).
#' @param parallel Whether to parallelize across cores (default: TRUE).
#' @param parallel_type Type of parallelization to use if parallel is TRUE:
#' "multi_session" (default), "multicore", or "sequential".
#' @param num_cores Number of CPU cores to use in parallelization (default: 2).
#' @param seed \code{numeric} seed value to be passed to all functions.
#' @param hn_trunc_thresh Truncation level for the clever covariate (default: 10).
#' @param adaptive_delta If TRUE, reduces the user-specified delta until
#' the Hn calculated for a shift does not have any observation greater
#' than hn_trunc_thresh (default: FALSE).
#' @param n_mc_sample Number of iterations to be used for the Monte Carlo integration
#' procedure when using continuous exposures (default: 1000).
#' @param exposure_quantized Whether the exposure has been discretized into bins,
#' in which case the integration procedure is skipped and weighted sums are used instead (default: FALSE).
#' @param mediator_quantized If the mediator is discretized, a multinomial ML function
#' is used in this regression to avoid density estimation (default: FALSE).
#' @param density_type Type of density estimation to be used: "sl" for Super Learner
#' (default) or "hal" for highly adaptive lasso.
#' @param n_bins Number of bins for quantizing the exposure if mediation is detected (default: 10).
#' @param max_degree Maximum degree of interactions used in the highly adaptive lasso
#' density estimator if used (default: 1).
#' @param integration_method Type of integration to be used in the continuous exposure
#' case: "MC" for Monte Carlo integration (default) or "AQ" for adaptive quadrature.
#' @param use_multinomial Whether to use multinomial regression for binned exposures
#' (default: FALSE).
#'
#' @return An S3 object of class \code{SuperNOVA} containing the results of the
#' procedure to compute a TML or one-step estimate of the counterfactual mean
#' under a modified treatment policy that shifts a continuous-valued exposure
#' by a scalar amount \code{delta}. These exposures are data-adaptively
#' identified using the CV-TMLE procedure.
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
                      zeta_learner = NULL,
                      n_folds = 2,
                      outcome_type = "continuous",
                      mediator_type = "continuous",
                      quantile_thresh = 0,
                      verbose = FALSE,
                      parallel = TRUE,
                      parallel_type = "multi_session",
                      num_cores = 2,
                      seed = seed,
                      hn_trunc_thresh = 10,
                      adaptive_delta = FALSE,
                      n_mc_sample = 1000,
                      exposure_quantized = TRUE,
                      mediator_quantized = FALSE,
                      density_type = "sl",
                      n_bins = 10,
                      max_degree = 2,
                      integration_method = "MC",
                      use_multinomial = FALSE,
                      discover_only = FALSE) {
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

  if (exposure_quantized == TRUE) {
    sls <- create_sls()
    quant_learner <- sls$quant_learner
  }

  if (exposure_quantized == TRUE) {
    if(use_multinomial == TRUE){
      exp_learner <- quant_learner
      exposure_type <- "categorical"
    }else{
      exp_learner <- pi_learner
      exposure_type <- "continuous"
    }

  } else {
    exp_learner <- pi_learner
    exposure_type <- "continuous"
  }

  if (mediator_quantized == TRUE) {
    med_learner <- quant_learner
    med_type <- "categorical"
  } else {
    med_learner <- pi_learner
    med_type <- "continuous"
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
  } else {
    data_internal$folds <- create_cv_folds(n_folds, data_internal$y)
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
          outcome_type = outcome_type,
          mediator_type = mediator_type,
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

  if (discover_only == TRUE) {
    return(basis_prop_in_fold)
  }

  fold_results_indiv <- list()
  fold_results_em <- list()
  fold_results_intxn <- list()
  fold_results_mediation <- list()
  joint_fold_results_mediation <- list()

  fold_SuperNOVA_results <- furrr::future_map(
    unique(data_internal$folds), function(fold_k) {

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

          at <- data_internal[data_internal$folds != fold_k, ]
          av <- data_internal[data_internal$folds == fold_k, ]

          exposure <- target

          delta <- deltas[[exposure]]

          lower_bound <- min(min(av[[exposure]]), min(at[[exposure]]))
          upper_bound <- max(max(av[[exposure]]), max(at[[exposure]]))


          ind_gn_exp_estim <- indiv_stoch_shift_est_g_exp(
            exposure = target,
            delta = delta,
            g_learner = pi_learner,
            covars = w_names,
            av = av,
            at = at,
            adaptive_delta = adaptive_delta,
            hn_trunc_thresh = hn_trunc_thresh,
            use_multinomial = use_multinomial,
            lower_bound = lower_bound,
            upper_bound = upper_bound,
            outcome_type = exposure_type,
            density_type = density_type,
            n_bins = n_bins,
            max_degree = max_degree
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
            Qn_unscaled = scale_to_original(ind_qn_estim$q_av, min_orig = min(av$y), max_orig = max(av$y)),
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

          at <- data_internal[data_internal$folds != fold_k, ]
          av <- data_internal[data_internal$folds == fold_k, ]

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

          lower_bound <- min(min(av[[exposure]]), min(at[[exposure]]))
          upper_bound <- max(max(av[[exposure]]), max(at[[exposure]]))

          gn_exp_estim <- indiv_stoch_shift_est_g_exp(
            exposure = exposure,
            delta = delta,
            g_learner = pi_learner,
            covars = covars,
            av = av,
            at = at,
            adaptive_delta = adaptive_delta,
            hn_trunc_thresh = hn_trunc_thresh,
            use_multinomial = use_multinomial,
            lower_bound = lower_bound,
            upper_bound = upper_bound,
            outcome_type = exposure_type,
            density_type = density_type,
            n_bins = n_bins,
            max_degree = max_degree
          )

          delta <- gn_exp_estim$delta

          covars <- c(a_names, w_names, z_names)

          qn_estim <- indiv_stoch_shift_est_Q(
            exposure = exposure,
            delta = delta,
            mu_learner = mu_learner,
            covars = covars,
            av = av,
            at = at,
            lower_bound = lower_bound,
            upper_bound = upper_bound
          )

          Hn_av <- gn_exp_estim$Hn_av
          Hn_at <- gn_exp_estim$Hn_at

          tmle_fit_av <- tmle_exposhift(
            data_internal = av,
            delta = delta,
            Qn_scaled = qn_estim$q_av,
            Qn_unscaled = scale_to_original(ind_qn_estim$q_av, min_orig = min(av$y), max_orig = max(av$y)),
            Hn = Hn_av,
            fluctuation = fluctuation,
            y = av$y
          )

          tmle_fit_at <- tmle_exposhift(
            data_internal = at,
            delta = delta,
            Qn_scaled = qn_estim$q_at,
            Qn_unscaled = scale_to_original(ind_qn_estim$q_at, min_orig = min(at$y), max_orig = max(at$y)),
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
            fold_k = fold_k
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

          at <- data_internal[data_internal$folds != fold_k, ]
          av <- data_internal[data_internal$folds == fold_k, ]


          exposures <- as.list(matches)
          delta <- deltas[unlist(exposures)]
          exposures[[3]] <- matches

          covars <- c(w_names)

          joint_gn_exp_estims <- joint_stoch_shift_est_g_exp(
            exposures,
            deltas,
            g_learner = pi_learner,
            covars = covars,
            av = av,
            at = at,
            adaptive_delta = adaptive_delta,
            hn_trunc_thresh = hn_trunc_thresh,
            use_multinomial = use_multinomial,
            density_type = density_type,
            max_degree = max_degree,
            n_bins = n_bins,
            outcome_type = exposure_type
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
            at
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

          at <- data_internal[data_internal$folds != fold_k, ]
          av <- data_internal[data_internal$folds == fold_k, ]

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


          # num_unique <- length(unique(data_internal[[exposure]]))
          # num_total <- length(data_internal[[exposure]])
          # exposure_numeric <- is.numeric(data_internal[[exposure]]) && !is.factor(data_internal[[exposure]]) && num_unique / num_total > 0.5

          # if (exposure_numeric == TRUE & exposure_quantized == FALSE) {
          #   data_internal_c <- data_internal
          #   quantile_breaks <- quantile(data_internal_c[[exposure]], probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE)
          #   data_internal_c[[exposure]] <- as.numeric(cut(data_internal_c[[exposure]], breaks = quantile_breaks, labels = FALSE, include.lowest = TRUE))
          #
          #   at <- data_internal_c[data_internal_c$folds != fold_k, ]
          #   av <- data_internal_c[data_internal_c$folds == fold_k, ]
          # }else{

          # }

          delta <- deltas[[exposure]]

          lower_bound <- min(min(av[[exposure]]), min(at[[exposure]]))
          upper_bound <- max(max(av[[exposure]]), max(at[[exposure]]))


          ## get g(A|W) under shifts and no shift
          gn_exp_estim <- indiv_stoch_shift_est_g_exp(
            exposure = exposure,
            delta = delta,
            g_learner = exp_learner,
            covars = w_names,
            av = av,
            at = at,
            adaptive_delta = adaptive_delta,
            hn_trunc_thresh = hn_trunc_thresh,
            use_multinomial = use_multinomial,
            lower_bound = lower_bound,
            upper_bound = upper_bound,
            outcome_type = exposure_type,
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
            g_learner = exp_learner,
            covars = c(w_names, mediator),
            av = av,
            at = at,
            adaptive_delta = adaptive_delta,
            hn_trunc_thresh = hn_trunc_thresh,
            use_multinomial = use_multinomial,
            lower_bound = lower_bound,
            upper_bound = upper_bound,
            outcome_type = exposure_type,
            density_type = density_type,
            n_bins = n_bins,
            max_degree = max_degree
          )

          e_model <- gn_exp_estim_z$model

          ## get r(Z|W) under shifts and no shift

          zn_exp_estim <- indiv_stoch_shift_est_g_exp(
            exposure = mediator,
            delta = delta,
            g_learner = med_learner,
            covars = w_names,
            av = av,
            at = at,
            adaptive_delta = adaptive_delta,
            hn_trunc_thresh = hn_trunc_thresh,
            use_multinomial = use_multinomial,
            lower_bound = min(av[[mediator]]),
            upper_bound = max(av[[mediator]]),
            outcome_type = med_type,
            density_type = density_type,
            n_bins = n_bins,
            max_degree = max_degree
          )


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
            covars = c(w_names, a_names, mediator),
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


          ## calculate dzw part of EIF
          if (exposure_quantized == TRUE & use_multinomial == TRUE) {
              d_z_w <- integrate_q_g_quant(
                av = av,
                at = at,
                covars = c(w_names, a_names, mediator),
                w_names = w_names,
                q_model = q_model,
                g_model = g_model,
                exposure = exposure,
                g_delta = 0,
                m_delta = delta_updated,
                n_bins = n_bins,
                density_type =  density_type,
                upper_bound = upper_bound
              )
          }

          if (exposure_quantized == TRUE & use_multinomial == FALSE) {
            d_z_w <- integrate_q_g_quant(
              av = av,
              at = at,
              covars = c(w_names, a_names, mediator),
              w_names = w_names,
              q_model = q_model,
              g_model = g_model,
              exposure = exposure,
              g_delta = 0,
              m_delta = delta_updated,
              n_bins = n_bins,
              density_type =  density_type,
              upper_bound = upper_bound
            )
            }
          if (exposure_quantized == FALSE) {
            d_z_w <- integrate_q_g(
              av = av,
              at = at,
              covars = c(w_names, a_names, mediator),
              w_names = w_names,
              q_model = q_model,
              g_model = g_model,
              exposure = exposure,
              g_delta = 0,
              m_delta = delta_updated,
              n_samples = n_mc_sample,
              density_type = density_type,
              lower_bound = lower_bound,
              upper_bound = upper_bound,
              integration_method = integration_method
            )
          }

          ## calculate da part of EIF using double integration

          if (exposure_quantized == TRUE & use_multinomial == TRUE) {
              d_a_int <- integrate_psi_g_discrete(
                av = av,
                at = at,
                covars = c(w_names, a_names, mediator),
                w_names = w_names,
                q_model = q_model,
                r_model = r_model,
                g_model = g_model,
                exposure = exposure,
                mediator = mediator,
                delta = delta_updated,
                n_samples = n_mc_sample,
                n_bins = n_bins,
                method = integration_method,
                mediator_quantized = mediator_quantized,
                density_type = density_type,
                upper_bound = upper_bound,
                use_multinomial = use_multinomial
              )
            }

          if (exposure_quantized == TRUE & use_multinomial == FALSE) {
              d_a_int <- integrate_psi_g_discrete(
                av = av,
                at = at,
                covars = c(w_names, a_names, mediator),
                w_names = w_names,
                q_model = q_model,
                r_model = r_model,
                g_model = g_model,
                exposure = exposure,
                mediator = mediator,
                delta = delta_updated,
                n_samples = n_mc_sample,
                n_bins = n_bins,
                method = integration_method,
                mediator_quantized = mediator_quantized,
                density_type = density_type,
                upper_bound = upper_bound,
                use_multinomial = use_multinomial
              )
            }

          if (exposure_quantized == FALSE){
            d_a_int <- integrate_psi_g_cont(
              av = av,
              at = at,
              covars = c(w_names, a_names, mediator),
              w_names = w_names,
              q_model = q_model,
              r_model = r_model,
              g_model = g_model,
              exposure = exposure,
              mediator = mediator,
              delta = delta_updated,
              n_samples = n_mc_sample,
              density_type = density_type,
              integration_method = integration_method
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


          # at$pseudo_outcome_scaled <- scale_to_unit(pseudo_regression_at)
          # av$pseudo_outcome_scaled <- scale_to_unit(pseudo_regression_av)

          at$pseudo_outcome <- pseudo_regression_at
          av$pseudo_outcome <- pseudo_regression_av

          sl_pseudo_task_at <- sl3::sl3_Task$new(
            data = at,
            outcome = "pseudo_outcome",
            covariates = c(a_names, w_names),
            outcome_type = "continuous"
          )

          sl_pseudo_task_av <- sl3::sl3_Task$new(
            data = av,
            outcome = "pseudo_outcome",
            covariates = c(a_names, w_names),
            outcome_type = "continuous"
          )

          sl <- sl3::Lrnr_sl$new(
            learners = mu_learner,
            metalearner = sl3::Lrnr_cv_selector$new(sl3::loss_loglik_binomial)
          )

          pseudo_model <- sl$train(sl_pseudo_task_at)
          psi_aw_av <- pseudo_model$predict(sl_pseudo_task_av)
          psi_aw_at <- pseudo_model$predict(sl_pseudo_task_at)

          # psi_aw_av <- scale_to_original(psi_aw_av, max_orig = max(av$pseudo_outcome), min_orig = min(av$pseudo_outcome))

          if (exposure_quantized == TRUE & use_multinomial == TRUE) {
            d_a_pseudo <- integrate_psi_aw_g_quant(
              at = at,
              av = av,
              covars = c(a_names, w_names),
              w_names,
              pseudo_model,
              g_model,
              exposure,
              psi_aw = psi_aw_av,
              n_bins = n_bins,
              use_multinomial = TRUE,
              density_type = density_type)
            }

          if (exposure_quantized == TRUE & use_multinomial == FALSE) {
              d_a_pseudo <- integrate_psi_aw_g_quant(
                at = at,
                av = av,
                covars = c(a_names, w_names),
                w_names,
                pseudo_model,
                g_model,
                exposure,
                psi_aw = psi_aw_av,
                use_multinomial = FALSE,
                density_type = density_type)
          }

          if (exposure_quantized == FALSE) {

            d_a_pseudo <- integrate_psi_aw_g(
              at = at,
              av = av,
              covars = c(a_names, w_names),
              w_names = w_names,
              pseudo_model,
              g_model,
              exposure,
              delta = 0,
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
            covars = c(w_names, a_names, z_names[z_names != mediator]),
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
            "Upper CI" = CI_nie_pseudo[2],
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
      fluctuation = fluctuation
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
