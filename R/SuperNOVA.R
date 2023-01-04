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
#' @param a \code{matrix}, \code{data.frame}, or similar containing an
#' individual or set of mediators
#' @param y \code{numeric} vector of the observed outcomes.
#' @param deltas A \code{numeric} value indicating the shift in the exposures to
#'  be used in defining the target parameter. This is defined with respect to
#'  the scale of the exposures (A).
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
#' @param n_folds Number of folds to use in cross-validation
#' @param family Outcome type family
#' @param quantile_thresh Threshold based on quantiles of the f-statistic used
#' to identify "important" basis functions in the data-adaptive procedure
#' @param verbose Whether to run verbosely
#' @param parallel TRUE/FALSE parallelize across cores
#' @param seed \code{numeric} seed value to be passed to all functions
#' @param hn_trunc_thresh Truncation level for the clever covariate
#'
#' @export
#' @importFrom MASS mvrnorm
#' @importFrom bindata rmvbin
#' @importFrom foreach %dopar%
#' @importFrom magrittr %>%
#' @importFrom stats as.formula glm p.adjust plogis predict qlogis qnorm qunif rnorm runif
#' @importFrom rlang :=
#' @importFrom stringr str_count
#' @importFrom data.table rbindlist
#' @return S3 object of class \code{SuperNOVA} containing the results of the
#'  procedure to compute a TML or one-step estimate of the counterfactual mean
#'  under a modified treatment policy that shifts a continuous-valued exposure
#'  by a scalar amount \code{delta}. These exposure are data-adaptively
#'  identified using the CV-TMLE procedure.
#' @examples
#' seed <- 429153
#' n_obs <- 300
#' W <- replicate(2, rbinom(n_obs, 1, 0.5))
#' A <- rnorm(n_obs, mean = 2 * W, sd = 1)
#' Y <- rbinom(n_obs, 1, plogis(A + W + rnorm(n_obs, mean = 0, sd = 1)))

#' example_results <- SuperNOVA(
#'  w = W,
#'  a = A,
#'  z = Z,
#'  y = Y,
#'  delta = 0.5,
#'  estimator = "tmle",
#'  seed = seed,
#'  n_folds = 2
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
                      n_folds,
                      family = "continuous",
                      quantile_thresh = 0.5,
                      verbose = FALSE,
                      parallel = TRUE,
                      parallel_type = "multi_session",
                      num_cores = 2,
                      seed = seed,
                      hn_trunc_thresh = 50,
                      alpha = 0.0001,
                      adaptive_delta = FALSE) {
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
  for (i in 1:dim(a)[2]) {
    check_a <- a[i]
    a_name <- a_names[i]
    if (length(unique(check_a)) < 20) {
      a[, a_name] <- a[, a_name] + runif(dim(a)[1], min = -0.05, max = 0.05)
    }
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

  if (is.null(em_learner)) {
    sls <- create_sls()
    em_learner <- sls$em_learner
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

  if (family == "binomial") {
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

        covars <- c(a_names, w_names, z_names)

        basis_results <- fit_basis_estimators(
          at = at,
          covars = covars,
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

          ind_gn_exp_estim <- indiv_stoch_shift_est_g_exp(
            exposure = target,
            delta = delta,
            pi_learner = pi_learner,
            covars = w_names,
            av,
            at,
            alpha = alpha,
            adaptive_delta = adaptive_delta
          )

          delta <- ind_gn_exp_estim$delta

          covars <- c(a_names, w_names)

          ind_qn_estim <- indiv_stoch_shift_est_Q(
            exposure = exposure,
            delta = delta,
            mu_learner = mu_learner,
            covars = covars,
            av = av,
            at = at
          )

          Hn <- ind_gn_exp_estim$Hn_av
          Hn_trunc <- as.data.frame(apply(Hn, 2, function(x) {
            ifelse(x >= hn_trunc_thresh, hn_trunc_thresh, x)
          }))

          tmle_fit <- tmle_exposhift(
            data_internal = av,
            delta = delta,
            Qn_scaled = ind_qn_estim$q_av,
            Hn = Hn_trunc,
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
            "Hn" = Hn_trunc,
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
            pi_learner = pi_learner,
            covars = covars,
            av = av,
            at = at,
            alpha = alpha,
            adaptive_delta = adaptive_delta
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
          Hn_av_trunc <- as.data.frame(apply(Hn_av, 2, function(x) {
            ifelse(x >= hn_trunc_thresh, hn_trunc_thresh, x)
          }))

          Hn_at <- gn_exp_estim$Hn_at
          Hn_at_trunc <- as.data.frame(apply(Hn_at, 2, function(x) {
            ifelse(x >= hn_trunc_thresh, hn_trunc_thresh, x)
          }))

          tmle_fit_av <- tmle_exposhift(
            data_internal = av,
            delta = delta,
            Qn_scaled = qn_estim$q_av,
            Hn = Hn_av_trunc,
            fluctuation = fluctuation,
            y = av$y
          )

          tmle_fit_at <- tmle_exposhift(
            data_internal = at,
            delta = delta,
            Qn_scaled = qn_estim$q_at,
            Hn = Hn_at_trunc,
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
            "Hn_av" = Hn_av_trunc,
            "Hn_at" = Hn_at_trunc,
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
            pi_learner,
            covars,
            av,
            at,
            delta_adapt = TRUE,
            alpha = alpha,
            adaptive_delta = adaptive_delta
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
            hn_estim_trunc <- as.data.frame(apply(hn_estim, 2, function(x) {
              ifelse(x >= hn_trunc_thresh, hn_trunc_thresh, x)
            }))

            qn_estim <- joint_qn_estims[[i]]
            delta <- deltas_updated[[i]]

            tmle_fit <- tmle_exposhift(
              data_internal = av,
              delta = delta,
              Qn_scaled = qn_estim,
              Hn = hn_estim_trunc,
              fluctuation = fluctuation,
              y = av$y
            )

            intxn_results_list[[i]] <- tmle_fit
            joint_hn_estims[[i]] <- hn_estim_trunc
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
          sum(stringr::str_count(matches, paste(c(z_names), collapse = "|"))) <= 2) {
          exposure <- stringr::str_extract(
            matches,
            paste(c(a_names), collapse = "|")
          )
          mediator <- stringr::str_extract(matches, paste(c(z_names),
            collapse = "|"
          ))
          exposure <- exposure[!is.na(exposure)]
          mediator <- mediator[!is.na(mediator)]

          delta <- deltas[[exposure]]
          covars <- c(w_names)

          gn_exp_estim <- indiv_stoch_shift_est_g_exp(
            exposure = exposure,
            delta = delta,
            pi_learner = pi_learner,
            covars = covars,
            av = av,
            at = at,
            alpha = alpha,
            adaptive_delta = adaptive_delta
          )

          delta_updated <- gn_exp_estim$delta

          covars <- c(w_names, a_names)

          zn_exp_estim <- indiv_stoch_shift_est_z_exp(
            exposure = exposure,
            mediator = mediator,
            delta = delta_updated,
            pi_learner = pi_learner,
            covars = covars,
            av = av,
            at = at
          )

          nde_gn_estim_av <- est_hn(gn_exp = gn_exp_estim$av)
          nie_gn_estim_av <- est_hn(gn_exp = gn_exp_estim$av * zn_exp_estim$av)

          nde_gn_estim_av_trunc <- as.data.frame(apply(nde_gn_estim_av, 2, function(x) {
            ifelse(x >= hn_trunc_thresh, hn_trunc_thresh, x)
          }))
          nie_gn_estim_av_trunc <- as.data.frame(apply(nie_gn_estim_av, 2, function(x) {
            ifelse(x >= hn_trunc_thresh, hn_trunc_thresh, x)
          }))

          covars <- c(a_names, w_names)

          zn_estim <- estimate_mediator(
            mediator = mediator,
            exposure = exposure,
            delta = delta_updated,
            mu_learner = mu_learner,
            covars = covars,
            av = av,
            at = at
          )

          covars <- c(w_names, a_names, z_names)

          qn_estim <- est_Q_w_shifted_mediation(
            exposure = exposure,
            mediator = mediator,
            delta = delta_updated,
            mu_learner = mu_learner,
            covars = covars,
            av = av,
            at = at,
            zn_estim = zn_estim$q_av
          )

          tmle_fit_a_shift <- tmle_exposhift(
            data_internal = av,
            delta = delta_updated,
            Qn_scaled = qn_estim$a_shifted,
            Hn = nde_gn_estim_av_trunc,
            fluctuation = fluctuation,
            y = av$y
          )

          tmle_fit_a_z_shift <- tmle_exposhift(
            data_internal = av,
            delta = delta_updated,
            Qn_scaled = qn_estim$a_z_shifted,
            Hn = nie_gn_estim_av_trunc,
            fluctuation = fluctuation,
            y = av$y
          )

          mediation_in_fold <- calc_mediation_param(
            tmle_fit_a_shift = tmle_fit_a_shift,
            tmle_fit_a_z_shift = tmle_fit_a_z_shift,
            exposure,
            mediator,
            y = av$y,
            fold_k = fold_k,
            delta = delta_updated
          )

          fold_results_mediation[[
            paste("Fold", fold_k, ":", paste(exposure, mediator, sep = ""))
          ]] <- list(
            "data" = av,
            "Qn_a_shift_scaled" = qn_estim$a_shifted,
            "Qn_a_z_shift_scaled" = qn_estim$a_z_shifted,
            "Hn_a_shift" = nde_gn_estim_av_trunc,
            "Hn_az_shift" = nie_gn_estim_av_trunc,
            "k_fold_result" = mediation_in_fold,
            "delta" = delta_updated
          )
        }
      }

      results_list <- list(
        fold_results_indiv,
        fold_results_em,
        fold_results_intxn,
        fold_results_mediation
      )

      names(results_list) <- c(
        "indiv_shift",
        "em_shift",
        "intxn_shift",
        "med_shift"
      )

      results_list
    },
    .options = furrr::furrr_options(seed = seed, packages = "SuperNOVA")
  )

  indiv_shift_results <- purrr::map(fold_SuperNOVA_results, c("indiv_shift"))
  em_shift_results <- purrr::map(fold_SuperNOVA_results, c("em_shift"))
  intxn_shift_results <- purrr::map(fold_SuperNOVA_results, c("intxn_shift"))
  med_shift_results <- purrr::map(fold_SuperNOVA_results, c("med_shift"))

  indiv_shift_results <- unlist(indiv_shift_results, recursive = FALSE)
  em_shift_results <- unlist(em_shift_results, recursive = FALSE)
  intxn_shift_results <- unlist(intxn_shift_results, recursive = FALSE)
  med_shift_results <- unlist(med_shift_results, recursive = FALSE)


  if (!is.null(indiv_shift_results)) {
    pooled_indiv_shift_results <- calc_pooled_indiv_shifts(
      basis_prop_in_fold = basis_prop_in_fold,
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
      basis_prop_in_fold = basis_prop_in_fold,
      em_shift_results = em_shift_results,
      estimator = estimator,
      w_names = w_names,
      a_names = a_names,
      y_name = y_name,
      fluctuation = fluctuation,
      em_learner = em_learner
    )
  } else {
    pooled_em_shift_results <- NULL
  }

  if (!is.null(intxn_shift_results)) {
    pooled_intxn_shift_results <- calc_pooled_intxn_shifts(
      basis_prop_in_fold = basis_prop_in_fold,
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
      basis_prop_in_fold = basis_prop_in_fold,
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
