#' @title Fit the zeta learner, a highly flexible estimator using splines
#' or basis functions in order to identify variable sets for the shift
#' interventions.
#'
#' @details Finds the best fitting flexible estimator using a discrete
#' Super Learner. This function non-parametrically does variable importance
#' of variable sets, such as identificaiton of two mixture components or
#' mixture components and baseline covariates which explain the outcome given
#' an ANOVA like decomposition. We use an F-statistic to threshold the
#' importance of variable sets and treat these variable sets as the data-
#' adaptive parameter.
#'
#' @param at Training dataframe
#' @param covars Covariates to be used as predictors in the Super Learner
#' @param outcome Variable name for the outcome
#' @param family Family of outcome variable
#' @param quantile_thresh Quantile level to set the f-statistic in determining
#' basis functions for estimation
#' @param zeta_learner Stack of algorithms made in SL 3 used in ensemble machine
#' learning to fit Y|A,W
#' @param fold Current fold in the cross-validation
#' @param seed Seed number for consistent results
#' @import sl3
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom dplyr group_by filter top_n
#' @importFrom stats lm anova model.matrix quantile
#' @importFrom purrr is_empty
#' @import polspline
#' @importFrom stringr str_extract_all
#' @return A list of results from fitting the best b-spline model to the data
#' this list includes the selected learner model, the name of the learner,
#' the anova fit on the model matrix of the basis functions, the basis used,
#' and the residual metrics

#' @export

fit_basis_estimators <- function(at,
                                 a_names,
                                 z_names,
                                 w_names,
                                 outcome,
                                 family,
                                 quantile_thresh,
                                 zeta_learner,
                                 fold,
                                 seed) {
  future::plan(future::sequential, gc = TRUE)
  set.seed(seed)


  if (!is_empty(z_names)) {
    mediator_exposures <- list()
    for (mediator_i in seq(z_names)) {
      mediator <- z_names[mediator_i]
      task <- sl3::make_sl3_Task(
        data = at,
        covariates = c(a_names, w_names),
        outcome = mediator,
        outcome_type = family
      )

      discrete_sl_metalrn <- sl3::Lrnr_cv_selector$new(sl3::loss_squared_error)

      mediator_discrete_sl <- sl3::Lrnr_sl$new(
        learners = zeta_learner,
        metalearner = discrete_sl_metalrn
      )

      sl_fit <- mediator_discrete_sl$train(task)
      selected_learner <- sl_fit$learner_fits[[which(sl_fit$coefficients == 1)]]
      exposure_drivers <- str_extract_all(rownames(selected_learner$fit_object$glm.coefficients), paste(a_names, collapse = "|"), simplify = FALSE)
      exposure_drivers <- exposure_drivers[lapply(exposure_drivers, length) > 0]
      exposure_drivers <- exposure_drivers[!is.na(exposure_drivers)]
      mediator_exposures[[mediator]] <- unique(unlist(exposure_drivers))
    }
  } else {
    mediator_exposures <- NULL
  }

  task <- sl3::make_sl3_Task(
    data = at,
    covariates = c(a_names, z_names, w_names),
    outcome = outcome,
    outcome_type = family
  )

  discrete_sl_metalrn <- sl3::Lrnr_cv_selector$new(sl3::loss_squared_error)

  covar_discrete_sl <- sl3::Lrnr_sl$new(
    learners = zeta_learner,
    metalearner = discrete_sl_metalrn
  )

  sl_fit <- covar_discrete_sl$train(task)

  selected_learner <- sl_fit$learner_fits[[which(sl_fit$coefficients == 1)]]
  learner_name <- selected_learner$name
  full_residual_SS <- (at[, get(outcome)] - sl_fit$predict())^2


  if (grepl("earth", learner_name)) {
    best_model_basis <- as.data.frame(model.matrix(
      selected_learner$fit_object
    ))
    full_lm_mod <- stats::lm(at[, get(outcome)] ~ .,
      data = as.data.frame(best_model_basis)
    )

    anova_fit <- stats::anova(full_lm_mod)

    if (quantile_thresh != 0) {
      anova_fit <- subset(anova_fit, `F value` >
        quantile(anova_fit$`F value`,
          quantile_thresh,
          na.rm = TRUE
        ))
    }

    anova_basis <- rownames(anova_fit)
  }

  if (grepl("polspline", learner_name)) {
    best_model_basis <- polspline::design.polymars(
      selected_learner$fit_object,
      x = at[, ..covars]
    )

    best_model_basis <- best_model_basis[, -1]

    pred1 <- selected_learner$fit_object$model$pred1
    pred1 <- pred1[-1]

    pred2 <- selected_learner$fit_object$model$pred2
    pred2 <- pred2[-1]

    for (i in 1:length(covars)) {
      pred1 <- sub(paste("^", i, sep = ""), covars[i], pred1)
      pred2 <- sub(paste("^", i, sep = ""), covars[i], pred2)
    }

    knot1 <- selected_learner$fit_object$model$knot1[-1]
    knot2 <- selected_learner$fit_object$model$knot2[-1]

    splines <- cbind(
      pred1,
      knot1,
      pred2,
      knot2
    )

    colnames(splines) <- c("pred1", "knot1", "pred2", "knot2")
    splines <- as.data.frame(splines)

    splines$knot1 <- round(as.numeric(splines$knot1), 2)
    splines$knot2 <- round(as.numeric(splines$knot2), 2)

    splines[is.na(splines)] <- ""
    splines[splines == 0] <- ""

    colnames(best_model_basis) <- paste(splines$pred1, splines$knot1,
      splines$pred2, splines$knot2,
      sep = ""
    )
    polymars_model <- lm(data[, get(outcome)] ~ .,
      data = as.data.frame(best_model_basis)
    )

    anova_fit <- stats::anova(polymars_model)

    if (quantile_thresh != 0) {
      anova_fit <- subset(anova_fit, `F value` > quantile(anova_fit$`F value`,
        quantile_thresh,
        na.rm = TRUE
      ))
    }
  }

  match_list <- list()
  i <- 1
  covars <- c(a_names, z_names, w_names)
  for (j in covars) {
    matches <- stringr::str_match(rownames(anova_fit), j)
    matches[is.na(matches)] <- ""
    match_list[[i]] <- matches
    i <- i + 1
  }

  matches <- do.call(cbind, match_list)
  basis_used <- unique(apply(matches, 1, paste, collapse = ""))
  basis_used <- basis_used[basis_used != ""]


  for (basis_name_i in seq(basis_used)) {
    basis_check <- basis_used[basis_name_i]
    if (basis_check %in% names(mediator_exposures)) {
      basis_used[basis_name_i] <- paste(append(mediator_exposures[[basis_check]], basis_check), collapse = "-")
    }
  }

  results <- list(
    "learner" = selected_learner,
    "learner name" = learner_name,
    "anova fit" = anova_fit,
    "basis" = basis_used,
    "fit residuals" = full_residual_SS
  )

  return(results)
}
