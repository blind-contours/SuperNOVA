fit_basis_estimators <- function(At,
                                 covars,
                                 exposures,
                                 outcome,
                                 family,
                                 quantile_thresh,
                                 Exposures_stack,
                                 Covariate_stack,
                                 fold,
                                 max_iter,
                                 verbose) {
  future::plan(future::sequential, gc = TRUE)

  task <- sl3::make_sl3_Task(
    data = At, covariates = covars,
    outcome = outcome,
    outcome_type = family
  )

  discrete_sl_metalrn <- sl3::Lrnr_cv_selector$new()

  covar_discrete_sl <- sl3::Lrnr_sl$new(
    learners = Covariate_stack,
    metalearner = discrete_sl_metalrn
  )
  # delayed_sl_fit <- delayed_learner_train(discrete_sl, task)

  sl_fit <- covar_discrete_sl$train(task)

  ## calculate remaining variance unexplained by W in residuals
  QbarW_initial <- sl_fit$predict()

  At <- cbind(At, QbarW_initial)

  task <- sl3::make_sl3_Task(
    data = At, covariates = exposures,
    outcome = outcome,
    outcome_type = family
  )

  exposure_discrete_sl <- sl3::Lrnr_sl$new(
    learners = Exposures_stack,
    metalearner = discrete_sl_metalrn,
  )

  sl_fit <- exposure_discrete_sl$train(task)

  QbarAW_initial <- sl_fit$predict()

  At <- cbind(At, QbarAW_initial)

  iter <- 0
  stop <- FALSE

  At_no_offset <- data.table::copy(At)
  At_no_offset$QbarW_initial <- 0
  At_no_offset$QbarAW_initial <- 0

  while (stop == FALSE) {
    iter <- iter + 1

    offset_task <- sl3::make_sl3_Task(
      data = At, covariates = covars,
      outcome = outcome,
      outcome_type = family,
      offset = "QbarAW_initial"
    )

    no_offset_task <- sl3::sl3_Task$new(
      data = At_no_offset,
      covariates = covars,
      outcome = outcome,
      outcome_type = family,
      offset = "QbarAW_initial"
    )

    sl_fit_backfit_offset <- covar_discrete_sl$train(offset_task)

    # preds_offset <- sl_fit_backfit$predict()
    g_preds_no_offset <- sl_fit_backfit_offset$predict(no_offset_task)
    g_preds_offset <- sl_fit_backfit_offset$predict(offset_task)

    At[, "QbarW_now"] <- g_preds_no_offset

    task_offset <- sl3::make_sl3_Task(
      data = At,
      covariates = exposures,
      outcome = outcome,
      outcome_type = family,
      offset = "QbarW_initial"
    )

    no_offset_task <- sl3::sl3_Task$new(
      data = At_no_offset,
      covariates = exposures,
      outcome = outcome,
      outcome_type = family,
      offset = "QbarW_initial"
    )

    sl_fit <- exposure_discrete_sl$train(task_offset)

    f_preds_offset <- sl_fit$predict(task_offset)
    f_preds_no_offset <- sl_fit$predict(no_offset_task)

    At[, "QbarAW_now"] <- f_preds_no_offset

    curr_diff <- abs(g_preds_offset - f_preds_offset)

    if (verbose) {
      if (iter == 1) {
        print(paste(
          "Fold: ", fold, "|",
          "Process: ", "Exposure Basis Backfitting", "|",
          "Iteration: ", iter, "|",
          "Delta: ", "None", "|",
          "Diff: ", mean(curr_diff)
        ))
      } else {
        print(paste(
          "Fold: ", fold, "|",
          "Process: ", "Exposure Basis Backfitting", "|",
          "Iteration: ", iter, "|",
          "Delta: ", mean(curr_diff - prev_diff), "|",
          "Diff: ", mean(curr_diff)
        ))
      }
    }

    At[, "QbarW_initial"] <- At$QbarW_now
    At[, "QbarAW_initial"] <- At$QbarAW_now

    if (iter == 1) {
      stop <- FALSE
      prev_diff <- curr_diff
    } else if (abs(mean(curr_diff - prev_diff)) <= 0.001) {
      stop <- TRUE
    } else if (iter >= max_iter) {
      stop <- TRUE
    } else {
      stop <- FALSE
      prev_diff <- curr_diff
    }
  }

  # task <- make_sl3_Task(
  #   data = data, covariates = covars,
  #   outcome = outcome, outcome_type = "continuous"
  # )
  #
  # discrete_sl_metalrn <- Lrnr_cv_selector$new()
  #
  # discrete_sl <- Lrnr_sl$new(
  #   learners = Q1_stack,
  #   metalearner = discrete_sl_metalrn
  # )
  #
  # sl_fit <- discrete_sl$train(task)
  selected_learner <- sl_fit$learner_fits[[which(sl_fit$coefficients == 1)]]
  learner_name <- selected_learner$name
  full_residual_SS <- (At[, get(outcome)] - f_preds_offset)^2

  if (grepl("earth", learner_name)) {
    best_model_basis <- as.data.frame(model.matrix(selected_learner$fit_object))
    # bx <- best_model_basis[, -1]
    # colnames(bx) <- colnames(best_model_basis)[-1]
    full_lm_mod <- lm(At[, get(outcome)] ~ ., data = as.data.frame(best_model_basis))

    anova_fit <- anova(full_lm_mod)

    if (quantile_thresh != 0) {
      anova_fit <- subset(anova_fit, `F value` > quantile(anova_fit$`F value`, quantile_thresh, na.rm = TRUE))
    }

    anova_basis <- rownames(anova_fit)
  }

  if (grepl("poly", learner_name)) {
    best_model_basis <- polspline::design.polymars(selected_learner$fit_object, x = data[, ..covars])
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

    colnames(best_model_basis) <- paste(splines$pred1, splines$knot1, splines$pred2, splines$knot2, sep = "")
    polymars_model <- lm(data[, get(outcome)] ~ ., data = as.data.frame(best_model_basis))
    anova_fit <- anova(polymars_model)

    if (quantile_thresh != 0) {
      anova_fit <- subset(anova_fit, `F value` > quantile(anova_fit$`F value`, quantile_thresh, na.rm = TRUE))
    }
  }

  match_list <- list()
  i <- 1
  for (j in exposures) {
    matches <- stringr::str_match(rownames(anova_fit), j)
    matches[is.na(matches)] <- ""
    match_list[[i]] <- matches
    i <- i + 1
  }

  matches <- do.call(cbind, match_list)
  basis_used <- unique(apply(matches, 1, paste, collapse = ""))
  basis_used <- basis_used[basis_used != ""]

  results <- list(
    "learner" = selected_learner,
    "learner name" = learner_name,
    "anova fit" = anova_fit,
    "basis" = basis_used,
    "fit residuals" = full_residual_SS
  )

  return(results)
}
