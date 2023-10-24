#' @title Fit the Zeta Learner: A Highly Flexible Estimator Using Splines or Basis Functions
#'
#' @description Identifies variable sets for shift interventions by fitting a flexible estimator
#' using a discrete Super Learner. This function non-parametrically performs variable importance
#' of variable sets, such as identifying two mixture components or mixture components and baseline
#' covariates that explain the outcome given an ANOVA-like decomposition. An F-statistic is used
#' to threshold the importance of variable sets, and these variable sets are treated as the data-
#' adaptive parameter.
#'
#' @param at Training dataframe
#' @param a_names Names of treatment variables
#' @param z_names Names of mediator  variables
#' @param w_names Names of baseline covariate variables
#' @param outcome Variable name for the outcome
#' @param outcome_type Type of the outcome variable (e.g., "continuous", "binary")
#' @param quantile_thresh Quantile level to set the F-statistic threshold for determining
#' basis functions for estimation
#' @param zeta_learner Stack of algorithms made in SL3 used in ensemble machine
#' learning to fit Y|A,W
#' @param fold Current fold in the cross-validation
#' @param seed Seed number for consistent results
#' @param outcome_type Variable type of the outcome
#' @param mediator_type Variable type of the mediator
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
#' the ANOVA fit on the model matrix of the basis functions, the basis used,
#' and the residual metrics
#' @export


fit_basis_estimators <- function(at,
                                 a_names,
                                 z_names,
                                 w_names,
                                 outcome,
                                 outcome_type,
                                 mediator_type,
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
        outcome_type = mediator_type
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
    outcome_type = outcome_type
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


  if (grepl("ranger", learner_name)) {
    # Initialize a list to hold terminal leaf information

    fit <- selected_learner$fit_object

    # Function to extract the path leading to a node
    get_path <- function(node_id, tree) {
      path <- c()

      while(node_id != 0) {
        parent_id <- which(tree$leftChild == node_id | tree$rightChild == node_id)
        if (length(parent_id) == 0) {
          break
        }

        path <- c(path, tree$splitvarName[parent_id])
        node_id <- parent_id
      }

      return(rev(path))
    }

    # Extract terminal nodes and their paths for a specific tree
    extract_paths_from_tree <- function(tree) {
      terminal_nodes <- tree$nodeID[tree$terminal]
      paths <- lapply(terminal_nodes, function(node) get_path(node, tree))
      names(paths) <- terminal_nodes
      return(paths)
    }

    # Using treeInfo on your fit to get the tree details for the first tree (as an example)
    tree <- treeInfo(fit, 1)

    # Extract terminal node paths for the first tree
    paths_for_tree_1 <- extract_paths_from_tree(tree)

    # Print the paths
    str(paths_for_tree_1)



  }

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

    get_pred_knot <- function(i, covars, pred, knot) {
      for (j in 1:length(covars)) {
        pred_match <- sub(paste("^", j, sep = ""), covars[j], pred)
      }
      return(list(pred = pred, knot = round(as.numeric(knot), 2)))
    }

    covars <- colnames(selected_learner$training_task$X)
    # n_covars <- length(covars)

    best_model_basis <- polspline::design.polymars(
      selected_learner$fit_object,
      x = selected_learner$fit_object$call$predictors
    )

    best_model_basis <- best_model_basis[, -1]

    if (is.null(dim(best_model_basis))) {
      single_pred <- covars[selected_learner$fit_object$model$pred1[-1]]

      anova_fit <- as.data.frame(single_pred)
      rownames(anova_fit) <- anova_fit

    } else {
      pred1_res <- get_pred_knot(1, covars, selected_learner$fit_object$model$pred1[-1],
                                 selected_learner$fit_object$model$knot1[-1])

      pred2_res <- get_pred_knot(2, covars, selected_learner$fit_object$model$pred2[-1],
                                 selected_learner$fit_object$model$knot2[-1])

      splines <- data.frame(
        pred1 = pred1_res$pred,
        knot1 = pred1_res$knot,
        pred2 = pred2_res$pred,
        knot2 = pred2_res$knot
      )

      splines[is.na(splines)] <- ""
      splines[splines == 0] <- ""

      colnames(best_model_basis) <- paste(splines$pred1, splines$knot1,
                                          splines$pred2, splines$knot2,
                                          sep = "")

      polymars_model <- lm(data[, get(outcome)] ~ ., data = as.data.frame(best_model_basis))

      anova_fit <- stats::anova(polymars_model)

      if (quantile_thresh != 0) {
        anova_fit <- subset(anova_fit, `F value` > quantile(anova_fit$`F value`,
                                                            quantile_thresh,
                                                            na.rm = TRUE))
      }
    }
  }

  # Extract variables from a basis function
  extract_vars_from_basis <- function(basis) {
    matches <- str_extract_all(basis, paste(covars, collapse = "|"))
    return(paste0(unique(unlist(matches)), collapse = "-"))
  }

  # Define covariates
  covars <- c(a_names, z_names, w_names)

  # Extract variables from each basis function and store them
  basis_used <- sapply(rownames(anova_fit), extract_vars_from_basis)

  # Filter out any empty strings
  basis_used <- basis_used[basis_used != ""]

  # Get unique variable sets
  unique_sets <- unique(basis_used)

  basis_used <- unique_sets

  exposure_mediator_pairs <- list()
  for (basis_name_i in seq(basis_used)) {
    basis_check <- basis_used[basis_name_i]
    if (basis_check %in% names(mediator_exposures)) {
      for (exposure_i in mediator_exposures[[basis_check]]) {
        exposure_mediator_pairs <- append(exposure_mediator_pairs, paste(append(exposure_i, basis_check), collapse = "-"))
      }
    }
  }

  basis_used <- c(basis_used, unlist(exposure_mediator_pairs))

  results <- list(
    "learner" = selected_learner,
    "learner name" = learner_name,
    "anova fit" = anova_fit,
    "basis" = basis_used,
    "fit residuals" = full_residual_SS
  )

  return(results)
}
