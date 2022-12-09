#' #' Print Method for Counterfactual Mean of Stochastic Shift Intervention
#' #'
#' #' @details The \code{print} method for objects of class \code{txshift}.
#' #'
#' #' @param x An object of class \code{txshift}.
#' #' @param ... Other options (not currently used).
#' #' @param ci_level A \code{numeric} indicating the level of the confidence
#' #'  interval to be computed.
#' #'
#' #' @method print txshift
#' #'
#' #' @importFrom stats confint
#' #' @importFrom scales percent
#' #'
#' #' @return None. Called for the side effect of printing an informative summary
#' #'  of slots of objects of class \code{txshift}.
#' #'
#' #' @export
#' #'
#' #' @examples
#' #' if (require("sl3")) {
#' #'   set.seed(429153)
#' #'   n_obs <- 100
#' #'   W <- replicate(2, rbinom(n_obs, 1, 0.5))
#' #'   A <- rnorm(n_obs, mean = 2 * W, sd = 1)
#' #'   Y <- rbinom(n_obs, 1, plogis(A + W + rnorm(n_obs, mean = 0, sd = 1)))
#' #'   txout <- txshift(
#' #'     W = W, A = A, Y = Y, delta = 0.5,
#' #'     estimator = "tmle",
#' #'     g_exp_fit_args = list(
#' #'       fit_type = "sl",
#' #'       sl_learners_density = Lrnr_density_hse$new(Lrnr_glm$new())
#' #'     ),
#' #'     Q_fit_args = list(
#' #'       fit_type = "glm",
#' #'       glm_formula = "Y ~ ."
#' #'     )
#' #'   )
#' #'   print(txout)
#' #' }
#' print.txshift <- function(x, ..., ci_level = 0.95) {
#'   # compute confidence interval
#'   ci <- stats::confint(x, level = ci_level)
#'
#'   # construct and print output
#'   message("Counterfactual Mean of Shifted Treatment")
#'   message("Intervention: ", "Treatment + ", x$.delta)
#'   message("txshift Estimator: ", x$estimator)
#'   message("Estimate: ", round(x$psi, 4))
#'   message("Std. Error: ", round(sqrt(x$var), 4))
#'   message(paste0(
#'     scales::percent(ci_level), " CI: [",
#'     round(ci[1], 4), ", ", round(ci[3], 4), "]"
#'   ))
#' }

#' ###############################################################################
#'
#' #' Print Method for Marginal Structural Models
#' #'
#' #' @details The \code{print} method for objects of class \code{txshift_msm}.
#' #'
#' #' @param x An object of class \code{txshift_msm}.
#' #' @param ... Other options (not currently used).
#' #'
#' #' @method print txshift_msm
#' #'
#' #' @importFrom scales percent
#' #'
#' #' @return None. Called for the side effect of printing an informative summary
#' #'  of slots of objects of class \code{txshift_msm}.
#' #'
#' #' @export
#' #'
#' #' @examples
#' #' if (require("sl3")) {
#' #'   set.seed(3287)
#' #'   n_obs <- 1000
#' #'   W <- as.numeric(replicate(1, rbinom(n_obs, 1, 0.5)))
#' #'   A <- as.numeric(rnorm(n_obs, mean = 2 * W, sd = 1))
#' #'   Y <- rbinom(n_obs, 1, plogis(2 * A - W))
#' #'   msm <- msm_vimshift(
#' #'     W = W, A = A, Y = Y, estimator = "tmle",
#' #'     g_exp_fit_args = list(
#' #'       fit_type = "sl",
#' #'       sl_learners_density = Lrnr_density_hse$new(Lrnr_glm$new())
#' #'     ),
#' #'     Q_fit_args = list(
#' #'       fit_type = "glm",
#' #'       glm_formula = "Y ~ ."
#' #'     ),
#' #'     delta_grid = seq(-1, 1, 0.25)
#' #'   )
#' #'   print(msm)
#' #' }
#' print.txshift_msm <- function(x, ...) {
#'   # construct and print output
#'   message("MSM (", x$.msm_type, ") for Grid of Shifted Treatments")
#'   message(
#'     "Intervention Grid: ", "Treatment + ",
#'     paste0("{", paste(x$.delta_grid, collapse = ", "), "}")
#'   )
#'   if (x$.msm_type == "piecewise") {
#'     message("Knot Point: Shift = ", x$.msm_knot)
#'   }
#'   message("txshift MSM Estimator: ", x$estimator)
#'   if (x$.msm_type == "piecewise") {
#'     message(
#'       "Estimated Slopes: ",
#'       round(x$msm_est$param_est[2], 4), ", ",
#'       round(x$msm_est$param_est[3], 4)
#'     )
#'     message(
#'       "Std. Errors: ",
#'       round(x$msm_est$param_se[2], 4), ", ",
#'       round(x$msm_est$param_se[3], 4)
#'     )
#'     message(
#'       scales::percent(x$.ci_level), " CIs: ",
#'       "[", round(x$msm_est$ci_lwr[2], 4), ", ",
#'       round(x$msm_est$ci_upr[2], 4), "]", ", ",
#'       "[", round(x$msm_est$ci_lwr[3], 4), ", ",
#'       round(x$msm_est$ci_upr[3], 4), "]"
#'     )
#'     message(
#'       "p-values (vs. no trend): ",
#'       round(x$msm_est$p_value[2], 4), ", ",
#'       round(x$msm_est$p_value[3], 4)
#'     )
#'   } else {
#'     message("Estimated Slope: ", round(x$msm_est$param_est[2], 4))
#'     message("Std. Error: ", round(x$msm_est$param_se[2], 4))
#'     message(
#'       scales::percent(x$.ci_level), " CI: [",
#'       round(x$msm_est$ci_lwr[2], 4), ", ",
#'       round(x$msm_est$ci_upr[2], 4), "]"
#'     )
#'     message(
#'       "p-value (vs. no trend): ",
#'       round(x$msm_est$p_value[2], 4)
#'     )
#'   }
#' }
#'
#'
###############################################################################
#' @title Calculate the Joint Parameter
#' @description Using the output results for the eif use the delta method
#' to simply subtract the joint estimate from the no shift estimate. Do
#' the same for the eif to get the variance for this new parameter.
#' @export
calc_joint_results <- function(results_element) {
  psi <- results_element$psi - results_element$noshift_psi
  psi_var <- var(results_element$eif - results_element$noshift_eif) /
    length(results_element$eif)
  se_ests <- sqrt(psi_var)
  CI <- calc_CIs(psi, se_ests)
  p_vals <- calc_pvals(psi, se_ests)

  joint_results <- c(psi, psi_var, se_ests, CI[1], CI[2], p_vals)

  return(joint_results)
}

###############################################################################
#' @title Calculates the Interaction Parameter
#' @description Using the output results for the eif use the delta method
#' to simply subtract the two additive estimates from the joint shift. Do
#' the same for the eif to get the variance for this new parameter.
#' @export

calc_intxn_results <- function(results_table, joint_shift_fold_results) {
  intxn_psi <- results_table[[3, 1]] -
    results_table[[2, 1]] -
    results_table[[1, 1]]

  psi_variance <- var((joint_shift_fold_results[[3]]$eif -
    joint_shift_fold_results[[3]]$noshift_eif) -
    (joint_shift_fold_results[[2]]$eif -
      joint_shift_fold_results[[2]]$noshift_eif) -
    (joint_shift_fold_results[[1]]$eif -
      joint_shift_fold_results[[1]]$noshift_eif)) /
    length(joint_shift_fold_results[[1]]$eif)

  psi_se <- sqrt(psi_variance)

  CI <- calc_CIs(intxn_psi, psi_se)
  p_vals <- calc_pvals(intxn_psi, psi_se)

  intxn_results <- c(intxn_psi, psi_variance, psi_se, CI[1], CI[2], p_vals)

  return(intxn_results)
}

###############################################################################
#' @title Calculate Confidence Interals
#' @description Gives information on which variable sets are used in the
#' basis function in the data-adaptive estimation section
#' @export
calc_CIs <- function(psi, psi_se) {
  psi_CIs <- c(
    round(psi + stats::qnorm(0.05 / 2, lower.tail = T) * psi_se, 4),
    round(psi + stats::qnorm(0.05 / 2, lower.tail = F) * psi_se, 4)
  )
  return(psi_CIs)
}


###############################################################################
#' @title Calculate the frequency basis functions are used in the folds
#' @description Gives information on which variable sets are used in the
#' basis function in the data-adaptive estimation section
#'
#' @export
calc_basis_freq <- function(fold_basis_results, n_folds) {
  fold_basis_table <- table(unlist(fold_basis_results))
  fold_basis_table <- round(fold_basis_table / n_folds, 2)

  return(fold_basis_table)
}



###############################################################################
#' @title Extract variables from the basis function search
#' @description This simple function creates a list of target variables
#' that were used in the data adpative procedure.
#'
#' @export
extract_vars_from_basis <- function(common_variables, i, a_names,
                                    w_names, z_names) {
  target <- common_variables[i]
  match_list <- list()

  for (j in 1:length(c(a_names, z_names, w_names))) {
    var <- c(a_names, z_names, w_names)[j]
    matches <- stringr::str_match(target, var)
    matches[is.na(matches)] <- ""
    match_list[j] <- matches[[1]]
  }

  matches <- unlist(match_list[!match_list %in% ""])
  return(list("matches" = matches, "target" = target))
}

###############################################################################
#' @title Calculate the p-values based on variance of the influence curve
#' @description Calculates the p-value using the standard formula
#'
#' @export
calc_pvals <- function(psi, variance) {
  p_value <- 2 * stats::pnorm(abs(psi / sqrt(variance)), lower.tail = F)
  return(p_value)
}

is.SuperNOVA <- function(x) {
  class(x) == "SuperNOVA"
}

is.SuperNOVA <- function(x) {
  class(x) == "SuperNOVA_msm"
}

###################################################################
#' @title Get rules from partykit object in rule fitting
#' @param x Partykit glmtree model object
#' @param i null
#' @param ... additional arguments
#' @return List of rules
#'
#' @export
# Copied from internal partykit function
list_rules_party <- function(x, i = NULL, ...) {
  if (is.null(i)) {
    i <- partykit::nodeids(x, terminal = TRUE)
  }
  if (length(i) > 1) {
    ret <- sapply(i, list_rules_party, x = x)
    names(ret) <- if (is.character(i)) {
      i
    } else {
      names(x)[i]
    }
    return(ret)
  }
  if (is.character(i) && !is.null(names(x))) {
    i <- which(names(x) %in% i)
  }
  stopifnot(length(i) == 1 & is.numeric(i))
  stopifnot(i <= length(x) & i >= 1)
  i <- as.integer(i)
  dat <- partykit::data_party(x, i)
  if (!is.null(x$fitted)) {
    findx <- which("(fitted)" == names(dat))[1]
    dat <- dat[, -(findx:ncol(dat)), drop = FALSE]
    if (ncol(dat) == 0) {
      dat <- x$data
    }
  } else {
    dat <- x$data
  }
  rule <- c()
  rec_fun <- function(node) {
    if (partykit::id_node(node) == i) {
      return(NULL)
    }
    kid <- sapply(partykit::kids_node(node), partykit::id_node)
    whichkid <- max(which(kid <= i))
    split <- partykit::split_node(node)
    ivar <- partykit::varid_split(split)
    svar <- names(dat)[ivar]
    index <- partykit::index_split(split)
    if (is.factor(dat[, svar])) {
      if (is.null(index)) {
        index <- ((1:nlevels(dat[, svar])) > partykit::breaks_split(split)) +
          1
      }
      slevels <- levels(dat[, svar])[index == whichkid]
      srule <- paste(svar, " %in% c(\"", paste(slevels,
        collapse = "\", \"", sep = ""
      ), "\")", sep = "")
    } else {
      if (is.null(index)) {
        index <- seq_along(kid)
      }
      breaks <- cbind(c(-Inf, partykit::breaks_split(split)), c(
        partykit::breaks_split(split),
        Inf
      ))
      sbreak <- breaks[index == whichkid, ]
      right <- partykit::right_split(split)
      srule <- c()
      if (is.finite(sbreak[1])) {
        srule <- c(srule, paste(svar, ifelse(right, ">",
          ">="
        ), sbreak[1]))
      }
      if (is.finite(sbreak[2])) {
        srule <- c(srule, paste(svar, ifelse(right, "<=",
          "<"
        ), sbreak[2]))
      }
      srule <- paste(srule, collapse = " & ")
    }
    rule <<- c(rule, srule)
    return(rec_fun(node[[whichkid]]))
  }
  node <- rec_fun(partykit::node_party(x))
  paste(rule, collapse = " & ")
}
