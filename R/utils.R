###############################################################################
#' @title Calculate the Joint Parameter
#' @description Using the output results for the eif use the delta method
#' to simply subtract the joint estimate from the no shift estimate. Do
#' the same for the eif to get the variance for this new parameter.
#' @param results_element Joint results for formatting
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
#' @param results_table Table of interaction results
#' @param joint_shift_fold_results Joint shift results
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
#' @param psi The Psi parameter calculated
#' @param psi_se Standard deviation for the Psi parameter
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
#' @param fold_basis_results List of lists of variable sets found in the basis
#' function data-adaptive parameter section
#' @param n_folds Number of folds used
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
#' @param common_variables Variable sets found in the data adaptive procedure
#' @param i iteration of search
#' @param a_names Names of the exposures
#' @param w_names Names of the covariates
#' @param z_names Mediator names
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
#' @param psi The estimated target parameter
#' @param variance The variance of that parameter
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
