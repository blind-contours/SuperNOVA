#' Estimate Probability of Censoring by Two-Phase Sampling
#'
#' @details Compute estimates of the sampling probability for inclusion in the
#'  the second-phase via the two-phase sampling mechanism. These estimates are
#'  used for the creation of inverse probability weights.
#'
#' @param V A \code{numeric} vector, \code{matrix}, \code{data.frame} or
#'  similar object giving the observed values of the covariates known to
#'  potentially inform the sampling mechanism.
#' @param C_samp A \code{numeric} vector of observed values of the indicator
#'  for inclusion in the second-phase sample.
#' @param fit_type A \code{character} indicating whether to perform the fit
#'  using GLMs or a Super Learner ensemble model. If use of Super Learner is
#'  desired, then the argument \code{sl_learners} must be provided.
#' @param sl_learners An \pkg{sl3} \code{Lrnr_sl} object, a Super Learner
#'  ensemble or learner instantiated externally using \pkg{sl3}.
#'
#' @importFrom stats glm predict binomial
#' @importFrom data.table as.data.table setnames
#' @importFrom assertthat assert_that
#'
#' @return A \code{numeric} vector of the estimated sampling mechanism.
est_samp <- function(V,
                     C_samp,
                     fit_type = c("sl", "glm"),
                     sl_learners = NULL) {
  # set defaults and check arguments
  fit_type <- match.arg(fit_type)
  if (fit_type == "sl") {
    assertthat::assert_that(!is.null(sl_learners),
      msg = paste(
        "`sl_learners` cannot be NULL",
        "when `fit_type = 'sl'`."
      )
    )
  }

  # make data objects from inputs
  fit_type <- match.arg(fit_type)
  data_in <- data.table::as.data.table(cbind(C_samp, V))
  if (!is.matrix(V)) V <- as.matrix(V)
  names_V <- paste0("V", seq_len(ncol(V)))
  data.table::setnames(data_in, c("C_samp", names_V))

  # make sl3 tasks from the data if fitting Super Learners
  if (fit_type == "sl" & !is.null(sl_learners)) {
    sl_task <- sl3::sl3_Task$new(
      data = data_in,
      covariates = names_V,
      outcome = "C_samp",
      outcome_type = "binomial"
    )
  }

  # fit logistic regression or binomial SL to get class probabilities for IPCW
  if (fit_type == "glm" & is.null(sl_learners)) {
    samp_reg <- stats::glm(
      as.formula("C_samp ~ ."),
      family = stats::binomial(),
      data = data_in
    )
    samp_probs <- unname(stats::predict(
      object = samp_reg,
      newdata = data_in,
      type = "response"
    ))
  } else if (fit_type == "sl" & !is.null(sl_learners)) {
    sl_fit <- sl_learners$train(sl_task)
    samp_probs <- as.numeric(sl_fit$predict())
  }

  # the estimated sampling weights Pi_n for creating inverse weights
  return(samp_probs)
}
