#' Integrate functions m and g when exposure is quantized
#'
#' @details Does the double integration as described in lemma 1
#'
#' @param data A \code{character} vector of exposures to be shifted.
#' @param covars The mediator variable
#' @param w_names Covariate names
#' @param q_model A \code{character} vector covariates to adjust for.
#' @param g_model The training data
#' @param exposure A \code{numeric} indicating the magnitude of the shift to be
#'  computed for the exposure \code{A}. This is passed to the internal
#' @param delta Object containing a set of instantiated learners from the
#'  \pkg{sl3}, to be used in fitting an ensemble model.
#'
#' @importFrom stats glm as.formula predict
#' @importFrom data.table as.data.table setnames copy set
#' @importFrom stringr str_detect
#' @importFrom assertthat assert_that
#' @export
#' @return A \code{data.table} with two columns, containing estimates of the
#'  outcome mechanism at the natural value of the exposure Q(A, W) and an
#'  upshift of the exposure Q(A + delta, W).
#'
integrate_q_g_quant <- function(av, at, covars, w_names, q_model, g_model, exposure, g_delta, m_delta, n_bins, density_type, upper_bound) {
  at <- as.data.frame(at)
  av <- as.data.frame(av)

  results <- numeric(nrow(av))

  for (i in 1:nrow(av)) {
    row_data <- av[i, ]

    bin_integrands <- integrand_q_g(sample_a = unique(av[[exposure]]), row_data, covars, q_model, g_model, exposure, g_delta, m_delta, upper_bound = upper_bound, density_type)

    sum_result <- sum(bin_integrands)
    results[i] <- sum_result
  }

  return(results)
}
