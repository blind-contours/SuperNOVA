#' Integrate functions m and g using monte carlo method
#'
#' @details Does the double integration as described in lemma 1
#'
#' @param at Training data
#' @param av Validation data
#' @param covars The mediator variable
#' @param w_names Covariate names
#' @param q_model A \code{character} vector covariates to adjust for.
#' @param g_model The training data
#' @param exposure A \code{numeric} indicating the magnitude of the shift to be
#'  computed for the exposure \code{A}. This is passed to the internal
#' @param g_delta Shift in exposure in the g model
#' @param m_delta Shift in exposure in the m model
#' @param n_samples Number of samples for MC integration
#' @param density_type Type of density estimation
#' @param lower_bound Lower bound of exposure
#' @param upper_bound Upper bound of exposure
#' @param integration_method Type of integration method to use
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
integrate_q_g <- function(av, at, covars, w_names, q_model, g_model, exposure, g_delta, m_delta, n_samples, density_type, lower_bound, upper_bound, integration_method = "MC") {
  at <- as.data.frame(at)
  av <- as.data.frame(av)

  if (integration_method == "MC") {
    sample_a <- runif(n_samples, lower_bound, upper_bound)

    results <- sapply(1:nrow(av), function(i) {
      row_data <- av[i, ]

      mc_integrands <- integrand_q_g(sample_a, row_data, covars, w_names, q_model, g_model, exposure, g_delta, m_delta, upper_bound, density_type)

      integral_result <- (max(sample_a) - min(sample_a)) * mean(mc_integrands)
      return(integral_result)
    })
  } else if (integration_method == "AQ") {
    results <- sapply(1:nrow(av), function(i) {
      row_data <- av[i, ]

      integral_result <- stats::integrate(
        function(sample_a) integrand_q_g(sample_a, row_data, covars, q_model, g_model, exposure, g_delta, m_delta, upper_bound, density_type),
        lower = lower_bound,
        upper = upper_bound,
        rel.tol = 0.001,
        subdivisions = 300,
        stop.on.error = FALSE
      )$value

      return(integral_result)
    })
  }

  return(results)
}
