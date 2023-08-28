#' Integrate Psi G for the Da part of the EIF for stochastic mediation
#'
#' @details Does the double integration as described in lemma 1
#'
#' @param at Training data
#' @param av Validation
#' @param covars Covariates used in the m model
#' @param w_names Covariate used in the g model
#' @param q_model A \code{character} vector covariates to adjust for.
#' @param r_model Mediator density estimator
#' @param g_model The training data
#' @param exposure A \code{numeric} indicating the magnitude of the shift to be
#'  computed for the exposure \code{A}. This is passed to the internal
#' @param mediator The mediator variable name
#'  \code{\link{shift_additive}} and is currently limited to additive shifts.
#' @param delta Object containing a set of instantiated learners from the
#'  \pkg{sl3}, to be used in fitting an ensemble model.
#' @param integration_method Type of integration method to be used
#' @param n_samples Number of MC samples for MC integration
#' @param density_type Type of density estimation
#'
#' @importFrom stats glm as.formula predict
#' @importFrom data.table as.data.table setnames copy set
#' @importFrom stringr str_detect
#' @importFrom assertthat assert_that
#' @importFrom pracma integral2
#' @export
#' @return A \code{data.table} with two columns, containing estimates of the
#'  outcome mechanism at the natural value of the exposure Q(A, W) and an
#'  upshift of the exposure Q(A + delta, W).

integrate_psi_g <- function(at, av, covars, w_names, q_model, r_model, g_model, exposure, mediator, delta, integration_method, n_samples, density_type) {
  av <- as.data.frame(av)

  lower_z <- min(av[[mediator]])
  upper_z <- max(av[[mediator]])

  lower_a <- min(av[[exposure]])
  upper_a <- max(av[[exposure]])

  results <- numeric(nrow(av))
  integral_inner_results <- numeric(nrow(av))

  for (i in 1:nrow(av)) {
    row_data <- av[i, ]

    if (integration_method == "MC") {
      sample_z <- runif(n_samples, lower_z, upper_z)
      sample_a <- runif(n_samples, lower_a, upper_a)

      mc_integrands_inner <- integrand_q_r(sample_z, row_data, covars, w_names, q_model, r_model, exposure, mediator, delta, upper_a, density_type)
      integral_inner <- (max(sample_z) - min(sample_z)) * mean(mc_integrands_inner)

      mc_integrands_outer <- mc_integrand_q_g_r(sample_a, sample_z, row_data, covars, w_names, q_model, g_model, r_model, exposure, mediator, delta, upper_a, density_type)
      integral_outer <- (max(sample_a) - min(sample_a)) * (max(sample_z) - min(sample_z)) * mean(mc_integrands_outer)
    } else {
      integral_inner <- stats::integrate(
        function(z) integrand_q_r(z, row_data, covars, w_names, q_model, r_model, exposure, mediator, delta, upper_a, density_type),
        lower = lower_z,
        upper = upper_z,
        rel.tol = 0.001,
        subdivisions = 300,
        stop.on.error = FALSE
      )$value

      integral_outer <- integral2(
        function(a, z) quad_integrand_q_g_r(a, z, row_data, covars, w_names, q_model, g_model, r_model, exposure, mediator, delta, upper_a, density_type),
        xmin = lower_a,
        ymin = lower_z,
        xmax = upper_a,
        ymax = upper_z,
        reltol = 0.001
      )$Q
    }
  }

  results[i] <- integral_inner - integral_outer
  integral_inner_results[i] <- integral_inner


  return(list("d_a" = results, "phi_aw" = integral_inner_results))
}
