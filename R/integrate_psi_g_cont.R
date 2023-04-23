#' Integrate Psi G for the Da part of the EIF for stochastic mediation using
#' Monte Carlo integration
#'
#' @details Computes the double integration as described in lemma 1 using
#' either Monte Carlo or adaptive quadrature methods.
#'
#' @param av A dataframe with the natural value of the exposure
#' @param at A dataframe with the upshift of the exposure
#' @param covars A vector of covariate names
#' @param w_names A vector of covariate names
#' @param q_model A fitted Q-model
#' @param r_model A fitted R-model (mediator density estimator)
#' @param g_model A fitted G-model (the training data)
#' @param exposure A string indicating the name of the exposure variable
#' @param mediator A string indicating the name of the mediator variable
#' @param delta A numeric indicating the magnitude of the shift for the exposure
#' @param n_samples A numeric specifying the number of samples for Monte Carlo integration
#' @param density_type A string specifying the type of density estimation ("sl" for Super Learner)
#' @param integration_method A string specifying the integration method to use ("MC" for Monte Carlo, "AQ" for Adaptive Quadrature)
#' @return A list with three elements: "d_a" containing the estimates of Da, "phi_aw" containing the
#'  estimates of the outcome mechanism at the natural value of the exposure Q(A, W), and "phi_aw_g" containing the
#'  estimates of the outcome mechanism at the upshift of the exposure Q(A + delta, W)
#'
#' @importFrom stats glm as.formula predict
#' @importFrom data.table as.data.table setnames copy set
#' @importFrom stringr str_detect
#' @importFrom assertthat assert_that
#' @import pracma
#' @import sl3
#' @export
integrate_psi_g_cont <- function(av, at, covars, w_names, q_model, r_model, g_model, exposure, mediator, delta, n_samples, density_type, integration_method) {
  av <- as.data.frame(av)
  at <- as.data.frame(at)

  lower_z <- floor(min(min(av[[mediator]]), min(at[[mediator]])))
  upper_z <- ceiling(max(max(av[[mediator]]), max(at[[mediator]])))

  lower_a <- min(min(av[[exposure]]), min(at[[exposure]]))
  upper_a <- max(max(av[[exposure]]), max(at[[exposure]]))

  results <- numeric(nrow(av))
  integral_inner_results <- numeric(nrow(av))
  integral_outer_results <- numeric(nrow(av))

  # Loop through each row of the av dataframe
  for (i in 1:nrow(av)) {
    row_data <- av[i, ]

    # Compute the integral_inner
    if (integration_method == "MC") {
      sample_z_inner <- runif(n_samples, lower_z, upper_z)
      sample_a <- runif(n_samples, lower_a, upper_a)
      sample_z_outer <- runif(n_samples, lower_z, upper_z)

      integrands_inner <- integrand_q_r(sample_z_inner, row_data, covars, w_names, q_model, r_model, exposure, mediator, delta, upper_a, density_type)
      integral_inner <- (max(sample_z_inner) - min(sample_z_inner)) * mean(integrands_inner)
    } else if (integration_method == "AQ") {
      integral_inner <- stats::integrate(
        function(z) integrand_q_r(z, row_data, covars, w_names, q_model, r_model, exposure, mediator, delta, upper_a, density_type),
        lower = lower_z,
        upper = upper_z,
        rel.tol = 0.001,
        subdivisions = 100,
        stop.on.error = FALSE
      )$value
    }

    # Compute the integral_outer
    if (integration_method == "MC") {
      mc_integrands_outer <- mc_integrand_q_g_r(sample_a, sample_z_outer, row_data, covars, w_names, q_model, g_model, r_model, exposure, mediator, delta, upper_a, density_type)
      average_integrand <- mean(mc_integrands_outer)
      integral_outer <- (max(sample_a) - min(sample_a)) * (max(sample_z_outer) - min(sample_z_outer)) * average_integrand
    } else if (integration_method == "AQ") {
      integral_outer <- pracma::integral2(
        function(a, z) quad_integrand_q_g_r(a, z, row_data, covars, w_names, q_model, g_model, r_model, exposure, mediator, delta, upper_a, density_type),
        xmin = lower_a,
        ymin = lower_z,
        xmax = upper_a,
        ymax = upper_z,
        reltol = 0.001
      )$Q
    }

    # Compute the results
    results[i] <- integral_inner - integral_outer
    integral_inner_results[i] <- integral_inner
    integral_outer_results[i] <- integral_outer
  }

  return(list("d_a" = results, "phi_aw" = integral_inner_results, "phi_aw_g" = integral_outer_results))
}

