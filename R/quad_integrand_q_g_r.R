#' Integrand Function for Q, g, r Model Adaptive Quadrature Double Integration in SuperNOVA
#'
#' Computes the integrand for adaptive quadrature double integration in SuperNOVA by combining the models for outcome given mediator and exposure (q_model), the density of the mediator given exposure and covariates (r_model), and the density of the exposure given covariates (g_model).
#'
#' @param a Numeric matrix, grid points for the exposure.
#' @param z Numeric matrix, grid points for the mediator.
#' @param row_data Data frame, containing the data for a single observation.
#' @param covars Character vector, names of the covariate columns in the data.
#' @param w_names Character vector, names of the baseline covariate columns in the data.
#' @param q_model Fitted Super Learner model, the outcome model given mediator and exposure.
#' @param g_model Fitted Super Learner or parametric model, the exposure density given covariates.
#' @param r_model Fitted Super Learner or parametric model, the mediator density given exposure and covariates.
#' @param exposure Character, name of the exposure variable in the data.
#' @param mediator Character, name of the mediator variable in the data.
#' @param delta Numeric, the amount to shift the exposure variable.
#' @param upper_a Numeric, the upper bound for the exposure variable.
#' @param density_type Character, either "sl" for Super Learner or "hal" for highly adaptive lasso
#'
#' @return Numeric matrix, the computed values of the integrand.
#' @export
quad_integrand_q_g_r <- function(a, z, row_data, covars, w_names, q_model, g_model, r_model, exposure, mediator, delta, upper_a, density_type) {
  rep <- dim(z)[2]
  row_data <- do.call("rbind", replicate(rep, row_data, simplify = FALSE))
  new_data_m <- new_data_g <- new_data_r <- row_data

  output_matrix <- matrix(0, nrow = rep, ncol = rep)

  for (col in seq(ncol(z))) {
    a_vec <- a[col, ]
    z_vec <- z[col, ]

    new_data_m[exposure] <- ifelse(a_vec + delta >= upper_a, upper_a, a_vec + delta)
    new_data_m[mediator] <- z_vec

    new_data_g[exposure] <- a_vec
    new_data_r[mediator] <- z_vec

    task_m <- sl3::sl3_Task$new(
      data = new_data_m,
      covariates = covars,
      outcome = "y"
    )

    m_val <- q_model$predict(task_m)

    if (density_type == "sl") {
      task_g <- sl3::sl3_Task$new(
        data = new_data_g,
        covariates = c(w_names),
        outcome = exposure,
      )

      task_r <- sl3::sl3_Task$new(
        data = new_data_r,
        covariates = c(w_names),
        outcome = mediator,
      )

      g_val <- g_model$predict(task_g)
      r_val <- r_model$predict(task_r)

      output <- m_val * g_val$likelihood * r_val$likelihood
    } else {
      r_val <- suppressMessages(predict(r_model, new_A = new_data_g[[exposure]], new_W = new_data_r[w_names]))
      g_val <- suppressMessages(predict(g_model, new_A = new_data_r[[mediator]], new_W = new_data_r[w_names]))
      output <- m_val * g_val * r_val
    }

    output_matrix[, col] <- output
  }

  return(output_matrix)
}
