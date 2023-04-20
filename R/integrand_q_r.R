#' Integrand Function for Mediation Analysis in SuperNOVA: Works for both Adaptive Quadrature and MC Methods for Integration
#'
#' Computes the integrand for mediation analysis in SuperNOVA by combining the models for outcome given mediator and exposure (q_model), and the density of the mediator given exposure and covariates (r_model).
#'
#' @param z Numeric vector, grid points for the mediator.
#' @param row_data Data frame, containing the data for a single observation.
#' @param covars Character vector, names of the covariate columns in the data.
#' @param w_names Character vector, names of the baseline covariate columns in the data.
#' @param q_model Fitted Super Learner model, the outcome model given mediator and exposure.
#' @param r_model Fitted Super Learner or parametric model, the mediator density given exposure and covariates.
#' @param exposure Character, name of the exposure variable in the data.
#' @param mediator Character, name of the mediator variable in the data.
#' @param delta Numeric, the amount to shift the exposure variable.
#' @param upper_a Numeric, the upper bound for the exposure variable.
#' @param density_type Character, either "sl" for Super Learner or "parametric" for parametric models.
#'
#' @return Numeric vector, the computed values of the integrand.
#' @export
integrand_q_r <- function(z, row_data, covars, w_names, q_model, r_model, exposure, mediator, delta, upper_a, density_type) {
  row_data <- do.call("rbind", replicate(length(z), row_data, simplify = FALSE))
  new_data_m <- new_data_r <- row_data

  new_data_m[mediator] <- z
  new_data_m[exposure] <- ifelse(new_data_m[[exposure]] + delta >= upper_a, upper_a, new_data_m[[exposure]] + delta)

  new_data_r[mediator] <- z

  task_m <- sl3::sl3_Task$new(
    data = new_data_m,
    covariates = covars,
    outcome = "y"
  )

  m_val <- q_model$predict(task_m)

  if (density_type == "sl") {
    task_r <- sl3::sl3_Task$new(
      data = new_data_r,
      covariates = c(w_names),
      outcome = mediator,
    )

    r_val <- r_model$predict(task_r)
    output <- m_val * r_val$likelihood
  } else {
    r_val <- suppressMessages(predict(r_model, new_A = new_data_r[[mediator]], new_W = new_data_r[w_names]))
    output <- m_val * r_val
  }

  return(output)
}
