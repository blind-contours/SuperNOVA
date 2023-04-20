#' Compute the Monte Carlo integrand for double integration in SuperNOVA
#'
#' This function computes the Monte Carlo integrand for double integration in the SuperNOVA package.
#' It combines the models for the outcome given mediator and exposure (q_model),
#' the density of the mediator given exposure and covariates (r_model),
#' and the density of the exposure given covariates (g_model).
#'
#' @param sample_a A vector of sampled exposure values.
#' @param sample_z A vector of sampled mediator values.
#' @param row_data A single row of data containing covariates and outcome.
#' @param covars A character vector of covariate names.
#' @param w_names A character vector of baseline covariate names.
#' @param q_model An SL object representing the outcome model.
#' @param g_model An SL object representing the density model of exposure given covariates.
#' @param r_model An SL object representing the density model of mediator given exposure and covariates.
#' @param exposure The name of the exposure variable.
#' @param mediator The name of the mediator variable.
#' @param delta The exposure shift for the counterfactual calculation.
#' @param upper_a The upper bound for the exposure.
#' @param density_type A string indicating the density estimation method, either "sl" or "nonparametric".
#'
#' @return A numeric vector representing the Monte Carlo integrand values for the given sample points.
mc_integrand_q_g_r <- function(sample_a, sample_z, row_data, covars, w_names, q_model, g_model, r_model, exposure, mediator, delta, upper_a, density_type) {
  rep <- length(sample_z)
  row_data <- do.call("rbind", replicate(rep, row_data, simplify = FALSE))
  new_data_m <- new_data_g <- new_data_r <- row_data

  new_data_m[exposure] <- ifelse(sample_a + delta >= upper_a, upper_a, sample_a + delta)
  new_data_m[mediator] <- sample_z

  new_data_g[exposure] <- sample_a
  new_data_r[mediator] <- sample_z

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
    r_val <- suppressMessages(predict(r_model, new_A = new_data_g[[mediator]], new_W = new_data_r[w_names]))
    g_val <- suppressMessages(predict(g_model, new_A = new_data_r[[exposure]], new_W = new_data_r[w_names]))
    output <- m_val * g_val * r_val
  }

  return(output)
}
