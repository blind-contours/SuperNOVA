#' Integrand function for q and g models
#'
#' @param sample_a vector of sampled exposure values
#' @param row_data dataframe containing the current row data
#' @param covars vector of covariate names
#' @param q_model sl3 learner object for the q model
#' @param g_model sl3 learner object for the g model
#' @param exposure character string indicating the exposure variable name
#' @param g_delta numeric value for the g model delta
#' @param m_delta numeric value for the m model delta
#' @param upper_bound numeric value specifying the upper bound for exposure values
#' @param density_type character string specifying the density type ("sl" or other)
#' @export
#'
#' @return output numeric vector of integrand values
integrand_q_g <- function(sample_a, row_data, covars, q_model, g_model, exposure, g_delta, m_delta, upper_bound, density_type) {
  row_data <- do.call("rbind", replicate(length(sample_a), row_data, simplify = FALSE))
  new_data_m <- new_data_g <- row_data
  new_data_m[exposure] <- ifelse(sample_a + m_delta >= upper_bound, upper_bound, sample_a + m_delta)
  new_data_g[exposure] <- ifelse(sample_a + g_delta >= upper_bound, upper_bound, sample_a + g_delta)

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

    g_val <- g_model$predict(task_g)$likelihood

    # g_val <- g_val$likelihood / sum(g_val$likelihood)
    output <- m_val * g_val
  } else {
    g_val <- suppressMessages(predict(g_model, new_A = new_data_g[[exposure]], new_W = new_data_g[w_names]))
    output <- m_val * g_val
  }

  return(output)
}
