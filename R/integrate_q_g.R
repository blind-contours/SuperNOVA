#' Integrate functions m and g using monte carlo method
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
integrate_q_g <- function(av, at, covars, w_names, q_model, g_model, exposure, g_delta, m_delta, n_samples, density_type, lower_bound, upper_bound, integration_method = "MC") {
  at <- as.data.frame(at)
  av <- as.data.frame(av)

  integrand <- function(sample_a, row_data, covars, q_model, g_model, exposure, g_delta, m_delta, upper_bound, density_type) {
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

      g_val <- g_model$predict(task_g)
      output <- m_val * g_val$likelihood
    } else {
      g_val <- suppressMessages(predict(g_model, new_A = new_data_g[[exposure]], new_W = new_data_g[w_names]))
      output <- m_val * g_val
    }

    return(output)
  }

  results <- numeric(nrow(av))

  if (integration_method == "MC") {
    sample_a <- runif(n_samples, lower_bound, upper_bound)

    for (i in 1:nrow(av)) {
      row_data <- av[i, ]

      mc_integrands <- integrand(sample_a, row_data, covars, q_model, g_model, exposure, g_delta, m_delta, upper_bound, density_type)

      integral_result <- (max(sample_a) - min(sample_a)) * mean(mc_integrands)
      results[i] <- integral_result
    }
  } else if (integration_method == "AQ") {
    for (i in 1:nrow(av)) {
      row_data <- av[i, ]

      aq_integrands <- function(a) {
        integrand(a, row_data, covars, q_model, g_model, exposure, g_delta, m_delta, upper_bound, density_type)
      }

      integral_result <- integrate(aq_integrands, lower = lower_bound, upper = upper_bound, rel.tol = 0.001)
      results[i] <- integral_result$value
    }
  }

  return(results)
}
