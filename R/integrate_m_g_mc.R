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
integrate_m_g_mc <- function(av, at, covars, w_names, q_model, g_model, exposure, g_delta, m_delta, n_samples) {
  at <- as.data.frame(at)
  av <- as.data.frame(av)

  lower <- floor(min(min(at[[exposure]]), min(av[[exposure]])))
  upper <- ceiling(max(max(at[[exposure]]), max(av[[exposure]])))

  integrand <- function(sample_a, row_data, covars, q_model, g_model, exposure, g_delta, m_delta, upper) {
    row_data <- do.call("rbind", replicate(length(sample_a), row_data, simplify = FALSE))
    new_data_m <- new_data_g <- row_data
    new_data_m[exposure] <- sample_a + m_delta #ifelse(sample_a + m_delta >= upper, upper, sample_a + m_delta)
    new_data_g[exposure] <- sample_a + g_delta #ifelse(sample_a + g_delta >= upper, upper, sample_a + g_delta)

    task_m <- sl3::sl3_Task$new(
      data = new_data_m,
      covariates = covars,
      outcome = "y")

    task_g <- sl3::sl3_Task$new(
      data = new_data_g,
      covariates = c(w_names),
      outcome = exposure,
    )

    m_val <- q_model$predict(task_m)

    g_val <- g_model$predict(task_g)
    output <- m_val * g_val$likelihood

    return(output)
  }

  results <- numeric(nrow(av))

  for (i in 1:nrow(av)) {
    row_data <- av[i, ]

    sample_a <- runif(n_samples, lower, upper)
    mc_integrands <- integrand(sample_a, row_data, covars, q_model, g_model, exposure, g_delta, m_delta, upper)

    integral_result <- (upper - lower) * mean(mc_integrands)
    results[i] <- integral_result
  }

  return(results)
}
