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
integrate_m_g_quant <- function(av, at, covars, w_names, q_model, g_model, exposure, g_delta, m_delta, n_bins) {
  at <- as.data.frame(at)
  av <- as.data.frame(av)

  integrand <- function(bin_a, row_data, covars, q_model, g_model, exposure, g_delta, m_delta, upper, av) {
    # row_data <- do.call("rbind", replicate(length(bin_a), row_data, simplify = FALSE))
    new_data_m <- new_data_g <- row_data
    new_data_m[exposure] <- ifelse(bin_a + m_delta >= upper, upper, bin_a + m_delta)
    new_data_g[exposure] <- ifelse(bin_a + g_delta >= upper, upper, bin_a + g_delta)

    task_m <- sl3::sl3_Task$new(
      data = new_data_m,
      covariates = covars,
      outcome = "y"
    )

    task_g <- sl3::sl3_Task$new(
      data = new_data_g,
      covariates = c(w_names),
      outcome = exposure,
    )

    m_val <- q_model$predict(task_m)
    g_val <- g_model$predict(task_g)

    index <- ifelse(bin_a + g_delta > upper, upper, bin_a + g_delta)
    g_val <- unlist(g_val)[[index]]
    g_val <- ifelse(g_val <= 1 / sqrt(nrow(av)), 1 / sqrt(nrow(av)), g_val)

    output <- m_val * g_val # Use the probability corresponding to the current bin_a
    return(output)
  }

  results <- numeric(nrow(av))

  for (i in 1:nrow(av)) {
    row_data <- av[i, ]

    bin_integrands <- sapply(1:n_bins, function(bin_a) {
      integrand(bin_a, row_data, covars, q_model, g_model, exposure, g_delta, m_delta, upper = n_bins, av)
    })

    sum_result <- sum(bin_integrands)
    results[i] <- sum_result
  }

  return(results)
}
