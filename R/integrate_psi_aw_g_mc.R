#' Pseudo regression integration using monte carlo approaches
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

integrate_psi_aw_g_mc <- function(at, av, covars, w_names, pseudo_model, g_model, exposure, delta, psi_aw, n_samples = 1000) {
  at <- as.data.frame(at)
  av <- as.data.frame(av)
  lower <- floor(min(min(av[exposure]), min(at[exposure])))
  upper <- ceiling(max(max(av[exposure]), max(at[exposure])))

  integrand <- function(a, row_data, covars, pseudo_model, g_model, exposure, delta, upper) {
    row_data <- do.call("rbind", replicate(length(a), row_data, simplify = FALSE))
    new_data_m <- new_data_g <- row_data
    new_data_m[exposure] <- a + delta #ifelse(a + delta >= upper, upper, a + delta)
    new_data_g[exposure] <- a

    task_m <- sl3::sl3_Task$new(
      data = new_data_m,
      covariates = covars,
      outcome = "pseudo_outcome")

    task_g <- sl3::sl3_Task$new(
      data = new_data_g,
      covariates = c(w_names),
      outcome = exposure,
    )

    m_val <- pseudo_model$predict(task_m)
    g_val <- g_model$predict(task_g)
    output <- m_val * g_val$likelihood
    return(output)
  }

  results <- numeric(nrow(av))

  for (i in 1:nrow(av)) {
    row_data <- av[i, ]
    sample_a <- runif(n_samples, lower, upper)
    integral_values <- integrand(sample_a, row_data, covars, pseudo_model, g_model, exposure, delta, upper)
    integral_result <- mean(integral_values) * (upper - lower)

    results[i] <- psi_aw[i] - integral_result
  }

  return(results)
}
