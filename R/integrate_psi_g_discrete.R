#' Integrate Psi G for the Da part of the EIF for stochastic mediation
#' when A is quartiles
#'
#' @details Does the double integration as described in lemma 1
#'
#' @param data A \code{character} vector of exposures to be shifted.
#' @param covars The mediator variable
#' @param w_names Covariate names
#' @param q_model A \code{character} vector covariates to adjust for.
#' @param r_model Mediator density estimator
#' @param g_model The training data
#' @param exposure A \code{numeric} indicating the magnitude of the shift to be
#'  computed for the exposure \code{A}. This is passed to the internal
#' @param mediator The mediator variable name
#'  \code{\link{shift_additive}} and is currently limited to additive shifts.
#' @param delta Object containing a set of instantiated learners from the
#'  \pkg{sl3}, to be used in fitting an ensemble model.
#' @param bins A \code{dataframe} of training data specific to the fold
#'
#' @importFrom stats glm as.formula predict
#' @importFrom data.table as.data.table setnames copy set
#' @importFrom stringr str_detect
#' @importFrom assertthat assert_that
#' @import pracma
#' @export
#' @return A \code{data.table} with two columns, containing estimates of the
#'  outcome mechanism at the natural value of the exposure Q(A, W) and an
#'  upshift of the exposure Q(A + delta, W)

integrate_psi_g_discrete <- function(av, at, covars, w_names, q_model, r_model, g_model, exposure, mediator, delta, n_bins, n_samples = 1000) {
  av <- as.data.frame(av)
  at <- as.data.frame(at)

  lower_z <- floor(min(min(av[[mediator]]), min(at[[mediator]])))
  upper_z <- ceiling(max(max(av[[mediator]]), max(at[[mediator]])))

  integrand_m_r <- function(sample_z_inner, row_data, covars, w_names, q_model, r_model, exposure, mediator, delta, upper_a) {
    row_data_rep <- do.call("rbind", replicate(length(sample_z_inner), row_data, simplify = FALSE))
    new_data_m <- new_data_r <- row_data_rep

    new_data_m[mediator] <- sample_z_inner
    new_data_m[exposure] <- ifelse(new_data_m[[exposure]] + delta >= upper_a, upper_a, new_data_m[[exposure]] + delta)

    new_data_r[mediator] <- sample_z_inner

    task_m <- sl3::sl3_Task$new(
      data = new_data_m,
      covariates = covars,
      outcome = "y"
    )

    task_r <- sl3::sl3_Task$new(
      data = new_data_r,
      covariates = c(w_names),
      outcome = mediator,
    )

    m_val <- q_model$predict(task_m)
    r_val <- r_model$predict(task_r)

    output <- m_val * r_val$likelihood
    return(output)
  }

  integrand_m_g_r_mc <- function(sample_a, sample_z, row_data, covars, w_names, q_model, g_model, r_model, exposure, mediator, delta, upper_a) {
    n_samples <- length(sample_z)
    row_data_rep <- do.call("rbind", replicate(n_samples, row_data, simplify = FALSE))
    new_data_m <- new_data_g <- new_data_r <- row_data_rep

    new_data_m[exposure] <- ifelse(sample_a + delta >= upper_a, upper_a, sample_a + delta)
    new_data_m[mediator] <- sample_z

    new_data_g[exposure] <- sample_a
    new_data_r[mediator] <- sample_z

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

    task_r <- sl3::sl3_Task$new(
      data = new_data_r,
      covariates = c(w_names),
      outcome = mediator,
    )

    m_val <- q_model$predict(task_m)
    g_val <- g_model$predict(task_g)
    r_val <- r_model$predict(task_r)

    g_vals <- unlist(g_val, recursive = FALSE)
    likelihood <- sapply(g_vals, function(x) x[[sample_a]])
    output <- m_val * likelihood * r_val$likelihood

    return(output)
  }

  results <- numeric(nrow(av))
  integral_inner_results <- numeric(nrow(av))
  integral_outer_results <- numeric(nrow(av))

  for (i in 1:nrow(av)) {
    row_data <- av[i, ]
    # integral_inner_results_i <- numeric(n_bins)
    integral_outer_results_i <- numeric(n_bins)

    sample_z_inner <- runif(n_samples, lower_z, upper_z)
    mc_integrands_inner <- integrand_m_r(sample_z_inner, row_data, covars, w_names, q_model, r_model, exposure, mediator, delta, upper_a = n_bins)
    integral_inner <- (upper_z - lower_z) * mean(mc_integrands_inner)

    for (a_val in 1:n_bins) {
      sample_z_outer <- runif(n_samples, lower_z, upper_z)
      integrand_val <- integrand_m_g_r_mc(sample_a = a_val, sample_z = sample_z_outer, row_data, covars, w_names, q_model, g_model, r_model, exposure, mediator, delta, upper_z)

      integral_outer <- (upper_z - lower_z) * mean(integrand_val)

      # integral_inner_results_i[a_val] <- integral_inner
      integral_outer_results_i[a_val] <- integral_outer
    }

    results[i] <- integral_inner - sum(integral_outer_results_i)
    integral_inner_results[i] <- integral_inner
    integral_outer_results[i] <- integral_outer
  }

  return(list("d_a" = results, "phi_aw" = integral_inner_results))
}
