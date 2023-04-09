#' Integrate Psi G for the Da part of the EIF for stochastic mediation
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
#'  upshift of the exposure Q(A + delta, W).

integrate_psi_g <- function(data, covars, w_names, q_model, r_model, g_model, exposure, mediator, delta) {
  data <- as.data.frame(data)
  # data[[exposure]] <- round(data[[exposure]])
  # data[[mediator]] <- round(data[[mediator]])

  lower_z <- min(data[[mediator]])
  upper_z <- max(data[[mediator]])

  lower_a <- min(data[[exposure]])
  upper_a <- max(data[[exposure]])


  integrand_m_r <- function(z, row_data, covars, w_names, q_model, r_model, exposure, mediator, delta, upper_a) {
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


  integrand_m_g_r <- function(a, z, row_data, covars, w_names, q_model, g_model, r_model, exposure, mediator, delta, upper_a) {
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

      output <- m_val * g_val$likelihood * r_val$likelihood
      output_matrix[, col] <- output
    }


    return(output_matrix)
  }

  results <- numeric(nrow(data))
  integral_inner_results <- numeric(nrow(data))

  for (i in 1:nrow(data)) {
    row_data <- data[i, ]

    integral_inner <- stats::integrate(
      function(z) integrand_m_r(z, row_data, covars, w_names, q_model, r_model, exposure, mediator, delta, upper_a),
      lower = lower_z,
      upper = upper_z,
      rel.tol = 0.0001,
      subdivisions = 100,
      stop.on.error = FALSE
    )$value


    integral_outer <- integral2(
      function(a, z) integrand_m_g_r(a, z, row_data, covars, w_names, q_model, g_model, r_model, exposure, mediator, delta, upper_a),
      xmin = lower_a,
      ymin = lower_z,
      xmax = upper_a,
      ymax = upper_z,
    )$Q

    results[i] <- integral_inner - integral_outer
    integral_inner_results[i] <- integral_inner
  }

  return(list("d_a" = results, "phi_aw" = integral_inner_results))
}
