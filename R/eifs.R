#' Compute the Shift Parameter Estimate and the Efficient Influence Function
#'
#' @details Estimate the value of the causal parameter alongside statistical
#'  inference for the parameter estimate based on the efficient influence
#'  function of the target parameter.
#'
#' @param y A \code{numeric} vector of the observed outcomes.
#' @param qn An object providing the value of the outcome evaluated after
#'  imposing a shift in the treatment. This object is passed in after being
#'  constructed by a call to the internal function \code{est_Q}.
#' @param hn An object providing values of the auxiliary ("clever") covariate,
#'  constructed from the treatment mechanism and required for targeted minimum
#'  loss-based estimation. This object object should be passed in after being
#'  constructed by a call to the internal function \code{est_Hn}.
#' @param estimator The type of estimator to be fit, either \code{"tmle"} for
#'  targeted maximum likelihood estimation or \code{"onestep"} for a one-step
#'  estimator.
#' @param fluc_mod_out An object giving values of the logistic tilting model
#'  for targeted minimum loss estimation. This type of object should be the
#'  output of the internal routines to perform this step of the TML estimation
#'  procedure, as given by \code{\link{fit_fluctuation}}.
#'
#' @importFrom stats var
#'
#' @return A \code{list} containing the parameter estimate, estimated variance
#'  based on the efficient influence function (EIF), the estimate of the EIF
#'  incorporating inverse probability of censoring weights, and the estimate of
#'  the EIF without the application of such weights.
eif <- function(y,
                qn,
                qn_unscaled,
                hn,
                estimator = c("tmle", "onestep"),
                fluc_mod_out = NULL,
                data) {
  # set TMLE as default estimator type
  estimator <- match.arg(estimator)

  # set Qn to use based on estimator type
  if (estimator == "tmle") {
    qn_shift <- fluc_mod_out$qn_shift_star
    qn_noshift <- fluc_mod_out$qn_noshift_star
  } else if (estimator == "onestep") {
    qn_shift <- qn_unscaled$upshift
    qn_noshift <- qn_unscaled$noshift
  }

  psi <- mean(qn_shift)
  noshift_psi <- mean(y)

  # compute the efficient influence function (EIF)
  eif <- rep(0, length(qn_shift))
  eif_no_shift <- rep(0, length(qn_shift))

  eif <-
    (hn$noshift * (y - qn_noshift) + (qn_shift - psi))

  eif_no_shift <- y - qn_noshift


  # add mean of EIF to parameter estimate if fitting one-step
  # NOTE: the estimate of psi is updated _after_ evaluating the EIF
  if (estimator == "onestep") {
    psi <- psi + mean(eif)
  }

  # compute the variance based on the EIF and scale by number of observations
  var_eif <- stats::var(eif) / length(y)
  var_noshift_eif <- stats::var(eif_no_shift) / length(y)

  # compute the confidence intervales based on the EIF
  se <- sqrt(var_eif)
  se_noshift <- sqrt(var_noshift_eif)

  CI <- c(
    round(psi + stats::qnorm(0.05 / 2, lower.tail = T) * se, 4),
    round(psi + stats::qnorm(0.05 / 2, lower.tail = F) * se, 4)
  )

  CI_noshift <- c(
    round(noshift_psi + stats::qnorm(0.05 / 2, lower.tail = T) * se_noshift, 4),
    round(noshift_psi + stats::qnorm(0.05 / 2, lower.tail = F) * se_noshift, 4)
  )

  p.value <- round(2 * stats::pnorm(abs(psi / se), lower.tail = F), 6)

  # return the variance and the EIF value at each observation
  out <- list(
    psi = psi,
    var = var_eif,
    se <- se,
    CI <- CI,
    p_value <- p.value,
    eif = eif,
    noshift_psi = noshift_psi,
    noshift_var = var_noshift_eif,
    noshift_se = se_noshift,
    noshift_CI = CI_noshift,
    noshift_eif = eif_no_shift
  )

  names(out) <- c(
    "psi", "var", "se", "CI", "p_value", "eif",
    "no shift psi", "no shift var", "no shift se", "no shift CI", "no shift eif"
  )
  return(out)
}
