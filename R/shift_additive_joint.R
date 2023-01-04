#' Simple Additive Modified Treatment Policy for Two Exposures
#'
#' @details A simple modified treatment policy that modifes the observed value
#'  of the exposure by shifting it by a value \code{delta}. Note that this
#'  shifting function assumes support of A|W across all strata of W. This function
#'  shifts two treatments or exposures.
#'
#' @param A A \code{numeric} vector of observed treatment values.
#' @param W A \code{numeric} matrix of observed baseline covariate values.
#' @param deltas A \code{numeric} list indicating the magnitude of the shift to be
#'  computed for the treatment \code{A}.
#'
#' @return A \code{numeric} vector containing the shifted exposure values.
shift_additive_joint <- function(a, w = NULL, deltas) {
  shifted_treatment_1 <- a[, 1] + deltas[[1]]
  shifted_treatment_2 <- a[, 2] + deltas[[2]]
  shifted_treatment <- cbind(shifted_treatment_1, shifted_treatment_2)
  return(shifted_treatment)
}
