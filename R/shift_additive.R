#' Simple Additive Modified Treatment Policy
#'
#' @details A simple modified treatment policy that modifes the observed value
#'  of the exposure by shifting it by a value \code{delta}. Note that this
#'  shifting function assumes support of A|W across all strata of W.
#'
#' @param a A \code{numeric} vector of observed treatment values.
#' @param w A \code{numeric} matrix of observed baseline covariate values.
#' @param delta A \code{numeric} indicating the magnitude of the shift to be
#'  computed for the treatment \code{A}.
#'
#' @return A \code{numeric} vector containing the shifted exposure values.
shift_additive <- function(a, w = NULL, delta) {
  shifted_treatment <- a + delta
  return(shifted_treatment)
}
