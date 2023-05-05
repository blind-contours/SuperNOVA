#' Joint Shift Additive Modified Treatment Policy
#'
#' @details A simple modified treatment policy that modifes the observed value
#'  of the exposure by shifting it by a value \code{delta}. Note that this
#'  shifting function assumes support of A|W across all strata of W. Same
#'  as shift additive but for two exposures found in the data-adaptive procedure
#'
#' @param a1 A \code{numeric} vector of observed treatment values.
#' @param a2 A \code{numeric} vector of observed treatment values.=
#' @param w A \code{numeric} matrix of observed baseline covariate values.
#' @param delta1 A \code{numeric} indicating the magnitude of the shift to be
#'  computed for the treatment \code{A}.
#' @param delta2 A \code{numeric} indicating the magnitude of the shift to be
#'  computed for the treatment \code{A}.
#' @param lower_bound1 Lower bound of expossure
#' @param lower_bound2 Lower bound of expossure
#' @param upper_bound1 Upper bound of exposure
#' @param upper_bound2 Upper bound of exposure

#'
#' @return A \code{numeric} vector containing the shifted exposure values.
joint_shift_additive <- function(a1, a2, w = NULL, delta1, delta2, lower_bound1, upper_bound1, lower_bound2, upper_bound2) {
  shifted_treatment1 <- ifelse(a1 + delta1 <= upper_bound1 & a1 + delta1 >= lower_bound1,
                               a1 + delta1,
                               ifelse(a1 + delta1 < lower_bound1, lower_bound1, upper_bound1)
  )

  shifted_treatment2 <- ifelse(a2 + delta2 <= upper_bound2 & a2 + delta2 >= lower_bound2,
                               a2 + delta2,
                               ifelse(a2 + delta2 < lower_bound2, lower_bound2, upper_bound2)
  )

  return(list(shifted_treatment1, shifted_treatment2))
}
