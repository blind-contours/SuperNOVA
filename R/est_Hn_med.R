#' Estimate Auxiliary Covariate of Full Data Efficient Influence Function
#' for the Mediation Parameter
#'
#' @details Compute an estimate of the auxiliary covariate of the efficient
#'  influence function required to update initial estimates through logistic
#'  tilting models for targeted minimum loss estimation.
#'
#' @param gn_exp An estimate of the exposure density (a generalized propensity
#'  score) using the output provided by \code{\link{est_g_exp}}.
#'
#' @importFrom data.table as.data.table setnames
#'
#' @return A \code{data.table} with two columns, containing estimates of the
#'  auxiliary covariate at the natural value of the exposure H(A, W) and at the
#'  shifted value of the exposure H(A + delta, W).
#' @export
est_hn_med <- function(gn_exp, zn_exp) {
  # set any g(a|w) = 0 values to a very small value above zero
  gn_exp$noshift <- bound_propensity(gn_exp$noshift)
  zn_exp$noshift <- bound_propensity(zn_exp$noshift)

  # compute the ratio of the propensity scores for Hn(A,W)
  ratio_g_noshift <- (gn_exp$downshift / gn_exp$noshift) +
    as.numeric(gn_exp$upshift == 0)

  # compute the ratio of the propensity scores for Hn(d(A,W),W)
  ratio_g_shift <- (gn_exp$noshift / gn_exp$upshift) *
    as.numeric(gn_exp$upshift != 0) + as.numeric(gn_exp$upupshift == 0)

  ratio_g_shift[is.na(ratio_g_shift)] <- 1

  # set up output table of auxiliary covariate under different shifts
  aux_covar <- data.table::as.data.table(cbind(ratio_g_noshift, ratio_g_shift))
  data.table::setnames(aux_covar, c("noshift", "shift"))
  return(aux_covar)
}
