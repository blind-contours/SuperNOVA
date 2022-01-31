

restrict_intxn_earth <- function(degree, pred, parents, namesx) {
  if (degree > 1) {
    predictor <- namesx[pred]
    parents <- namesx[parents != 0]

    if ((grepl("W", predictor) && any(grepl("M|V", parents))) ||
      (grepl("W", parents) && any(grepl("M|V", predictor)))) {
      return(FALSE)
    }
  }
  return(TRUE)
}

restrict_intxn_polymars <- function(data, W_names, exposure_names) {
  w_indices <- match(W_names, colnames(data))
  exposure_indices <- match(exposure_names, colnames(data))
  no_intxn_combns <- expand.grid(exposure_indices, w_indices)

  return(no_intxn_combns)
}

no.int <- function(degree, pred, parents, namesx) {
  if (degree > 1) {
    predictor <- namesx[pred]
    parents <- namesx[parents != 0]
    if ((any(predictor %in% PREDICTORS) && any(parents %in% PARENTS)) ||
      (any(predictor %in% PARENTS) && any(parents %in% PREDICTORS))) {
      return(FALSE)
    }
  }
  TRUE
}
