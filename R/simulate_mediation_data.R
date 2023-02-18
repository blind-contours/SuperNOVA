#' DGP for testing SuperNOVA with mediation
#'
#' @param n_obs Number of observations
#' @param delta Delta for M1 Exposure
#'
#' @return A data frame of simulated mediation data
#' @export
#' @importFrom MASS mvrnorm
#' @importFrom stats rnorm rbinom

simulate_mediation_data <- function(n_obs = 100000,
                                    delta = 1) {
  # simulate some baseline covariates
  w_1 <- rnorm(n_obs, mean = 20, 2)
  w_2 <- rbinom(n_obs, 1, 0.5)
  w_3 <- rbinom(n_obs, 1, 0.5)

  a_1 <- rnorm(n_obs, 1 + 0.5 * w_1, 1)
  a_1_shift <- a_1 + delta

  a_2 <- rnorm(n_obs, mean = 2 * w_2 * w_3, sd = 1)
  a_2_shift <- a_2 + delta

  z_1 <- rnorm(n_obs, 2 * a_1 + w_1, 1)
  z_1_shift <- rnorm(n_obs, 2 * a_1_shift + w_1, 1)

  y <- (a_1^2 * z_1) + ifelse(w_2 == 1, a_2 * 400, a_2 * 2) - w_3
  y_shift_a_1 <- (a_1_shift^2 * z_1) + ifelse(w_2 == 1, a_2 * 400, a_2 * 2) - w_3
  y_shift_a_2 <- (a_1^2 * z_1) + ifelse(w_2 == 1, a_2_shift * 400, a_2_shift * 2) - w_3
  y_shift_a_1_z_1 <- (a_1_shift^2 * z_1_shift) + ifelse(w_2 == 1, a_2 * 400, a_2 * 2) - w_3

  nde <- mean(y_shift_a_1 - y)
  nie <- mean(y_shift_a_1_z_1 - y_shift_a_1)
  ate_a_1 <- nde + nie

  ate_a_2 <- mean(y_shift_a_2 - y)

  data <- as.data.frame(cbind(
    w_1, w_2, w_3, a_1, a_1_shift, a_2, a_2_shift, z_1, z_1_shift, y,
    y_shift_a_1, y_shift_a_2, y_shift_a_1_z_1
  ))

  y_shift_a_2_w_2_1 <- mean(subset(data, w_2 == 1)$y_shift_a_2)
  y_shift_a_2_w_2_0 <- mean(subset(data, w_2 == 0)$y_shift_a_2)

  return(list(
    "data" = data,
    "nde" = nde,
    "nie" = nie,
    "ate_a_1" = ate_a_1,
    "ate_a_2" = ate_a_2,
    "effect_mod_a2w2_lvl_1" = y_shift_a_2_w_2_1,
    "effect_mod_a2w2_lvl_0" = y_shift_a_2_w_2_0
  ))
}
