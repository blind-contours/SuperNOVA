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
  w_4 <- rnorm(n_obs, 30, 3)
  w_5 <- rpois(n_obs, 1.2)

  a_1 <- rnorm(n_obs, 1 + 0.5 * w_1, 1)
  a_1_quant <- cut(a_1, quantile(a_1), labels = FALSE, include.lowest = TRUE)

  a_1_shift <- ifelse(a_1 + delta <= max(a_1), a_1 + delta, max(a_1))
  a_1_quant_shift <- ifelse(a_1_quant + delta <= max(a_1_quant), a_1_quant + delta, max(a_1_quant))

  # a_2 <- rnorm(n_obs, mean = 2 * w_2 * w_3, sd = 1)
  # a_2_shift <- a_2 + delta
  #
  # a_3 <- rnorm(n_obs, mean = 1.5 * w_4 / 20 * w_1 / 3, sd = 2)
  # a_3_shift <- a_3 + delta
  #
  # a_4 <- rnorm(n_obs, mean = 2 * w_2 * w_5, sd = 1)
  # a_4_shift <- a_4 + delta

  z_1 <- rnorm(n_obs, 2 * a_1 + w_1, 1)
  z_1_quant <- rnorm(n_obs, 2 * a_1_quant + w_1, 1)

  z_1_shift <- rnorm(n_obs, 2 * a_1_shift + w_1, 1)
  z_1_shift_quant <- rnorm(n_obs, 2 * a_1_quant_shift + w_1, 1)

  # z_2 <- rnorm(n_obs, 2 * a_2 + w_2, 1)
  # z_2_shift <- rnorm(n_obs, 2 * a_2_shift + w_1, 1)
  #
  # z_3 <- rnorm(n_obs, 5 * a_3 * a_4 + w_3, 1)
  # z_3_shift <- rnorm(n_obs, 5 * a_3_shift * a_4_shift + w_1, 1)

  y <- 10 * z_1 + 40 * a_1 #+ z_3 + 15 * a_3 * a_4 - w_3
  y_quant <- 10 * z_1_quant + 40 * a_1_quant #+ z_3 + 15 * a_3 * a_4 - w_3

  y_shift_a_1 <- 10 * z_1 + 40 * a_1_shift #+ z_3 + 15 * a_3 * a_4 - w_3
  y_shift_a_1_quant <- 10 * z_1_quant + 40 * a_1_quant_shift #+ z_3 + 15 * a_3 * a_4 - w_3

  # y_shift_a_3 <- 10 * z_1 + 40 * a_1 #+ z_3 + 15 * a_3_shift * a_4 - w_3
  # y_shift_a_4 <- 10 * z_1 + 40 * a_1 #+z_3 + 15 * a_3 * a_4_shift - w_3


  y_shift_a_1_z_1 <- 10 * z_1_shift + 40 * a_1_shift #+ z_3 + 15 * a_3 * a_4 - w_3
  y_shift_a_1_z_1_quant <- 10 * z_1_shift_quant + 40 * a_1_quant_shift #+ z_3 + 15 * a_3 * a_4 - w_3
  # y_shift_a_3_z_3 <- 10 * z_1 + 40 * a_1 #+ z_3_shift + 15 * a_3_shift * a_4 - w_3
  # y_shift_a_4_z_3 <- 10 * z_1 + 40 * a_1 #+ z_3_shift + 15 * a_3 * a_4_shift - w_3

  # y_shift_a_34 <- 10 * z_1 + 40 * a_1 #+ z_3 + 15 * a_3_shift * a_4_shift - w_3
  # y_shift_a_34_z_3 <- 10 * z_1 + 40 * a_1 #+ z_3_shift + 15 * a_3_shift * a_4_shift - w_3

  nde_a1 <- mean(y_shift_a_1 - y)
  nie_a1 <- mean(y_shift_a_1_z_1 - y_shift_a_1)
  ate_a1 <- nde_a1 + nie_a1

  nde_a1_quant <- mean(y_shift_a_1_quant - y_quant)
  nie_a1_quant <- mean(y_shift_a_1_z_1_quant - y_shift_a_1_quant)
  ate_a1_quant <- nde_a1_quant + nie_a1_quant

  # nde_a3 <- mean(y_shift_a_3 - y)
  # nie_a3 <- mean(y_shift_a_3_z_3 - y_shift_a_3)
  # ate_a3 <- nde_a3 + nie_a3
  #
  # nde_a4 <- mean(y_shift_a_4 - y)
  # nie_a4 <- mean(y_shift_a_4_z_3 - y_shift_a_4)
  # ate_a4 <- nde_a4 + nie_a4
  #
  # nde_a3a4 <- mean(y_shift_a_34 - y)
  # nie_a3a4 <- mean(y_shift_a_34_z_3 - y_shift_a_34)
  # ate_a34 <- nde_a3a4 + nie_a3a4

  data <- as.data.frame(cbind(
    w_1, w_2, w_3, w_4, w_5,
    a_1, a_1_shift, a_1_quant, a_1_quant_shift,
    z_1, z_1_shift, z_1_quant, z_1_shift_quant,
    y, y_quant, y_shift_a_1, y_shift_a_1_quant, y_shift_a_1_z_1, y_shift_a_1_z_1_quant
  ))


  return(list(
    "data" = data,
    "nde_a1" = nde_a1,
    "nie_a1" = nie_a1,
    "ate_a1" = ate_a1,
    "nde_a1_quant" = nde_a1_quant,
    "nie_a1_quant" = nie_a1_quant,
    "ate_a1_quant" = ate_a1_quant
  ))
}
