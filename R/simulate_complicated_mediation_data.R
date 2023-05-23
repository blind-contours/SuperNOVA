#' DGP for testing SuperNOVA with mediation (complicated!)
#'
#' @param n_obs Number of observations
#' @param delta Delta for M1 Exposure
#' @param breaks Number of bins to discretize exposure
#'
#' @return A data frame of simulated mediation data
#' @export
#' @importFrom MASS mvrnorm
#' @importFrom stats rnorm rbinom rpois

simulate_complicated_mediation_data <- function(n_obs = 100000,
                                    delta = 1,
                                    n_bins = 10)  {
  # simulate some baseline covariates
  w_1 <- rnorm(n_obs, mean = 20, 2)
  w_2 <- rbinom(n_obs, 1, 0.5)
  w_3 <- rbinom(n_obs, 1, 0.5)
  w_4 <- rnorm(n_obs, 30, 3)
  w_5 <- rpois(n_obs, 1.2)

  # generate exposures and discretize them
  a_1 <- rnorm(n_obs, 1 + 0.5 * w_1, 1)
  a_1_breaks <- quantile(a_1, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE)
  a_1_quant <- as.numeric(cut(a_1, breaks = a_1_breaks, labels = FALSE, include.lowest = TRUE))
  a_1_quant_shift <- ifelse(a_1_quant + delta <= max(a_1_quant), a_1_quant + delta, max(a_1_quant))

  a_2 <- rnorm(n_obs, mean = 2 * w_2 * w_3, sd = 1)
  a_2_breaks <- quantile(a_2, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE)
  a_2_quant <- as.numeric(cut(a_2, breaks = a_2_breaks, labels = FALSE, include.lowest = TRUE))
  a_2_quant_shift <- ifelse(a_2_quant + delta <= max(a_2_quant), a_2_quant + delta, max(a_2_quant))


  a_3 <- rnorm(n_obs, mean = 1.5 * w_4 / 20 * w_1 / 3, sd = 2)
  a_3_breaks <- quantile(a_3, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE)
  a_3_quant <- as.numeric(cut(a_3, breaks = a_3_breaks, labels = FALSE, include.lowest = TRUE))
  a_3_quant_shift <- ifelse(a_3_quant + delta <= max(a_3_quant), a_3_quant + delta, max(a_3_quant))

  a_4 <- rnorm(n_obs, mean = 3 * w_4 / 2 * w_2 / 3, sd = 3)
  a_4_breaks <- quantile(a_4, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE)
  a_4_quant <- as.numeric(cut(a_4, breaks = a_4_breaks, labels = FALSE, include.lowest = TRUE))
  a_4_quant_shift <- ifelse(a_4_quant + delta <= max(a_4_quant), a_4_quant + delta, max(a_4_quant))

  # generate z from discrete a and w and discrete a shifted z given a shift in a

  z_1_quant <- rnorm(n_obs, 2 * a_1_quant + w_1, 1)
  z_1_shift_quant <- rnorm(n_obs, 2 * a_1_quant_shift + w_1, 1)

  z_2_quant <- rnorm(n_obs, 2 * a_2_quant + w_2, 1)
  z_2_shift_quant <- rnorm(n_obs, 2 * a_2_quant_shift + w_2, 1)

  z_3_quant <- rnorm(n_obs, 5 * a_3_quant * a_4 + w_3, 1)
  z_3_shift_quant <- rnorm(n_obs, 5 * a_3_quant_shift * a_4_quant_shift + w_1, 1)

  # y <- 10 * z_1 + 40 * a_1 + 15 - w_3 + a_2 + z_2
  y_quant <- 10 * z_1_quant + 40 * a_1_quant + 15 - w_3  + 6 * a_2_quant + 7 * z_2_quant

  # y_shift_a_1 <- 10 * z_1 + 40 * a_1_shift + 15 - w_3 + a_2 + z_2
  y_shift_a_1_shift <- 10 * z_1_quant + 40 * a_1_quant_shift + 15 - w_3  + 6 * a_2_quant + 7 * z_2_quant

  # y_shift_a_1_z_1 <- 10 * z_1_shift + 40 * a_1_shift + 15 - w_3 + a_2 + z_2
  y_shift_a_1_z_1_shift <- 10 * z_1_shift_quant + 40 * a_1_quant_shift + 15 - w_3  + 6 *  a_2_quant + 7 * z_2_quant

  # y_shift_a_2 <- 10 * z_1 + 40 * a_1 + 15 - w_3 + a_2_shift + z_2
  y_shift_a_2_shift <- 10 * z_1_quant + 40 * a_1_quant + 15 - w_3 + 6 * a_2_quant_shift + 7 *  z_2_quant

  # y_shift_a_2_z_2 <- 10 * z_1 + 40 * a_1 + 15 - w_3 + a_2_shift + z_2_shift
  y_shift_a_2_z_2_shift <- 10 * z_1_quant + 40 * a_1_quant + 15 - w_3 + 6 * a_2_quant_shift + 7 * z_2_shift_quant

  nde_a1 <- mean(y_shift_a_1_shift - y_quant)
  nie_a1 <- mean(y_shift_a_1_z_1_shift - y_shift_a_1_shift)
  ate_a1 <- nde_a1 + nie_a1

  nde_a2 <- mean(y_shift_a_2_shift - y_quant)
  nie_a2 <- mean(y_shift_a_2_z_2_shift - y_shift_a_2_shift)
  ate_a2 <- nde_a2 + nie_a2

  data <- as.data.frame(cbind(
    w_1, w_2, w_3, w_4, w_5,
    a_1, a_2, a_3, a_4,
    a_1_quant, a_1_quant_shift,
    a_2_quant, a_2_quant_shift,
    a_3_quant, a_3_quant_shift,
    a_4_quant, a_4_quant_shift,
    z_1_quant, z_1_shift_quant,
    z_2_quant, z_2_shift_quant,
    z_3_quant, z_3_shift_quant,
    y_quant,
    y_shift_a_1_shift, y_shift_a_1_z_1_shift,
    y_shift_a_2_shift, y_shift_a_2_z_2_shift
  ))

  glm_model <- glm(as.formula(y_quant~z_1_quant + a_1 + w_3 + a_2_quant + z_2_quant), data = data)
  y_quant <- predict.glm(glm_model, newdata = as.data.frame(cbind(z_1_quant, a_1_quant, w_3, a_2_quant, z_2_quant)))
  y_shift_a_1_shift <- predict.glm(glm_model, newdata = as.data.frame(cbind(z_1_quant, a_1_quant_shift, w_3, a_2_quant, z_2_quant)))
  y_shift_a_1_z_1_shift <- predict.glm(glm_model, newdata = as.data.frame(cbind(z_1_shift_quant, a_1_quant_shift, w_3, a_2_quant, z_2_quant)))


  return(list(
    "data" = data,
    "nde_a1" = nde_a1,
    "nie_a1" = nie_a1,
    "ate_a1" = ate_a1,
    "nde_a2" = nde_a2,
    "nie_a2" = nie_a2,
    "ate_a2" = ate_a2
  ))
}
