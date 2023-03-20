#' DGP for testing SuperNOVA without mediation
#'
#' @param n_obs Number of observations
#' @param sigma_mod Sigma matrix of exposures
#' @param delta Amount to shift exposure by
#' @param shift_var_index for the exposure vector that indicates which
#' exposure to shift
#'
#' @return A dataframe of simulated data
#' @export
#' @importFrom MASS mvrnorm

simulate_data <- function(n_obs = 100000,
                          sigma_mod = matrix(
                            c(
                              1, 0.5, 0.8, 0.5,
                              1, 0.7, 0.8, 0.7, 1
                            ),
                            nrow = 3,
                            ncol = 3
                          ),
                          delta = 0.5) {
  covars <- MASS::mvrnorm(
    n = n_obs, mu = c(6, 7),
    Sigma = matrix(c(1, 0.4, 0.4, 1), 2, 2)
  )

  covars <- as.data.frame(covars)
  covars[, 3] <- rbinom(n_obs, size = 1, prob = 0.5)

  mu_1 <- mean(exp(covars[, 1] / 2))
  mu_2 <- mean(covars[, 2] / 2)
  mu_3 <- 5

  mixtures <- MASS::mvrnorm(
    n = n_obs,
    mu = c(mu_1, mu_2, mu_3),
    Sigma = sigma_mod
  )

  mixtures <- as.data.frame(mixtures)
  mixtures[, 4] <- rnorm(n_obs, mean = 4, sd = 2)

  mixtures_m14 <- mixtures_m1 <- mixtures_m2 <- mixtures_m3 <- mixtures_m4 <- mixtures
  covars_0 <- covars_1 <- covars

  covars_0[, 3] <- 0
  covars_1[, 3] <- 1

  y_mean <- function(mixtures, covars) {
    1.3 * mixtures[, 4] +
      ifelse(covars[, 3] == 1, mixtures[, 3]^2, mixtures[, 3]) +
      0.4 * mixtures[, 4] * mixtures[, 1] +
      0.1 * covars[, 1] + 0.3 * covars[, 2]
  }

  y <- y_mean(mixtures, covars)

  # shift mixture 1
  mixtures_m1[, 1] <- mixtures[, 1] + delta
  m1_y_shifted <- y_mean(mixtures_m1, covars)

  m1_effect <- mean(
    m1_y_shifted -
      y
  )

  # shift mixture 2
  mixtures_m2[, 2] <- mixtures[, 2] + delta
  m2_y_shifted <- y_mean(mixtures_m2, covars)

  m2_effect <- mean(
    m2_y_shifted -
      y
  )

  # shift mixture 3
  mixtures_m3[, 3] <- mixtures[, 3] + delta
  m3_y_shifted <- y_mean(mixtures_m3, covars)

  diff <- m3_y_shifted - y

  m3_y_shifted_covar0 <- mean(diff[covars$V3 == 0])
  m3_y_shifted_covar1 <- mean(diff[covars$V3 == 1])

  effect_mod_results <- list(
    "Level 0 Shift Diff in W3 <= 0" = m3_y_shifted_covar1,
    "Level 1 Shift Diff in W3 <= 0" = m3_y_shifted_covar0
  )

  m3_effect <- mean(
    m3_y_shifted -
      y
  )

  # shift mixture 4
  mixtures_m4[, 4] <- mixtures[, 4] + delta
  m4_y_shifted <- y_mean(mixtures_m4, covars)

  m4_effect <- mean(
    m4_y_shifted -
      y
  )

  # shift mixture 1 and 2
  mixtures_m14[, c(1, 4)] <- mixtures[, c(1, 4)] + delta
  m14_y_shifted <- y_mean(mixtures_m14, covars)

  m14_effect <- mean(
    m14_y_shifted -
      y
  )

  m14_intxn <- m14_effect - (m1_effect + m4_effect)

  data <- as.data.frame(cbind(mixtures, covars, y))
  colnames(data) <- c("M1", "M2", "M3", "M4", "W1", "W2", "W3", "Y")


  return(list(
    "data" = data,
    "m1_effect" = m1_effect,
    "m2_effect" = m2_effect,
    "m3_effect" = m3_effect,
    "m4_effect" = m4_effect,
    "m14_effect" = m14_effect,
    "m14_intxn" = m14_intxn,
    "effect_mod" = effect_mod_results
  ))
}
