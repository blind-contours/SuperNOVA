#' DGP for testing SuperNOVA
#'
#' @param n_obs Number of observations
#' @param mu Vector of means for exposures
#' @param sigma_mod Sigma matrix of exposures
#' @param effect_mod_binary TRUE/FALSE Whether the effect modifier is binary
#' @param delta_1 Delta for M1 Exposure
#' @param delta_2 Delta for M3 Exposure
#'
#' @return
#' @export
#' @importFrom MASS mvrnorm
#' @importFrom bindata rmvbin

#' @examples
simulate_data <- function(n_obs = 100000,
                          shift_var_index = 1,
                          sigma_mod = matrix(
                            c(
                              1, 0.5, 0.8, 0.5,
                              1, 0.7, 0.8, 0.7, 1
                            ),
                            nrow = 3,
                            ncol = 3
                          ),
                          delta = 1) {

  covars <- MASS::mvrnorm(
    n = n_obs, mu = c(6, 7),
    Sigma = matrix(c(1, 0.4, 0.4, 1), 2, 2)
  )

  mu_1 <- mean(exp(covars[,1] / 4))
  mu_2 <- mean(covars[,2]/5)
  mu_3 <- 2

  mixtures <- MASS::mvrnorm(n = n_obs,
                            mu = c(mu_1, mu_2, mu_3),
                            Sigma = sigma_mod)


  y_mean = function(mixtures, covars) {
    0.1 * log(4 * pi * mixtures[, 1]) + 0.6 * plogis(2 * mixtures[, 3]) *
      covars[, 1] + 0.4 * (plogis(2*mixtures[, 3]) * exp(mixtures[, 1]) / 100) + 0.1 *
      covars[, 1] + 0.3 * covars[, 2]
  }

  y <- y_mean(mixtures, covars)

  mixtures[, shift_var_index] <- mixtures[, shift_var_index]+delta

  y_shifted <- y_mean(mixtures, covars)

  effect = mean(
    y_shifted -
      y
  )


  data <- as.data.frame(cbind(mixtures, covars, y, y_shifted))
  colnames(data) <- c("M1", "M2", "M3", "W1", "W2", "Y", "Y_shifted")

  return(list("data" = data, "effect" = effect))
}
