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
                                    delta = 1)  {
  # simulate some baseline covariates
  w_1 <- rnorm(n_obs, mean = 20, sd = 2)
  w_2 <- rbinom(n_obs, 1, 0.5)
  w_3 <- rbinom(n_obs, 1, 0.5)
  w_4 <- rnorm(n_obs, mean = 30, sd = 3)
  w_5 <- rpois(n_obs, 1.2)

  # mean and covariance matrix for A variables
  mu <- c(mean(1 + 0.5 * w_1),
          mean(2 * w_2 * w_3),
          mean(1.5 * w_4 / 20 * w_1 / 3),
          mean(3 * w_4 / 2 * w_2 / 3),
          mean(2 * w_5))

  Sigma <- matrix(c(1, 0.8, 0.3, 0.3, 0.2,
                    0.8, 1, 0.3, 0.3, 0.2,
                    0.3, 0.3, 1, 0.8, 0.2,
                    0.3, 0.3, 0.8, 1, 0.2,
                    0.2, 0.2, 0.2, 0.2, 1), ncol = 5, byrow = TRUE)

  # generate exposures
  a <- MASS::mvrnorm(n_obs, mu, Sigma)
  a_shift <- pmin(a + delta, apply(a, 2, max)) # Apply the max for each column (exposure)

  # generate mediator variables and shifted versions
  z_1 <- rnorm(n_obs, 2 * a[,1] + w_1, 1)
  z_1_shift <- rnorm(n_obs, 2 * a_shift[,1] + w_1, 1)

  z_2 <- rnorm(n_obs, 2 * a[,2] + w_2, 1)
  z_2_shift <- rnorm(n_obs, 2 * a_shift[,2] + w_2, 1)

  z_3 <- rnorm(n_obs, 5 * a[,3] * a[,4] + w_3, 1)
  z_3_shift <- rnorm(n_obs, 5 * a_shift[,3] * a_shift[,4] + w_3, 1)

  z_4 <- rnorm(n_obs, 3 * a[,4] * w_4, 1)
  z_4_shift <- rnorm(n_obs, 3 * a_shift[,4] * w_4, 1)

  z_5 <- rnorm(n_obs, 4 * a[,5] * w_5, 1)
  z_5_shift <- rnorm(n_obs, 4 * a_shift[,5] * w_5, 1)

  z <- cbind(z_1, z_2, z_3, z_4, z_5)
  z_shift <- cbind(z_1_shift, z_2_shift, z_3_shift, z_4_shift, z_5_shift)

  # y and shifted versions
  y <- 10 * z[,1] + 40 * a[,1] + 15 - w_3  + 6 * a[,2] + 7 * z[,2]

  data <- as.data.frame(cbind(
    w_1, w_2, w_3, w_4, w_5,
    a,
    z,
    y
  ))

  names(data) <- c(paste0("w_", 1:5), paste0("a_", 1:5), paste0("z_", 1:5), "y")

  return(
   data
  )
}
