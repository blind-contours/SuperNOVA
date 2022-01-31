#' Title
#'
#' @param n_obs
#' @param mu
#' @param sigma_mod
#'
#' @return
#' @export
#' @importFrom MASS mvrnorm
#' @importFrom bindata rmvbin

#' @examples
simulate_data <- function(n_obs = 2000,
                          mu = c(2.5, 3, 4),
                          sigma_mod = matrix(c(1, 0.5, 0.8, 0.5, 1, 0.7, 0.8, 0.7, 1),
                            nrow = 3,
                            ncol = 3
                          )) {
  mixtures <- MASS::mvrnorm(n = n_obs, mu = mu, Sigma = sigma_mod)

  mixtures[mixtures < 0] <- 0.001

  ## Construct a binary correlation matrix
  # m <- matrix(c(1,rho,rho,1), ncol=2)
  #
  # ## Simulate 10000 x-y pairs, and check that they have the specified
  # ## correlation structure
  # w <- bindata::rmvbin(n_obs, margprob = c(0.5, 0.5), bincorr = m)
  # cor(w)

  w_em <- rbinom(n_obs, 1, 0.48)
  covars <- MASS::mvrnorm(n = n_obs, mu = c(6, 7), Sigma = matrix(c(1, 0.4, 0.4, 1), 2, 2))

  y <- 0.2 * log(2 * pi * mixtures[, 1]) + 0.6 * (((exp(mixtures[, 1]) / 100) * plogis(mixtures[, 3])) * w_em) + 0.4 * (plogis(mixtures[, 3])) + 0.1 * (covars[, 1] * covars[, 2]) + rnorm(dim(mixtures)[1], mean = 0, sd = 0.01)

  data <- as.data.frame(cbind(mixtures, w_em, covars, y))
  colnames(data) <- c("M1", "M2", "M3", "W1", "W2", "W3", "Y")

  m1_w1_plot <- ggplot2::ggplot(data = data, ggplot2::aes(x = data$M1, y = data$Y, group = data$W1, color = as.factor(data$W1))) +
    ggplot2::geom_point() +
    ggplot2::xlab("Mixture 1") +
    ggplot2::ylab("Outcome") +
    ggplot2::labs(color = "Covariate 1")

  m2_w1_plot <- ggplot2::ggplot(data = data, ggplot2::aes(x = data$M2, y = data$Y, group = data$W1, color = as.factor(data$W1))) +
    ggplot2::geom_point() +
    ggplot2::xlab("Mixture 3") +
    ggplot2::ylab("Outcome") +
    ggplot2::labs(color = "Covariate 1")

  m3_w1_plot <- ggplot2::ggplot(data = data, ggplot2::aes(x = data$M3, y = data$Y, group = data$W1, color = as.factor(data$W1))) +
    ggplot2::geom_point() +
    ggplot2::xlab("Mixture 3") +
    ggplot2::ylab("Outcome") +
    ggplot2::labs(color = "Covariate 1")


  return(list("data" = data, "plot 1" = m1_w1_plot, "plot 2" = m2_w1_plot, "plot 3" = m3_w1_plot))
}
