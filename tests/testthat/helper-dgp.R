# data generating mechanism function factory
make_dgp <- function() {
  # treatment mechanism
  g_mech <- function(w) {
    p_score <- (rowSums(w) / 4) + 0.1
    return(p_score)
  }

  # mediation mechanism
  z_mech <- function(w, a) {
    w <- data.frame(w)
    z1_prob <- 1 - plogis((a + w[, 1]) / (a + w[, 1]^3 + 0.5))
    z2_prob <- plogis((a - 1) + w[, 2] / (w[, 3] + 3))
    z3_prob <- plogis((a - 1) + 2 * w[, 1]^3 - 1 / (2 * w[, 1] + 0.5))
    return(list(z1_prob, z2_prob, z3_prob))
  }

  # outcome mechanism for continuous
  m_mech_cont <- function(w, a, z, eps_sd = 0.05) {
    w <- data.frame(w)
    z <- data.frame(z)
    y <- ifelse(w[1] == 1, (z[, 3]^2 + exp(a + z[, 3] / (1 + rowSums(w)^2))),
      (z[, 2] + exp(a + z[, 2] / (1 + rowSums(w)^2)))
    )

    y <- y + rnorm(length(a), mean = 0, sd = eps_sd)

    return(y)
  }


  # outcome mechanism for binary
  m_mech_binary <- function(w, a, z_probs, eps_sd = 0.5) {
    w <- data.frame(w)
    z_probs <- data.frame(z_probs)
    y_probs <- plogis(z_probs[, 1]^2 + z_probs[, 2]^2 - z_probs[, 3] + exp(a + z_probs[, 3] / (1 + rowSums(w)^2)) +
      rnorm(length(a), mean = 0, sd = eps_sd))
    return(y_probs)
  }


  # return DGP functions
  return(list(g_mech = g_mech, z_mech = z_mech, m_mech_cont = m_mech_cont, m_mech_binary = m_mech_binary))
}

make_simulated_W <- function(n_obs = 10000) {
  # baseline covariate -- simple, binary
  W_1 <- rbinom(n_obs, 1, prob = 0.50)
  W_2 <- rbinom(n_obs, 1, prob = 0.65)
  W_3 <- rbinom(n_obs, 1, prob = 0.35)
  W <- cbind(W_1, W_2, W_3)
  return(W)
}

make_simulated_data <- function(n_obs = 10000, binary_outcome = FALSE) { # no. observations
  # get data generating process functions
  dgp <- make_dgp()

  W <- make_simulated_W(n_obs)

  # get probabilities of treatment
  g_mech <- dgp$g_mech(W)

  # create treatment based on baseline W
  A <- as.numeric(rbinom(n_obs, 1, prob = g_mech))

  # get mediator probabilities
  z_mech <- dgp$z_mech(W, A)

  ## 1st mediator (binary)
  Z_1 <- rbinom(n_obs, 1, prob = z_mech[[1]])
  ## 2nd mediator (binary)
  Z_2 <- rbinom(n_obs, 1, prob = z_mech[[2]])
  ## 3rd mediator (binary)
  Z_3 <- rbinom(n_obs, 1, prob = z_mech[[3]])
  ## build matrix of mediators
  Z <- cbind(Z_1, Z_2, Z_3)

  # create outcome as a function of A, Z, W + white noise
  # if y needs to be binary
  if (binary_outcome) {
    y_probs <- dgp$m_mech_binary(W, A, z_probs = cbind(z_mech[[1]], z_mech[[2]], z_mech[[3]]), eps_sd = 0.5)
    Y <- rbinom(n_obs, 1, prob = y_probs)
  } else {
    Y <- dgp$m_mech_cont(W, A, Z, eps_sd = 0.5)
    Y_test <- dgp$m_mech_cont(W, A, z_mech, eps_sd = 0.5)
  }

  # full data structure
  data <- as.data.table(cbind(Y, Z, A, W))
  setnames(data, c(
    "Y", paste("Z", 1:3, sep = "_"), "A",
    paste("W", seq_len(dim(W)[2]), sep = "_")
  ))
  return(data)
}
