get_sim_truth_NIE_NDE <- function(n_obs = 1e7, binary_outcome = FALSE, EIC = FALSE) {
  dgp <- make_dgp()

  # compute large set of data
  W <- make_simulated_W(n_obs)

  g1 <- dgp$g_mech(W)
  g0 <- 1 - g1
  A <- rbinom(n_obs, 1, prob = g1)
  g <- A * g1 + (1 - A) * g0

  z_probs <- dgp$z_mech(W, A)
  Z_1 <- rbinom(n_obs, 1, prob = z_probs[[1]])
  ## 2nd mediator (binary)
  Z_2 <- rbinom(n_obs, 1, prob = z_probs[[2]])
  ## 3rd mediator (binary)
  Z_3 <- rbinom(n_obs, 1, prob = z_probs[[3]])
  ## build matrix of mediators
  Z <- cbind(Z_1, Z_2, Z_3)

  Z_1_probs <- dgp$z_mech(W, 1)
  Z_0_probs <- dgp$z_mech(W, 0)

  # Z1 counterfactuals
  Z1_1 <- rbinom(n_obs, 1, prob = Z_1_probs[[1]])
  Z1_0 <- rbinom(n_obs, 1, prob = Z_0_probs[[1]])

  # Z2 counterfactuals
  Z2_1 <- rbinom(n_obs, 1, prob = Z_1_probs[[2]])
  Z2_0 <- rbinom(n_obs, 1, prob = Z_0_probs[[2]])

  # Z3 counterfactuals
  Z3_1 <- rbinom(n_obs, 1, prob = Z_1_probs[[3]])
  Z3_0 <- rbinom(n_obs, 1, prob = Z_0_probs[[3]])

  Z1 <- cbind(Z_1 = Z1_1, Z_2 = Z2_1, Z_3 = Z3_1)
  Z0 <- cbind(Z_1 = Z1_0, Z_2 = Z2_0, Z_3 = Z3_0)


  # compute Y counterfactuals
  if (binary_outcome) {
    # Y(1, 1)
    EY_A1_Z1 <- dgp$m_mech_binary(W, 1, Z_1_probs, eps_sd = 0)
    # Y(0, 0)
    EY_A0_Z0 <- dgp$m_mech_binary(W, 0, Z_0_probs, eps_sd = 0)
    # Y(1, 0)
    EY_A1_Z0 <- dgp$m_mech_binary(W, 1, Z_0_probs, eps_sd = 0)
    # Y(0, 1)
    EY_A0_Z1 <- dgp$m_mech_binary(W, 0, Z_1_probs, eps_sd = 0)

    # compute TRUE M under counterfactual regimes
    m_Ais1 <- dgp$m_mech_binary(W, 1, z_probs, eps_sd = 0)
    m_Ais0 <- dgp$m_mech_binary(W, 0, z_probs, eps_sd = 0)

    psi_Z_NDE <- EY_A1_Z0 - EY_A0_Z0
    psi_Z_NIE <- EY_A1_Z1 - EY_A1_Z0
  } else {
    # Y(1, 1)
    EY_A1_Z1 <- dgp$m_mech_cont(W, 1, Z1, eps_sd = 0)
    # Y(0, 0)
    EY_A0_Z0 <- dgp$m_mech_cont(W, 0, Z0, eps_sd = 0)
    # Y(1, 0)
    EY_A1_Z0 <- dgp$m_mech_cont(W, 1, Z0, eps_sd = 0)
    # Y(0, 1)
    EY_A0_Z1 <- dgp$m_mech_cont(W, 0, Z1, eps_sd = 0)

    # compute TRUE M under counterfactual regimes
    m_Ais1 <- dgp$m_mech_cont(W, 1, Z, eps_sd = 0)
    m_Ais0 <- dgp$m_mech_cont(W, 0, Z, eps_sd = 0)

    psi_Z_NDE <- EY_A1_Z0 - EY_A0_Z0
    psi_Z_NIE <- EY_A1_Z1 - EY_A1_Z0
  }


  NDE <- mean(psi_Z_NDE)
  NIE <- mean(psi_Z_NIE)

  if (EIC) {
    df1 <- cbind(Z1, W) %>% as.data.frame()
    df0 <- cbind(Z0, W) %>% as.data.frame()
    Q_Z1 <- df1 %>%
      dplyr::group_by(
        Z_1, Z_1, Z_2, Z_3,
        W_1, W_1, W_2, W_3
      ) %>%
      summarize(
        prob = n() / n_obs
      )
    Q_Z0 <- df0 %>%
      dplyr::group_by(
        Z_1, Z_1, Z_2, Z_3,
        W_1, W_1, W_2, W_3
      ) %>%
      summarize(
        prob = n() / n_obs
      )

    Q_Z_Ais1 <- df1 %>%
      dplyr::left_join(
        Q_Z1
      ) %>%
      pull(prob)

    Q_Z_Ais0 <- df0 %>%
      dplyr::left_join(
        Q_Z0
      ) %>%
      pull(prob)

    Q_Z_ratio <- Q_Z_Ais0 / Q_Z_Ais1

    if (binary_outcome) {
      Y <- dgp$m_mech_binary(W, A, z_probs)
      Qbar_Y <- dgp$m_mech_binary(W, A, Z, eps_sd = 0)
    } else {
      Y <- dgp$m_mech_cont(W, A, Z)
      Qbar_Y <- dgp$m_mech_cont(W, A, Z, eps_sd = 0)
    }

    g0 <- 1 - g1

    ###################################################
    ## NDE

    HY <- (A / g1) * (Q_Z_Ais0 / Q_Z_Ais1) - ((1 - A) / g0)

    HZ <- 1 / g0

    # compute individual scores for DY, DA, DW
    D_Y_NDE <- HY * (Y - Qbar_Y)
    D_Z_NDE <- (1 - A) * HZ * (m_Ais1 - m_Ais0 - psi_Z_NDE)
    D_W_NDE <- psi_Z_NDE - NDE

    EIC_NDE <- D_Y_NDE + D_Z_NDE + D_W_NDE

    ###################################################
    ## NIE

    EIC_NIE <- ((A / g1) * (Y - EY_A1_Z1 - Q_Z_ratio) * (Y - m_Ais1)) -
      ((1 - A) / g0) * (m_Ais1 - EY_A1_Z0) +
      psi_Z_NIE - NIE

    # output: true values of nuisance parameters
    return(list(
      EY_A1_Z1 = EY_A1_Z1,
      EY_A1_Z0 = EY_A1_Z0,
      EY_A0_Z1 = EY_A0_Z1,
      EY_A0_Z0 = EY_A0_Z0,
      NDE = NDE,
      NIE = NIE,
      EIC_NDE = EIC_NDE,
      EIC_NIE = EIC_NIE
    ))
  } else {
    # output: true values of nuisance parameters
    return(list(
      EY_A1_Z1 = EY_A1_Z1,
      EY_A1_Z0 = EY_A1_Z0,
      EY_A0_Z1 = EY_A0_Z1,
      EY_A0_Z0 = EY_A0_Z0,
      NDE = NDE,
      NIE = NIE
    ))
  }
}
