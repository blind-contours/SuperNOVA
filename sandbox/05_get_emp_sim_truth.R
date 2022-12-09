# load scripts, parallelization, PRNG
source(here("sandbox/01_setup_data.R"))
source(here("sandbox/02_fit_estimators.R"))
source(here("/tests/testthat/helper-dgp.R"))

timestamp <- '2020-04-20_16:57:52'
sample_size <- 2000
data_sim <- make_simulated_data(n_obs = sample_size)

################################################################################
get_sim_truth_NIE_NDE <- function(data = data_sim, n_obs) { # number of observations
  # compute large data set for true values

  w_names <- str_subset(colnames(data), "W")
  z_names <- str_subset(colnames(data), "Z")
  W <- subset(data, select =w_names)
  Z <- subset(data, select = z_names)

  dgp <- make_dgp()
  Z_1_probs <- dgp$z_mech(W, 1)
  Z_0_probs <- dgp$z_mech(W, 0)

  WZ_0_probs <- dgp$z_mech(W, 0)

  # Z1 counterfactuals
  Z1_1 <- rbinom(n_obs, 1, prob = Z_1_probs[[1]])
  Z1_0 <- rbinom(n_obs, 1, prob = Z_0_probs[[1]])

  # Z2 counterfactuals
  Z2_1 <- rbinom(n_obs, 1, prob = Z_1_probs[[2]])
  Z2_0 <- rbinom(n_obs, 1, prob = Z_0_probs[[2]])

  # Z3 counterfactuals
  Z3_1 <- rbinom(n_obs, 1, prob = Z_1_probs[[3]])
  Z3_0 <- rbinom(n_obs, 1, prob = Z_0_probs[[3]])

  Z1 <- cbind(Z1_1, Z2_1, Z3_1)
  Z0 <- cbind(Z1_0, Z2_0, Z3_0)

  # compute Y counterfactuals
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

  # output: true values of nuisance parameters
  return(list(
    EY_A1_Z1 = EY_A1_Z1,
    EY_A1_Z0 = EY_A1_Z0,
    EY_A0_Z1 = EY_A0_Z1,
    EY_A0_Z0 = EY_A0_Z0
  ))
}

# simulate data and extract components for computing true parameter value
sim_truth <- get_sim_truth_NIE_NDE(data = sim_data, n_obs = 10000)

EY_A1_Z1 <- sim_truth$EY_A1_Z1
EY_A1_Z0 <- sim_truth$EY_A1_Z0
EY_A0_Z1 <- sim_truth$EY_A0_Z1
EY_A0_Z0 <- sim_truth$EY_A0_Z0

# compute true NIE via empirical substitution estimator
psi_NDE_true <- mean(EY_A1_Z0 - EY_A0_Z0)
psi_NDE_true
psi_NIE_true <- mean(EY_A1_Z1 - EY_A1_Z0)
psi_NIE_true

# load simulated data

sim_results <- readRDS(file = here("data", paste0("tmle3mediate_", timestamp, ".rds")))




