# packages
library(here)
library(devtools)
library(dplyr)
library(magrittr)
library(stringr)
source(here("sandbox/02_fit_estimators.R"))
source(here("R/simulate_complicated_mediation_data.R"))
load_all()

# simulation parameters
n_sim <- 10 # number of simulations
n_obs <- c(250, 500, 1000, 1500, 2500, 3000)
p0_obs <- 100000

n_core <- 20
n_fold <- 10

# Generate simulated data -----------------

full_data <- simulate_complicated_mediation_data(n_obs = 100000,
                                     delta = 1,
                                     n_bins = 10)
p0_data <- full_data$data

nde_a1 <- full_data$nde_a1
nie_a1 <- full_data$nie_a1
ate_a1 <- full_data$ate_a1

nde_a2 <- full_data$nde_a2
nie_a2 <- full_data$nie_a2
ate_a2 <- full_data$ate_a2

covars <- c("w_1", "w_2", "w_4", "w_5")
exposures <- c("a_1_quant", "a_2_quant", "a_3_quant", "a_4_quant")
mediators <- c("z_1_quant", "z_2_quant", "z_3_quant")
outcome <- "y_quant"

# perform simulation across sample sizes
sim_results_df <- data.frame()

for (sample_size in n_obs) {
  # get results in parallel
  results <- list()
  print(sample_size)

  for (this_iter in seq_len(n_sim)) {
    print(this_iter)

    seed <- sample(1:10000, 1)
    set.seed(seed)

    data_sim <- p0_data %>%
      dplyr::slice_sample(n = sample_size)

    w <- data_sim[, covars]
    a <- data_sim[, exposures]
    z <- data_sim[, mediators]
    y <- data_sim[, outcome]

    sim_results <- SuperNOVA(
      w = w,
      a = a,
      z = z,
      y = y,
      deltas = deltas,
      estimator = "tmle",
      fluctuation = "standard",
      n_folds = n_fold,
      outcome_type = "continuous",
      quantile_thresh = 0,
      verbose = TRUE,
      parallel = TRUE,
      seed = seed,
      num_cores = n_core,
      n_mc_sample = n_mc_sample,
      discover_only = TRUE
    )

    est_out_discrete_e$n_obs <- sample_size

    results[[this_iter]] <- est_out_discrete_e
  }
  # concatenate iterations
  results_out <- bind_rows(results, .id = "sim_iter")
  sim_results_df <- rbind(sim_results_df, results_out)
}


# save results to file
saveRDS(
  object = sim_results_df,
  file = here("sandbox/data", paste0("SuperNOVA_", "mediation_quant_1", ".rds"))
)
