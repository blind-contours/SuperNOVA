# packages
library(here)
library(devtools)
library(dplyr)
library(magrittr)
library(stringr)
source(here("sandbox/02_fit_estimators.R"))
source(here("sandbox/simulate_discrete_data.R"))
load_all()

# simulation parameters
n_sim <- 10 # number of simulations
n_obs <- c(250, 500, 1000, 1500, 2500, 3000)
p0_obs <- 100000
n_mc_sample <- 2000

n_core <- 20
n_fold <- 10

# Generate simulated data -----------------

full_data <- simulate_mediation_data(n_obs = 100000,
                                     delta = 1)
p0_data <- full_data$data

nde_a1_quant <- full_data$nde_a1_quant
nie_a1_quant <- full_data$nie_a1_quant
ate_a1_quant <- full_data$ate_a1_quant

covars <- c("w_1", "w_2")
exposures <- c("a_1_quant")
mediators <- c("z_1_quant")
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

    est_out_discrete_m <- fit_estimators_mediation(
      w = w,
      a = a,
      z = z,
      y = y,
      seed = seed,
      nde_effects = c(nde_a1_quant),
      nie_effects = c(nie_a1_quant),
      ate_effects = c(ate_a1_quant),
      deltas = list("a" = 1),
      cv_folds = n_fold,
      num_cores = n_core,
      var_sets = "a-z",
      exposure_quantized = TRUE,
      mediator_quantized = FALSE,
      n_mc_sample = n_mc_sample,
      density_type = "sl",
      integration_method = "MC"
    )

    est_out_discrete_m$n_obs <- sample_size

    results[[this_iter]] <- est_out_discrete_m
  }
  # concatenate iterations
  results_out <- bind_rows(results, .id = "sim_iter")
  sim_results_df <- rbind(sim_results_df, results_out)
}


# save results to file
saveRDS(
  object = sim_results_df,
  file = here("sandbox/data", paste0("SuperNOVA_", "mediation_quant_4", ".rds"))
)
