# packages
library(here)
library(devtools)
library(dplyr)
library(SuperNOVA)
library(magrittr)
library(stringr)
source(here("sandbox/02_fit_estimators.R"))

# simulation parameters
n_sim <- 50 # number of simulations
n_obs <- c(250, 500, 1000, 1500, 2000, 2500) # sample sizes at root-n scale

# Generate simulated data -----------------

full_data <- simulate_data(n_obs = 100000, shift_var_index = 1)
p0_data <- full_data$data
effect <- full_data$effect

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

    est_out <- fit_estimators(
      data = as.data.frame(data_sim),
      covars = c("W1", "W2"),
      exposures = c("M1", "M2", "M3"),
      outcome = "Y",
      seed = seed,
      effect_truth = effect,
      deltas = list("M1" = 1, "M2" = 0, "M3" = 0),
      shift_var = "M1",
      cv_folds = 2
    )

    est_out$n_obs <- sample_size

    results[[this_iter]] <- est_out
  }
  # concatenate iterations
  results_out <- bind_rows(results, .id = "sim_iter")
  sim_results_df <- rbind(sim_results_df, results_out)
}


# save results to file
timestamp <- str_replace_all(Sys.time(), " ", "_")
saveRDS(
  object = sim_results_df,
  file = here("sandbox/data", paste0("SuperNOVA_", "sim", ".rds"))
)
