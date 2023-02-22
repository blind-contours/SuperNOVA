# packages
library(here)
library(devtools)
library(dplyr)
library(magrittr)
library(stringr)
source(here("sandbox/02_fit_estimators.R"))
library(SuperNOVA)

# simulation parameters
n_sim <- 10 # number of simulations
n_obs <- c(250, 500, 1000, 1500, 2000, 2500, 3000, 5000) # sample sizes at root-n scale
p0_obs <- 100000

# Generate simulated data -----------------

full_data <- simulate_data(n_obs = 100000,
                           delta = 1)
p0_data <- full_data$data

m1_effect <- full_data$m1_effect
m2_effect <- full_data$m2_effect
m3_effect <- full_data$m3_effect
m4_effect <- full_data$m4_effect
m1m4_effect <- full_data$m14_effect
m1m4_intxn <- full_data$m14_intxn
effect_modification <- full_data$effect_mod

covars <- c("W1", "W2", "W3")
exposures <- c("M1", "M2", "M3", "M4")
outcome <- "Y"


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
    y <- data_sim[, outcome]

    est_out <- fit_estimators(w = w,
                              a = a,
                              y = y,
                              seed = seed,
                              true_effects = c(m1_effect, m2_effect, m3_effect, m4_effect),
                              m14_effect_truth = m1m4_effect,
                              m14_intxn_truth = m1m4_intxn,
                              true_em_effects = effect_modification,
                              deltas = list("M1" = 1, "M2" = 1, "M3" = 1, "M4" = 1),
                              cv_folds = 2,
                              var_sets = c("M1", "M2", "M3", "M4", "M1M4", "M3W3")
    )

    est_out$n_obs <- sample_size

    results[[this_iter]] <- est_out
  }
  # concatenate iterations
  results_out <- bind_rows(results, .id = "sim_iter")
  sim_results_df <- rbind(sim_results_df, results_out)
}


# save results to file
saveRDS(
  object = sim_results_df,
  file = here("sandbox/data", paste0("SuperNOVA_", "4_sim", ".rds"))
)
