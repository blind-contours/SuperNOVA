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
n_sim <- 50 # number of simulations
n_obs <- c(250, 500, 1000, 1500, 2500, 3000)
p0_obs <- 100000

n_core <- 20
n_fold <- 10

# Generate simulated data -----------------

full_data <- simulate_complicated_mediation_data(n_obs = 100000,
                                     delta = 1)
p0_data <- full_data

covars <- c("w_1", "w_2", "w_4", "w_5")
exposures <- c("a_1", "a_2", "a_3", "a_4", "a_5")
mediators <- c("z_1", "z_2", "z_3", "z_4", "z_5")
outcome <- "y"

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

    sim_results <- as.data.frame(sim_results)

    sim_results$n_obs <- sample_size

    results[[this_iter]] <- sim_results
  }
  # concatenate iterations
  results_out <- bind_rows(results, .id = "sim_iter")
  sim_results_df <- rbind(sim_results_df, results_out)
}

avg_freq_df <- sim_results_df %>%
  group_by(Var1, n_obs) %>%
  summarize(avg_freq = mean(Freq), .groups = 'drop')


# save results to file
saveRDS(
  object = avg_freq_df,
  file = here("sandbox/data", paste0("SuperNOVA_", "mediation_quant_1", ".rds"))
)
