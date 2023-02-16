# be sure to set this env variable by `export R_LIBDIR=/path/to/your/R/libs`
r_libdir <- Sys.getenv("R_LIBDIR")

# set user-specific package library
if (grepl("savio2", Sys.info()["nodename"])) {
  .libPaths(r_libdir)
  Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
}

# packages
library(here)
library(foreach)
library(future)
library(doFuture)
library(doRNG)
library(data.table)
library(tidyverse)
library(hal9001)
library(origami)
library(sl3)
library(tmle3)
library(SuperNOVA)
devtools::load_all(here())
source(here("sandbox", "02_fit_estimators.R"))


# simulation parameters
set.seed(7259)
n_sim <- 5 # number of simulations
n_obs <- c(500, 1000, 1500) # sample sizes at root-n scale

# Generate simulated data -----------------

data <- simulate_mediation_data(n_obs = 100000)
nde_truth <- data$nde
nie_truth <- data$nie
ate_a_1 <- data$ate_a_1
ate_a_2 <- data$ate_a_2
a2_effect_in_w2_1 <- data$effect_mod_a2w2_lvl_1
a2_effect_in_w2_0 <- data$effect_mod_a2w2_lvl_0

p0_data <- data$data

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
      covars = c("w_1", "w_2", "w_3"),
      exposures = c("a_1", "a_2"),
      mediators = c("z_1"),
      outcome = "y",
      seed = seed,
      nde_truth = nde_truth,
      nie_truth = nie_truth,
      ate_a_1 = ate_a_1,
      ate_a_2 = ate_a_2,
      a2_effect_in_w2_1 = a2_effect_in_w2_1,
      a2_effect_in_w2_0 = a2_effect_in_w2_0,
      deltas = list("a_1" = 1, "a_2" = 1),
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
  file = here("sandbox/data", paste0("CVtreeMLE_", "sim", ".rds"))
)
