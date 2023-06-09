# Load the necessary libraries
library(tidyverse)
library(here)
library(ggpubr)
library(viridis)
library(cividis)
library(hrbrthemes)

# Load data from given paths
load_data <- function(paths) {
  purrr::map_df(paths, readRDS)
}

compute_sim_stats <- function(data) {
  data %>%
    group_by(n_obs) %>%
    summarize(
      mean_nde_bias_pseudo = mean(pseudo_reg_nde_bias),
      nde_est_sd_pseudo = sd(pseudo_reg_nde_est),
      nde_est_MSE_pseudo = mean_nde_bias_pseudo^2 + nde_est_sd_pseudo^2,
      nde_pseudo_CI_coverage = mean(pseudo_reg_nde_cov),
      #
      mean_nde_bias_int = mean(integrated_nde_bias),
      nde_est_sd_int = sd(integrated_nde_est),
      nde_est_MSE_int = mean_nde_bias_int^2 + nde_est_sd_int^2,
      nde_int_CI_coverage = mean(integrated_nde_cov),
      #
      mean_nie_bias_pseudo = mean(pseudo_reg_nie_bias),
      nie_est_sd_pseudo = sd(pseudo_reg_nie_est),
      nie_est_MSE_pseudo = mean_nie_bias_pseudo^2 + nie_est_sd_pseudo^2,
      nie_pseudo_CI_coverage = mean(pseudo_reg_nie_cov),
      #
      mean_nie_bias_int = mean(integrated_nie_bias),
      nie_est_sd_int = sd(integrated_nie_est),
      nie_est_MSE_int = mean_nie_bias_int^2 + nie_est_sd_int^2,
      nie_int_CI_coverage = mean(integrated_nie_cov),
      #
      mean_ate_bias = mean(ate_bias),
      ate_est_sd = sd(ate_est),
      ate_est_MSE = mean_ate_bias^2 + ate_est_sd^2,
      ate_CI_coverage = mean(ate_cov)
      #
    ) %>%
    mutate(
      abs_nde_pseudo_bias = abs(mean_nde_bias_pseudo),
      abs_nde_int_bias = abs(mean_nde_bias_int),
      abs_nie_pseudo_bias = abs(mean_nie_bias_pseudo),
      abs_nie_int_bias = abs(mean_nie_bias_int),
      abs_ate_bias = abs(mean_ate_bias),
      #
      nde_pseudo_sqrt_n_abs_bias = sqrt(n_obs) * abs_nde_pseudo_bias,
      nde_int_sqrt_n_abs_bias = sqrt(n_obs) * abs_nde_int_bias,
      nie_pseudo_sqrt_n_abs_bias = sqrt(n_obs) * abs_nie_pseudo_bias,
      nie_int_sqrt_n_abs_bias = sqrt(n_obs) * abs_nie_int_bias,
      ate_sqrt_n_abs_bias = sqrt(n_obs) * abs_ate_bias,
      #
      nde_pseudo_n_MSE = n_obs * nde_est_MSE_pseudo,
      nde_int_n_MSE = n_obs * nde_est_MSE_int,
      nie_pseudo_n_MSE = n_obs * nie_est_MSE_pseudo,
      nie_int_n_MSE = n_obs * nie_est_MSE_int,
      ate_n_MSE = n_obs * ate_est_MSE
    )
}

# Define paths
quant_paths <- paste0(here("sandbox/data/SuperNOVA_mediation_quant_"), 2:5, ".rds")
cont_paths <- paste0(here("sandbox/data/SuperNOVA_mediation_cont_mc_"), 1:10, ".rds")

# Load the data
quant_sim_results <- load_data(quant_paths)
cont_sim_results <- load_data(cont_paths)

# Compute statistics
quant_sim_stats <- compute_sim_stats(quant_sim_results)
cont_sim_stats <- compute_sim_stats(cont_sim_results)

# Convert to long form for plotting
sim_statistics_long <- bind_rows(list(quant_sim_stats = quant_sim_stats, cont_sim_stats = cont_sim_stats), .id = "data_set") %>%
  pivot_longer(-c(n_obs, data_set), names_to = "statistic", values_to = "value") %>%
  mutate(type = ifelse(str_detect(statistic, "pseudo"), "pseudo_regression", "integration"))

# Plot the statistics
plot_sim_stats <- function(data, stats, data_set_type, title) {
  data %>%
    filter(data_set_type == data_set, statistic %in% stats) %>%
    ggplot(aes(x = n_obs, y = value, group = statistic, color = statistic)) +
    geom_line() +
    facet_wrap(~type, ncol = 2, scales = "free_y") +
    theme_few() +
    ggtitle(title) +
    ylab("Value") +
    xlab("Number Observations") +
    scale_colour_viridis_d(option = "plasma", labels = stats)
}

# Define the statistics for plotting
pseudo_bias_stats <- c("abs_ate_bias",
                       "abs_nde_pseudo_bias",
                       "abs_nie_pseudo_bias")

int_bias_stats <- c("abs_ate_bias",
                    "abs_nde_int_bias",
                    "abs_nie_int_bias")

bias_stats <- c(pseudo_bias_stats,int_bias_stats)
sqrt_n_bias_stats <- c("nde_pseudo_sqrt_n_abs_bias", "nie_pseudo_sqrt_n_abs_bias", "ate_sqrt_n_abs_bias")
mse_stats <- c("nde_pseudo_n_MSE", "nie_pseudo_n_MSE", "ate_n_MSE")
sd_stats <- c("nde_est_sd_pseudo", "nie_est_sd_pseudo", "ate_est_sd")

# Plot the defined statistics
plot_sim_stats(data = sim_statistics_long,
               stats = bias_stats,
               data_set_type = "quant_sim_stats",
               title = "Bias Plot")

plot_sim_stats(data = sim_statistics_long,
               stats = sqrt_n_bias_stats,
               data_set = "quant_sim_stats",
               title = "Sqrt N Bias Plot")

plot_sim_stats(data = sim_statistics_long,
               stats = mse_stats,
               data_set = "quant_sim_stats",
               title = "MSE Plot")

plot_sim_stats(data = sim_statistics_long,
               stats = sd_stats,
               data_set = "quant_sim_stats",
               title = "SD Plot")
