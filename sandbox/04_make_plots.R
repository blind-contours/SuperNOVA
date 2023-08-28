library(tidyverse)
library(here)
library(ggpubr)
library(viridis)
library(cividis)
library(hrbrthemes)

## simulation results under discrete exposure
# quant_sim_results_1 <- readRDS(
#   here("sandbox/data/SuperNOVA_mediation_quant_1.rds")
# )

quant_sim_results_2 <- readRDS(
  here("sandbox/data/SuperNOVA_mediation_quant_2.rds")
)

quant_sim_results_3 <- readRDS(
  here("sandbox/data/SuperNOVA_mediation_quant_3.rds")
)
#
quant_sim_results_4 <- readRDS(
  here("sandbox/data/SuperNOVA_mediation_quant_4.rds")
)

quant_sim_results_5 <- readRDS(
  here("sandbox/data/SuperNOVA_mediation_quant_5.rds")
)

quant_sim_results_6 <- readRDS(
  here("sandbox/data/SuperNOVA_mediation_quant_6.rds")
)

quant_sim_results_7 <- readRDS(
  here("sandbox/data/SuperNOVA_mediation_quant_7.rds")
)

quant_sim_results_8 <- readRDS(
  here("sandbox/data/SuperNOVA_mediation_quant_8.rds")
)

quant_sim_results_9 <- readRDS(
  here("sandbox/data/SuperNOVA_mediation_quant_9.rds")
)

quant_sim_results_10 <- readRDS(
  here("sandbox/data/SuperNOVA_mediation_quant_10.rds")
)

## simulation results under continuous exposure using MC integration

cont_sim_mc_results_1 <- readRDS(
  here("sandbox/data/SuperNOVA_mediation_cont_mc_1.rds")
)

cont_sim_mc_results_2 <- readRDS(
  here("sandbox/data/SuperNOVA_mediation_cont_mc_2.rds")
)

cont_sim_mc_results_3 <- readRDS(
  here("sandbox/data/SuperNOVA_mediation_cont_mc_3.rds")
)
#
cont_sim_mc_results_4 <- readRDS(
  here("sandbox/data/SuperNOVA_mediation_cont_mc_4.rds")
)

cont_sim_mc_results_5 <- readRDS(
  here("sandbox/data/SuperNOVA_mediation_cont_mc_5.rds")
)

cont_sim_mc_results_6 <- readRDS(
  here("sandbox/data/SuperNOVA_mediation_cont_mc_6.rds")
)

cont_sim_mc_results_7 <- readRDS(
  here("sandbox/data/SuperNOVA_mediation_cont_mc_7.rds")
)

cont_sim_mc_results_8 <- readRDS(
  here("sandbox/data/SuperNOVA_mediation_cont_mc_8.rds")
)

cont_sim_mc_results_9 <- readRDS(
  here("sandbox/data/SuperNOVA_mediation_cont_mc_9.rds")
)

cont_sim_mc_results_10 <- readRDS(
  here("sandbox/data/SuperNOVA_mediation_cont_mc_10.rds")
)


cont_sim_mc_results <- rbind(cont_sim_mc_results_1, cont_sim_mc_results_2, cont_sim_mc_results_3, cont_sim_mc_results_4, cont_sim_mc_results_5, cont_sim_mc_results_6, cont_sim_mc_results_7, cont_sim_mc_results_8, cont_sim_mc_results_9, cont_sim_mc_results_10)
quant_results <- rbind(quant_sim_results_2, quant_sim_results_3, quant_sim_results_8, quant_sim_results_9)
data_to_plot <- quant_results


sim_statistics <- data_to_plot %>%
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
  )

sim_statistics <- sim_statistics %>%
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

sim_z_statistics <- data_to_plot %>%
  group_by(n_obs) %>%
  mutate(
    nde_std_bias_pseudo = pseudo_reg_nde_bias / sd(pseudo_reg_nde_est),
    nde_std_bias_int = integrated_nde_bias / sd(integrated_nde_est),
    nie_std_bias_pseudo = pseudo_reg_nie_bias / sd(pseudo_reg_nie_est),
    nie_std_bias_int = integrated_nie_bias / sd(integrated_nie_est),
    ate_std_bias = ate_bias / sd(ate_est)
  ) %>%
  select(n_obs, nde_std_bias_pseudo, nde_std_bias_int, nie_std_bias_pseudo, nie_std_bias_int, ate_std_bias)


sim_statistics_long <- sim_statistics %>%
  gather(statistic, value, -c(n_obs))

sim_z_statistics_long <- sim_z_statistics %>%
  gather(statistic, value, -c(n_obs))


sim_statistics_long <- sim_statistics_long %>%
  mutate(method = ifelse(grepl("pseudo", statistic), "pseudo-regression", "integration"))

sim_z_statistics_long <- sim_z_statistics_long %>%
  mutate(method = ifelse(grepl("pseudo", statistic), "pseudo-regression", "integration"))



make_density_plot <- function(sim_statistics_long, z_stat, label) {

  filtered_sim_stats <- sim_statistics_long %>%
    mutate(effect_type = str_replace_all(statistic, "_std_bias_pseudo|_std_bias_int", ""))

  ate_pseudo_rows <- filtered_sim_stats %>%
    filter(effect_type == "ate_std_bias", method == "integration") %>%
    mutate(method = "pseudo-regression")

  filtered_sim_stats <- filtered_sim_stats %>%
    bind_rows(ate_pseudo_rows)

  plot_titles <- c("integration" = "Integration",
                   "pseudo-regression" = "Pseudo-Regression")

  legend_labels <- c("ate_std_bias" = "Total Effect",
                     "nde" = "Natural Direct Effect",
                     "nie" = "Natural Indirect Effect")

  labeler <- labeller(
    effect_type = legend_labels,
    method = plot_titles
  )

  ggplot(filtered_sim_stats, aes(x=value, color=factor(n_obs), fill=factor(n_obs))) +
    geom_density(alpha=0.4) +
    scale_fill_viridis_d(name="Sample Size", option = "plasma") +
    scale_color_viridis_d(name="Sample Size", option = "plasma") +
    theme_minimal() +
    theme(
      legend.position="left",
      panel.spacing = unit(0.1, "lines"),
      strip.text.x = element_text(size = 8)
    ) +
    ggtitle("Quantized Exposure Sampling Distribution for In(Direct) Effect") +
    xlab("Z-Score Standardized Bias (Bias/SD of Estimate)") +
    ylab("Density") +
    facet_grid(effect_type~method, scales = "free", labeller = labeler) +
    coord_cartesian(xlim = c(-1, 1))
}


make_bias_plot <- function(sim_statistics_long, stats, title, legend_labels) {

  filtered_sim_stats <- sim_statistics_long %>%
    filter(statistic %in% stats) %>%
    mutate(effect_type = str_replace_all(statistic, "_pseudo_bias|_int_bias", ""))

  ate_pseudo_rows <- filtered_sim_stats %>%
    filter(effect_type == "abs_ate_bias", method == "integration") %>%
    mutate(method = "pseudo-regression")

  filtered_sim_stats <- filtered_sim_stats %>%
    bind_rows(ate_pseudo_rows)

  filtered_sim_stats <- filtered_sim_stats %>%
    group_by(effect_type, method) %>%
    mutate(
      theoretical_bias = ifelse(n_obs == min(n_obs), value, NA)
    ) %>%
    fill(theoretical_bias) %>%
    ungroup() %>%
    mutate(
      theoretical_bias = theoretical_bias / sqrt(n_obs / min(n_obs))
    )

  plot_titles <- c("integration" = "Integration Method",
                   "pseudo-regression" = "Pseudo-Regression Method")

  legend_labels <- c("abs_ate_bias" = "Total Effect Bias",
                     "abs_nde" = "Natural Direct Effect Bias",
                     "abs_nie" = "Natural Indirect Effect Bias")


  filtered_sim_stats %>%
    ggplot(aes(x = n_obs, y = value, group = effect_type, color = effect_type)) +
    geom_line() +
    geom_line(aes(y = theoretical_bias), linetype = "dashed") +
    facet_wrap(~method, scales = "free", labeller = as_labeller(plot_titles)) +
    scale_color_manual(values = viridis(3, option = "plasma"), labels = legend_labels) +
    labs(
      title = title,
      subtitle = "Dashed lines represent the expected bias decrease if the estimators were sqrt(n)-consistent.",
      x = "Number of Observations",
      y = "Absolute Bias",
      color = "Effect"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.title = element_text(face = "bold"),
      legend.position = "bottom",
      strip.text = element_text(face = "bold")
    )
}


make_sim_coverage_plot <- function(sim_statistics_long, stats, title) {

  filtered_sim_stats <- sim_statistics_long %>%
    filter(statistic %in% stats) %>%
    mutate(effect_type = str_replace_all(statistic, "_pseudo_CI_coverage|_int_CI_coverage", ""))

  ate_pseudo_rows <- filtered_sim_stats %>%
    filter(effect_type == "ate_CI_coverage", method == "integration") %>%
    mutate(method = "pseudo-regression")

  filtered_sim_stats <- filtered_sim_stats %>%
    bind_rows(ate_pseudo_rows)

  plot_titles <- c("integration" = "Integration Method",
                   "pseudo-regression" = "Pseudo-Regression Method")

  legend_labels <- c("ate_CI_coverage" = "Total Effect Coverage",
                     "nde" = "Natural Direct Effect Coverage",
                     "nie" = "Natural Indirect Effect Coverage")

  filtered_sim_stats %>%
    ggplot(aes(x = n_obs, y = value, group = effect_type, color = effect_type)) +
    geom_point(position = position_jitter(width = 50, height = 0)) + # jitter points to avoid overlap
    geom_line() +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") + # add a dashed line at 0.95
    facet_wrap(~method, scales = "free", labeller = as_labeller(plot_titles)) +
    scale_color_manual(values = viridis(3, option = "plasma"), labels = legend_labels) +
    labs(
      title = title,
      subtitle = "Coverage is the proportion of confidence intervals that contain the true value",
      x = "Number of Observations",
      y = "Coverage",
      color = "Effect"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.title = element_text(face = "bold"),
      legend.position = "bottom",
      strip.text = element_text(face = "bold")
    )
}





# Statistics for plotting

bias_stats <- c(abs_ate_bias = "ATE Bias",
                       abs_nde_pseudo_bias = "NDE Bias using Pseudo-Regression",
                       abs_nie_pseudo_bias = "NIE Bias using Pseudo-Regression",
                       abs_nde_int_bias = "NDE Bias using Integration",
                       abs_nie_int_bias = "NIE Bias using Integration")

coverage_stats <- c(nde_pseudo_CI_coverage = "NDE Coverage using Pseudo Regression",
                           nie_pseudo_CI_coverage = "NIE Coverage using Pseudo Regression",
                           nde_int_CI_coverage = "NDE Coverage using Integration",
                           nie_int_CI_coverage = "NIE Coverage using Integration",
                           ate_CI_coverage = "ATE Coverage")

bias_z_stats <- c(nde_std_abs_bias_pseudo = "NDE Standardized Bias using Pseudo Regression",
                  nde_std_abs_bias_int = "NDE Standardized Bias using Integration",
                  nie_std_abs_bias_pseudo = "NIE Standardized Bias using Pseudo Regression",
                  nie_std_abs_bias_int = "NIE Standardized Bias using Integration",
                  ate_std_abs_bias = "Total Effect Standardized Bias")

# Bias
bias_plot <- make_bias_plot(sim_statistics_long = sim_statistics_long,
                                             stats = names(bias_stats),
                                             title = "Continuous Exposure Bias Results for In(Direct) Effects",
                                             legend_labels = bias_stats)

# Coverage
ci_coverage_plot <- make_coverage_plot(sim_statistics_long = sim_statistics_long,
                                        stats = names(coverage_stats),
                                        title = "Continuous Exposure Coverage Results for In(Direct) Effects",
                                        legend_labels = coverage_stats)

# Standardized Bias
z_bias_plot <- make_density_plot(sim_statistics_long = sim_z_statistics_long,
                                       stats = names(bias_z_stats),
                                       title = "Quantized Exposure Standardized Bias Results for In(Direct) Effects",
                                       legend_labels = bias_z_stats)



# Save the plots to files
# Save the plots to files
ggsave(here("sandbox/plots/Pseudo_bias.png"), pseudo_bias_plot, width = 13, height = 7, units = "in", device = "png")
ggsave(here("sandbox/plots/Integration_bias.png"), integration_bias_plot, width = 13, height = 7, units = "in", device = "png")
ggsave(here("sandbox/plots/Pseudo_cov.png"), pseudo_cov_plot, width = 13, height = 7, units = "in", device = "png")
ggsave(here("sandbox/plots/Pseudo_regression_sqrt_n_bias.png"), pseudo_regression_sqrt_n_bias, width = 13, height = 7, units = "in", device = "png")
ggsave(here("sandbox/plots/Int_regression_sqrt_n_bias.png"), int_regression_sqrt_n_bias, width = 13, height = 7, units = "in", device = "png")
