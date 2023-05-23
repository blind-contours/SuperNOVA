library(tidyverse)
library(here)
library(ggpubr)
library(viridis)
library(cividis)
library(hrbrthemes)

sim_results_1 <- readRDS(
  here("sandbox/data/SuperNOVA_mediation_quant_1.rds")
)

sim_results_2 <- readRDS(
  here("sandbox/data/SuperNOVA_mediation_quant_2.rds")
)

sim_results_3 <- readRDS(
  here("sandbox/data/SuperNOVA_mediation_quant_3.rds")
)
#
sim_results_4 <- readRDS(
  here("sandbox/data/SuperNOVA_mediation_quant_4.rds")
)

sim_results_5 <- readRDS(
  here("sandbox/data/SuperNOVA_mediation_quant_5.rds")
)


sim_results <- rbind(sim_results_1, sim_results_2, sim_results_3, sim_results_4)


sim_statistics <- sim_results %>%
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


sim_statistics_long <- sim_statistics %>%
  gather(statistic, value, -c(n_obs))

n_obs <- c(250, 500, 1000, 1500, 2500, 3000)
inverse_sqrt <- 1 /sqrt(n_obs)
comp_condition <- "1/sqrt(n)"
ref_group <- bind_cols(n_obs, comp_condition, inverse_sqrt)
colnames(ref_group) <- c("n_obs", "statistic", "value")

sim_statistics_long <- rbind(sim_statistics_long, ref_group)


make_density_plot <- function(sim_statistics,
                              z_stat,
                              label) {


  sim_statistics$n_obs <- as.factor(sim_statistics$n_obs)

  ggplot(sim_statistics, aes(x=eval(parse(text = "joint_z")), color=n_obs, fill=n_obs)) +
    geom_density(alpha=0.4) +
    scale_fill_viridis(discrete=TRUE , option = "plasma") +
    scale_color_viridis(discrete=TRUE, option = "plasma") +
    theme_ipsum() +
    theme(
      legend.position="left",
      panel.spacing = unit(0.1, "lines"),
      strip.text.x = element_text(size = 8)
    ) +
    ggtitle(label) +
    xlab("Z-Score Standardized Bias (Bias/SD)") +
    ylab("Probability") +
    stat_function(fun = dnorm,
                  args = list(mean = 0,
                              sd = 1),
                  col = "#1b98e0",
                  size = 1
    )

  }

sim_statistics_long <- rbind(sim_statistics_long, ref_group)

make_sim_statistics_plot <- function(sim_statistics_long, stats, title, legend_labels) {
  filtered_sim_stats <- sim_statistics_long %>%
    filter(statistic %in% stats)

  # Create the plot with the "inferno" color scale
  ggplot(filtered_sim_stats,
         aes(x = n_obs, y = value, group = statistic, color = statistic)) +
    geom_line() +
    # stat_smooth(aes(y=value, x=n_obs), method = lm, formula = y ~ poly(x, 2), se = FALSE) +

    scale_colour_viridis_d(option = "F", labels = legend_labels) +
    ggtitle(title) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = "black")) +
    ylab("Bias") +
    xlab("Number Observations")
}

# Statistics for plotting

pseudo_bias_stats <- c(abs_ate_bias = "ATE Bias",
                       abs_nde_pseudo_bias = "NDE Bias using Pseudo-Regression",
                       abs_nie_pseudo_bias = "NIE Bias using Pseudo-Regression")

int_bias_stats <- c(abs_ate_bias = "ATE Bias",
                    abs_nde_int_bias = "NDE Bias using Integration",
                    abs_nie_int_bias = "NIE Bias using Integration")

pseudo_sqrt_n_bias_stats <- c(nde_pseudo_sqrt_n_abs_bias = "Sqrt(n) * Bias using Pseudo-Regression for NDE",
                              nie_pseudo_sqrt_n_abs_bias = "Sqrt(n) * Bias using Pseudo-Regression for NIE",
                              ate_sqrt_n_abs_bias = "Sqrt(n) * Bias for ATE")

int_sqrt_n_bias_stats <- c(nde_int_sqrt_n_abs_bias = "Sqrt(n) * Bias using Integration for NDE",
                          nie_int_sqrt_n_abs_bias = "Sqrt(n) * Bias using Integration for NIE",
                          ate_sqrt_n_abs_bias = "Sqrt(n) * Bias for ATE")


pseudo_mse_stats <- c(nde_pseudo_n_MSE = "NDE n * MSE using Pseudo-Regression",
                      nie_pseudo_n_MSE = "NIE n * MSE using Pseudo-Regression",
                      ate_n_MSE = "ATE n * MSE")

int_mse_stats <- c(nde_int_n_MSE = "NDE MSE using Integration",
                   nie_int_n_MSE = "NIE MSE using Integration",
                   ate_n_mSE = "ATE n * MSE")

pseudo_coverage_stats <- c(nde_pseudo_CI_coverage = "NDE Coverage using Pseudo Regression",
                    nie_pseudo_CI_coverage = "NIE Coverage using Pseudo Regression",
                    ate_CI_coverage = "ATE Coverage")

int_coverage_stats <- c(nde_int_CI_coverage = "NDE Coverage using Integration",
                        nie_int_CI_coverage = "NIE Coverage using Integration",
                        ate_CI_coverage = "ATE Coverage")


# Bias

# Pseudo-Regression
pseudo_bias_plot <- make_sim_statistics_plot(sim_statistics_long, names(pseudo_bias_stats),
                                             "Bias Pseudo Regressions Results for In(Direct) Effects", pseudo_bias_stats)

# Integration
integration_bias_plot <- make_sim_statistics_plot(sim_statistics_long, names(int_bias_stats),
                                             "Bias Integration Results for In(Direct) Effects", int_bias_stats)

# Sqrt(N) * bias

# Pseudo-Regression
pseudo_regression_sqrt_n_bias <- make_sim_statistics_plot(sim_statistics_long, stats = names(pseudo_sqrt_n_bias_stats),
                                                  title = "Sqrt(n) * Bias Pseudo Regression Results for In(Direct) Effects", legend_labels = pseudo_sqrt_n_bias_stats)

# Integration
int_regression_sqrt_n_bias <- make_sim_statistics_plot(sim_statistics_long, stats = names(int_sqrt_n_bias_stats),
                                                       title = "Sqrt(n) * Bias Integration Results for In(Direct) Effects", legend_labels = int_sqrt_n_bias_stats)


# Coverage

# Pseudo-Regression
pseudo_cov_plot <- make_sim_statistics_plot(sim_statistics_long, names(pseudo_coverage_stats),
                                     "MSE for In(Direct) Effects", pseudo_coverage_stats)

# Integration
integration_cov_plot <- make_sim_statistics_plot(sim_statistics_long, names(int_coverage_stats),
                                            "Coverage for In(Direct) Effects using Integration", int_coverage_stats)

# Save the plots to files
# Save the plots to files
ggsave(here("sandbox/plots/Pseudo_bias.png"), pseudo_bias_plot, width = 13, height = 7, units = "in", device = "png")
ggsave(here("sandbox/plots/Integration_bias.png"), integration_bias_plot, width = 13, height = 7, units = "in", device = "png")
ggsave(here("sandbox/plots/Pseudo_cov.png"), pseudo_cov_plot, width = 13, height = 7, units = "in", device = "png")
ggsave(here("sandbox/plots/Pseudo_regression_sqrt_n_bias.png"), pseudo_regression_sqrt_n_bias, width = 13, height = 7, units = "in", device = "png")
ggsave(here("sandbox/plots/Int_regression_sqrt_n_bias.png"), int_regression_sqrt_n_bias, width = 13, height = 7, units = "in", device = "png")

