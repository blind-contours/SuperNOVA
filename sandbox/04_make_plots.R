library(tidyverse)
library(here)
library(ggpubr)
library(viridis)
library(cividis)
library(hrbrthemes)

sim_results_1 <- readRDS(
  here("sandbox/data/SuperNOVA_mediation_cont_mc.rds")
)

# sim_results_2 <- readRDS(
#   here("sandbox/data/SuperNOVA_2_sim.rds")
# )
#
# sim_results_3 <- readRDS(
#   here("sandbox/data/SuperNOVA_3_sim.rds")
# )
#
# sim_results_4 <- readRDS(
#   here("sandbox/data/SuperNOVA_5_sim.rds")
# )
#
# sim_results_5 <- readRDS(
#   here("sandbox/data/SuperNOVA_6_sim.rds")
# )
#
# sim_results_6 <- readRDS(
#   here("sandbox/data/SuperNOVA_7_sim.rds")
# )


sim_results <- rbind(sim_results_1)


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
    nie_int_CI_coverage = mean(integrated_nie_cov)
    #
  ) %>%
  mutate(
    abs_nde_pseudo_bias = abs(mean_nde_bias_pseudo),
    abs_nde_int_bias = abs(mean_nde_bias_int),
    abs_nie_pseudo_bias = abs(mean_nie_bias_pseudo),
    abs_nie_int_bias = abs(mean_nie_bias_int),
    #
    nde_pseudo_sqrt_n_abs_bias = sqrt(n_obs) * abs_nde_pseudo_bias,
    nde_int_sqrt_n_abs_bias = sqrt(n_obs) * abs_nde_int_bias,
    nie_pseudo_sqrt_n_abs_bias = sqrt(n_obs) * abs_nie_pseudo_bias,
    nie_int_sqrt_n_abs_bias = sqrt(n_obs) * abs_nie_int_bias,
    #
    nde_pseudo_n_MSE = n_obs * nde_est_MSE_pseudo,
    nde_int_n_MSE = n_obs * nde_est_MSE_int,
    nie_pseudo_n_MSE = n_obs * nie_est_MSE_pseudo,
    nie_int_n_MSE = n_obs * nie_est_MSE_int
  ) %>%
  ungroup()


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

make_sim_statistics_plot <- function(sim_statistics_long, stats, title) {
  filtered_sim_stats <- sim_statistics_long %>%
    filter(statistic %in% stats)

  # Create the plot with the "inferno" color scale
  ggplot(filtered_sim_stats,
         aes(x = n_obs, y = value, group = statistic, color = statistic)) +
    geom_line(position = position_jitter(w = 0.0, h = 0.0)) +
    scale_colour_viridis_d(option = "inferno") +
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

# NDE bias plot
pseudo_stats <- c("mean_nde_bias_pseudo", "mean_nie_bias_pseudo", "1/sqrt(n)")
pseudo_bias_plot <- make_sim_statistics_plot(sim_statistics_long, pseudo_stats, "Pseudo Regressions")

# NIE bias plot
integration_stats <- c("mean_nde_bias_int", "mean_nie_bias_int", "1/sqrt(n)")
nie_bias_plot <- make_sim_statistics_plot(sim_statistics_long, integration_stats, "NIE Bias")


# Integration CI plot
int_cov_stats <- c("nde_int_CI_coverage", "nie_int_CI_coverage", "1/sqrt(n)")
int_cov_plot <- make_sim_statistics_plot(sim_statistics_long, int_cov_stats, "Coverage with Integration")


# Save the plots to files
ggsave(
  here("sandbox/plots/NDE_bias.png"),
  nde_bias_plot,
  width = 13,
  height = 7,
  device = png()
)

ggsave(
  here("sandbox/plots/NIE_bias.png"),
  nie_bias_plot,
  width = 13,
  height = 7,
  device = png()
)
