library(tidyverse)
library(here)
library(ggpubr)
library(viridis)
library(cividis)
library(hrbrthemes)

sim_results_1 <- readRDS(
  here("sandbox/data/SuperNOVA_1_sim.rds")
)

sim_results_2 <- readRDS(
  here("sandbox/data/SuperNOVA_2_sim.rds")
)

sim_results_3 <- readRDS(
  here("sandbox/data/SuperNOVA_3_sim.rds")
)

sim_results_4 <- readRDS(
  here("sandbox/data/SuperNOVA_5_sim.rds")
)

sim_results_5 <- readRDS(
  here("sandbox/data/SuperNOVA_6_sim.rds")
)

sim_results_6 <- readRDS(
  here("sandbox/data/SuperNOVA_7_sim.rds")
)


sim_results <- rbind(sim_results_1, sim_results_2, sim_results_3, sim_results_4, sim_results_5, sim_results_6)


sim_statistics <- sim_results %>%
  group_by(n_obs) %>%
  summarize(
    mean_indiv_bias = mean(indiv_bias),
    indiv_est_sd = sd(indiv_est),
    indiv_est_MSE = mean_indiv_bias^2 + indiv_est_sd^2,
    indiv_CI_coverage = mean(indiv_cov),
    indiv_z = indiv_bias / indiv_est_sd,
    mean_effect_mod_bias = mean(em_bias),
    effect_mod_sd = sd(em_bias),
    effect_est_MSE = mean_effect_mod_bias^2 + effect_mod_sd^2,
    effect_mod_z = em_bias / effect_mod_sd,
    em_cov = mean(em_cov),
    mean_joint_bias = mean(joint_bias),
    joint_est_sd =  sd(joint_est),
    joint_est_MSE = mean_joint_bias^2 + joint_est_sd^2,
    joint_cov = mean(joint_cov),
    joint_z = joint_bias / joint_est_sd,
    mean_intxn_bias = mean(intxn_bias),
    intxn_est_sd = sd(intxn_est),
    intxn_est_MSE = mean_intxn_bias^2 + intxn_est_sd^2,
    intxn_cov = mean(intxn_cov),
    intxn_z = intxn_bias/intxn_est_sd
  ) %>%
  mutate(
    abs_indiv_bias = abs(mean_indiv_bias),
    abs_effect_mod_bias = abs(mean_effect_mod_bias),
    abs_joint_abs_bias =  abs(mean_joint_bias),
    abs_intxn_bias = abs(mean_intxn_bias),
    indiv_sqrt_n_abs_bias = sqrt(n_obs) * abs_indiv_bias,
    effect_mod_sqrt_n_abs_bias = sqrt(n_obs) * abs_effect_mod_bias,
    joint_sqrt_n_abs_bias = sqrt(n_obs) * abs_joint_abs_bias,
    intxn_sqrt_n_abs_bias = sqrt(n_obs) * abs_intxn_bias,
    indiv_n_MSE = n_obs * indiv_est_MSE,
    em_n_MSE = n_obs * effect_est_MSE,
    joint_n_MSE = n_obs * joint_est_MSE,
    intxn_n_MSE = n_obs * intxn_est_MSE
    ) %>%
  ungroup()

sim_statistics_long <- sim_statistics %>%
  gather(statistic, value, -c(n_obs))

n_obs <- c(250, 500, 1000, 1500, 2000, 2500, 3000, 5000)
inverse_sqrt <- 1 /sqrt(n_obs)
comp_condition <- "1/sqrt(n)"
ref_group <- bind_cols(n_obs, comp_condition, inverse_sqrt)
colnames(ref_group) <- c("n_obs", "statistic", "value")

sim_statistics_long <- rbind(sim_statistics_long, ref_group)

make_sim_statistics_plot <- function(sim_statistics_long, stats) {

  filtered_sim_stats <- sim_statistics_long %>%
    filter(
      statistic %in% stats

    )

  # Create the plot with the "cividis" color scale
  ggplot(filtered_sim_stats,
         aes(x=n_obs, y=value, group=statistic, color=statistic)) +
    geom_line(position = position_jitter(w=0.0, h=0.0)) +
    scale_colour_viridis_d(option = "plasma") +
    ggtitle("Simulation Results for Individual Shift") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = "black")) +
    ylab("Bias") +
    xlab("Number Observations")
}


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

stats <- c("abs_indiv_bias", "abs_effect_mod_bias",
           "abs_joint_abs_bias", "abs_intxn_bias", "1/sqrt(n)")

abs_bias_stats_plot <- make_sim_statistics_plot(
  sim_statistics_long,
  stats = stats
)

stats <- c("indiv_est_MSE", "effect_est_MSE",
           "joint_est_MSE", "intxn_est_MSE")

MSE_stats_plot <- make_sim_statistics_plot(
  sim_statistics_long,
  stats = stats
)

stats <- c("indiv_est_MSE", "effect_est_MSE",
           "joint_est_MSE", "intxn_est_MSE")

MSE_stats_plot <- make_sim_statistics_plot(
  sim_statistics_long,
  stats = stats
)

stats <- c("indiv_CI_coverage", "em_cov",
           "joint_cov", "intxn_cov")

cov_stats_plot <- make_sim_statistics_plot(
  sim_statistics_long,
  stats = stats
)

stats <- c("indiv_est_sd", "effect_mod_sd",
           "joint_est_sd", "intxn_est_sd")

sd_stats_plot <- make_sim_statistics_plot(
  sim_statistics_long,
  stats = stats
)

stats <- c("indiv_sqrt_n_abs_bias", "effect_mod_sqrt_n_abs_bias",
           "joint_sqrt_n_abs_bias", "intxn_sqrt_n_abs_bias")

sd_stats_plot <- make_sim_statistics_plot(
  sim_statistics_long,
  stats = stats
)

stats <- c("indiv_sqrt_n_abs_bias", "effect_mod_sqrt_n_abs_bias",
           "joint_sqrt_n_abs_bias", "intxn_sqrt_n_abs_bias")

sd_stats_plot <- make_sim_statistics_plot(
  sim_statistics_long,
  stats = stats
)

stats <- c("indiv_n_MSE", "em_n_MSE", "joint_n_MSE", "intxn_n_MSE")

sd_stats_plot <- make_sim_statistics_plot(
  sim_statistics_long,
  stats = stats
)



ggsave(
  here("sandbox/plots/NDE_combined.png"),
  NDE_combined_stats_plot,
  width = 13,
  height = 7,
  device = png()
)

ggsave(
  here("sandbox/plots/NIE_combined.png"),
  NIE_combined_stats_plot,
  width = 13,
  height = 7,
  device = png()
)

ggsave(
  here("sandbox/plots/NDE_steps.png"),
  NDE_steps_plot,
  width = 13,
  height = 7,
  device = png()
)

ggsave(
  here("sandbox/plots/NIE_steps.png"),
  NIE_steps_plot,
  width = 13,
  height = 7,
  device = png()
)
