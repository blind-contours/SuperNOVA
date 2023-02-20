library(tidyverse)
library(here)
library(ggpubr)


sim_results <- readRDS(
  here("sandbox/data/SuperNOVA_sim.rds")
)

sim_statistics <- sim_results %>%
  group_by(n_obs) %>%
  summarize(
    est_bias = mean(pooled_bias),
    est_sd = sd(pooled_bias),
    est_MSE = est_bias^2 + est_sd^2,
    CI_coverage = mean(exposure_cov)
  ) %>%
  mutate(
    abs_bias = abs(est_bias),
    sqrt_n_abs_bias = sqrt(n_obs) * abs(est_bias),
    n_MSE = n_obs * est_MSE,
  ) %>%
  ungroup()

sim_statistics_long <- sim_statistics %>%
  gather(statistic, value, -c(n_obs))

n_obs <- c(250, 500, 1000, 1500, 2000, 2500, 3000, 5000)

make_sim_statistics_plot <- function(sim_statistics_long, stats) {

  filtered_sim_stats <- sim_statistics_long %>%
    filter(
      statistic %in% stats
    )

  # Create the plot with the "cividis" color scale
  ggplot(filtered_sim_stats, aes(x=n_obs, y=value, group=statistic, color=statistic)) +
    geom_line(position = position_jitter(w=0.0, h=0.0)) +
    scale_color_cividis(discrete = TRUE) +
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

bias_stats_plot <- make_sim_statistics_plot(
  sim_statistics_long,
  stats = c("abs_bias", "est_bias")
)

MSE_stats_plot <- make_sim_statistics_plot(
  sim_statistics_long,
  stats = c("est_MSE", "n_MSE")
)

NDE_pathological_stats_plot <- make_sim_statistics_plot(
  sim_statistics_long %>% filter(sim_name %in% pathalogical),
  est_type = "NDE",
  stats = c("abs_bias", "sqrt_n_abs_bias", "n_MSE", "CI_coverage")
)

NIE_pathological_stats_plot <- make_sim_statistics_plot(
  sim_statistics_long %>% filter(sim_name %in% pathalogical),
  est_type = "NIE",
  stats = c("abs_bias", "sqrt_n_abs_bias", "n_MSE", "CI_coverage")
)

NDE_combined_stats_plot <- ggarrange(
  NDE_stats_plot, NDE_pathological_stats_plot,
  ncol = 2
) %>% annotate_figure(
  top = text_grob("NDE", size = 14)
)

NIE_combined_stats_plot <- ggarrange(
  NIE_stats_plot, NIE_pathological_stats_plot,
  ncol = 2
) %>% annotate_figure(
  top = text_grob("NIE", size = 14)
)

NDE_steps_plot <- make_sim_statistics_plot(
  sim_statistics_long,
  est_type = "NDE",
  stats = c("mean_steps", "median_steps", "sd_steps")
) + ggtitle("NDE")

NIE_steps_plot <- make_sim_statistics_plot(
  sim_statistics_long,
  est_type = "NIE",
  stats = c("mean_steps", "median_steps", "sd_steps")
) + ggtitle("NIE")

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
