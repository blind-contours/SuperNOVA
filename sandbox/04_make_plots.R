library(tidyverse)
library(here)
library(ggpubr)

if (TRUE) {
  source(here("tests/testthat/helper-dgp.R"))
  source(here("tests/testthat/helper-NDE_NIE.R"))
  sim_truth_binary_y <- get_sim_truth_NIE_NDE(n_obs = 1e7, binary_outcome = TRUE, EIC = TRUE)
  sim_truth_cont_y <- get_sim_truth_NIE_NDE(n_obs = 1e7, binary_outcome = FALSE, EIC = TRUE)

  saveRDS(sim_truth_binary_y, here("sandbox/data/sim_truth_binary_y.rds"))
  saveRDS(sim_truth_cont_y, here("sandbox/data/sim_truth_cont_y.rds"))
} else {
  sim_truth_binary_y <- readRDS(here("sandbox/data/sim_truth_binary_y.rds"))
  sim_truth_cont_y <- readRDS(here("sandbox/data/sim_truth_cont_y.rds"))
}

## set sim truth here for binary or continuous plotting, calcs
sim_truth <- sim_truth_cont_y

EY_A1_Z1 <- sim_truth$EY_A1_Z1
EY_A1_Z0 <- sim_truth$EY_A1_Z0
EY_A0_Z1 <- sim_truth$EY_A0_Z1
EY_A0_Z0 <- sim_truth$EY_A0_Z0
# compute true parameters via empirical substitution estimators
psi_NDE_true <- mean(EY_A1_Z0 - EY_A0_Z0)
psi_NIE_true <- mean(EY_A1_Z1 - EY_A1_Z0)
n_MSE_NDE_true <- var(sim_truth$EIC_NDE)
n_MSE_NIE_true <- var(sim_truth$EIC_NIE)

sim_results <- readRDS(
  here("sandbox/data/tmle3mediate_2020-05-10_21:19:18.rds")
)

# combine simulation results into one big df with n_obs column
sim_results <- names(sim_results) %>% map_dfr(
  function(this_n) {
    df <- sim_results[[this_n]]
    n_num <- as.numeric(stringr::str_replace(this_n, "n_", ""))
    cbind(n_obs = n_num, df)
  }
)

## TODO: fix labeling
rename_sim_types <- function(sym_types) {
  sapply(
    sym_types, switch,
    corr_NDE = "Efficient",
    mis_i_NDE = "Eff.(g & e--mis.)",
    mis_ii_NDE = "Eff.(psi_Z & e--mis.)",
    mis_iii_NDE = "Eff.(psi_Z & m--mis.)",
    mis_e_NDE = "Eff.(e--mis.)",
    mis_m_NDE = "Eff.(m--mis.)",
    mis_g_NDE = "Eff.(g--mis.)",
    mis_psi_Z_NDE = "Eff.(psi_Z--mis.)",
    corr_NIE = "Efficient",
    mis_i_NIE = "Eff.(g & e--mis.)",
    mis_ii_NIE = "Eff.(psi_Z & e--mis.)",
    mis_iii_NIE = "Eff.(psi_Z & m--mis.)",
    mis_e_NIE = "Eff.(e--mis.)",
    mis_m_NIE = "Eff.(m--mis.)",
    mis_g_NIE = "Eff.(g--mis.)",
    mis_psi_Z_NIE = "Eff.(psi_Z--mis.)"
  )
}

sim_statistics <- sim_results %>%
  mutate(
    covers = ifelse(
      (lower <= psi_NDE_true & psi_NDE_true <= upper) |
        (lower <= psi_NIE_true & psi_NIE_true <= upper), 1, 0
    )
  ) %>%
  group_by(n_obs, type, sim_type) %>%
  summarize(
    est_bias = if (first(type) == "NDE") {
      mean(tmle_est) - psi_NDE_true
    } else {
      mean(tmle_est) - psi_NIE_true
    },
    est_sd = sd(tmle_est),
    est_MSE = est_bias^2 + est_sd^2,
    mean_steps = mean(steps),
    median_steps = median(steps),
    sd_steps = sd(steps),
    CI_coverage = mean(covers)
  ) %>%
  mutate(
    sim_name = rename_sim_types(sim_type),
    ## stats for plotting
    abs_bias = abs(est_bias),
    sqrt_n_abs_bias = sqrt(n_obs) * abs(est_bias),
    n_MSE = n_obs * est_MSE,
  ) %>%
  ungroup()

sim_statistics_long <- sim_statistics %>%
  gather(statistic, value, -c(n_obs, type, sim_name, sim_type))

n_obs <- (cumsum(rep(sqrt(100), 8))^2)[-1]

make_sim_statistics_plot <- function(sim_statistics_long, est_type, stats) {
  filtered_sim_stats <- sim_statistics_long %>%
    filter(
      type == est_type,
      statistic %in% stats
    )

  ggplot(filtered_sim_stats, aes(n_obs, value)) +
    geom_point() +
    geom_line(linetype = "dashed") +
    facet_grid(statistic ~ sim_name, scales = "free_y") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
    ) +
    scale_x_sqrt(breaks = n_obs)
}

pathalogical <- c(
  "Eff.(psi_Z & e--mis.)",
  "Eff.(psi_Z & m--mis.)",
  "Eff.(m--mis.)",
  "Eff.(psi_Z--mis.)"
)

NDE_stats_plot <- make_sim_statistics_plot(
  sim_statistics_long %>% filter(!sim_name %in% pathalogical),
  est_type = "NDE",
  stats = c("abs_bias", "sqrt_n_abs_bias", "n_MSE", "CI_coverage")
)

NIE_stats_plot <- make_sim_statistics_plot(
  sim_statistics_long %>% filter(!sim_name %in% pathalogical),
  est_type = "NIE",
  stats = c("abs_bias", "sqrt_n_abs_bias", "n_MSE", "CI_coverage")
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
