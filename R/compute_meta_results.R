
compute_meta_results <- function(supernova_results, parameter) {
  if (parameter == "Indiv Shift") {
    stratified_results <- supernova_results$`Indiv Shift Results`
  } else if (parameter == "Joint Shift") {
    stratified_results <-
      supernova_results$`Joint Shift Results` %>%
      dplyr::group_by(Condition, Variables)
    stratified_results <- dplyr::group_split(stratified_results)
  } else {
    stratified_results <- split(
      supernova_results$`Effect Mod Results`,
      supernova_results$`Effect Mod Results`$Condition
    )
  }

  plot_list <- list()
  pooled_results_list <- list()

  if (parameter == "Joint Shift") {
    names_list <- list()
  }

  for (i in seq_along(stratified_results)) {
    results_df <- stratified_results[[i]]
    title <- paste(unique(results_df$Condition), unique(
      results_df$Variables
    ), sep = "-")

    if (parameter == "Joint Shift") {
      names_list[i] <- title
    }

    results_df$Psi <- round(results_df$Psi, 3)
    results_df$Variance <- round(results_df$Variance, 3)
    results_df$SE <- round(results_df$SE, 3)
    results_df$`Lower CI` <- round(results_df$`Lower CI`, 3)
    results_df$`Upper CI` <- round(results_df$`Upper CI`, 3)
    results_df$`P-value` <- round(results_df$`P-value`, 3)


    weighted_mean <- sum(results_df$Psi * (1 / results_df$SE^2)) /
      sum((1 / results_df$SE^2))
    pooled_se <- sqrt(1 / (1 / sum(results_df$SE^2)))

    pooled_p_val <- round(2 * stats::pnorm(abs(weighted_mean / pooled_se),
      lower.tail = F
    ), 5)

    pooled_ci <- c(
      round(weighted_mean + stats::qnorm(0.05 / 2, lower.tail = TRUE) *
        pooled_se, 4),
      round(weighted_mean + stats::qnorm(0.05 / 2, lower.tail = FALSE) *
        pooled_se, 4)
    )


    pooled_results <- c(
      unique(results_df$Condition),
      round(as.numeric(weighted_mean), 3),
      round(as.numeric(pooled_se^2), 3),
      round(pooled_se, 3),
      round(pooled_ci[1], 3),
      round(pooled_ci[2], 3),
      round(pooled_p_val, 3),
      "Pooled",
      unique(results_df$Type),
      unique(results_df$Variables),
      sum(results_df$N),
      1
    )

    pooled_results <- rbind(results_df, pooled_results)

    pooled_results$Psi <- round(as.numeric(pooled_results$Psi), 3)
    pooled_results$`Lower CI` <- round(as.numeric(pooled_results$`Lower CI`), 3)
    pooled_results$`Upper CI` <- round(as.numeric(pooled_results$`Upper CI`), 3)

    text_size <- 12
    line_size <- 2
    point_size <- 3
    text_theme <- ggplot2::element_text(size = text_size, color = "black")
    axis_text_theme <- ggplot2::element_text(size = text_size, color = "black")

    plot <- ggplot2::ggplot(
      pooled_results,
      ggplot2::aes(
        x = Psi, y = Fold,
        xmin = `Lower CI`, xmax = `Upper CI`,
        color = Fold
      )
    ) +
      ggplot2::geom_errorbarh(size = line_size) +
      ggplot2::geom_point(size = point_size) +
      ggplot2::geom_vline(
        xintercept = 0, alpha = .25, linetype = "dotted",
        size = line_size
      ) +
      ggplot2::labs(x = "Psi", y = "Fold", color = "") +
      ggplot2::ggtitle(title) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        text = text_theme, axis.text = axis_text_theme,
        legend.position = "none"
      ) +
      ggplot2::scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9"))

    plot_list[[i]] <- plot
    pooled_results_list[[i]] <- pooled_results
  }

  if (parameter == "Joint Shift") {
    names(plot_list) <- names_list
    names(pooled_results_list) <- names_list
  } else {
    names(plot_list) <- names(stratified_results)
    names(pooled_results_list) <- names(stratified_results)
  }


  results <- list("Plots" = plot_list, "Pooled Results" = pooled_results_list)
  return(results)
}
