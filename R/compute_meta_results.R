
compute_meta_results <- function(SuperNOVA_results, parameter) {
  if (parameter == "Indiv Shift") {
    stratified_results <- split(SuperNOVA_results$`Indiv Shift Results`, SuperNOVA_results$`Indiv Shift Results`$Variables)
  } else if (parameter == "Joint Shift") {
    stratified_results <-
      SuperNOVA_results$`Joint Shift Results` %>% dplyr::group_by(Condition, Variables)
    stratified_results <- dplyr::group_split(stratified_results)
  } else {
    stratified_results <- split(SuperNOVA_results$`Effect Mod Results`, SuperNOVA_results$Condition)
  }

  plot_list <- list()
  pooled_results_list <- list()

  if (parameter == "Joint Shift") {
    names_list <- list()
  }

  for (i in 1:length(stratified_results)) {
    results_df <- stratified_results[[i]]
    title <- paste(unique(results_df$Condition), unique(results_df$Variables), sep = "-")

    if (parameter == "Joint Shift") {
      names_list[i] <- title
    }

    weighted_mean <- sum(results_df$Psi * (1 / results_df$SE^2)) / sum((1 / results_df$SE^2))
    pooled_se <- sqrt(1 / (1 / sum(results_df$SE^2)))

    pooled_P_val <- round(2 * stats::pnorm(abs(weighted_mean / pooled_se), lower.tail = F), 5)

    pooled_CI <- c(
      round(weighted_mean + stats::qnorm(0.05 / 2, lower.tail = T) * pooled_se, 4),
      round(weighted_mean + stats::qnorm(0.05 / 2, lower.tail = F) * pooled_se, 4)
    )


    pooled_results <- c(
      unique(results_df$Condition),
      as.numeric(weighted_mean),
      as.numeric(pooled_se^2),
      pooled_se,
      pooled_CI[1],
      pooled_CI[2],
      pooled_P_val,
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
    plot_width <- 10
    plot_height <- 8
    color_vals <- c("#311885", "#E69F00", "#56B4E9")
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
      ggplot2::geom_vline(xintercept = 0, alpha = .25, linetype = "dotted", size = line_size) +
      ggplot2::labs(x = "Psi", y = "Fold", color = "") +
      ggplot2::ggtitle(title) +
      ggplot2::theme_classic() +
      ggplot2::theme(text = text_theme, axis.text = axis_text_theme, legend.position = "none") +
      ggplot2::scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

    # ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45))

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

# meta_results_Ave <- lapply(split(SuperNOVA_results, SuperNOVA_results$Variables),
#                           function(dd) (sum(dd$Psi *(1/dd$SE)))/sum((1/dd$SE)))
#
#
#
# meta_results_SE <- lapply(split(SuperNOVA_results, SuperNOVA_results$Variables),
#         function(dd) sqrt(sum(dd$SE^2 * (dd$n-1) )/(sum(dd$n-1)-nrow(dd))))
