data <- NIEHS_data_1
data$Z_2 <- rnorm(nrow(data))

w <- data[, c("Z", "Z_2")]
a <- data[, c("X1", "X2", "X3", "X4", "X5", "X6", "X7")]
y <- data$Y

deltas <- list("X1" = 0.5, "X2" = 0.5, "X3" = 0.5, "X4" = 0.5,
               "X5" = 0.5, "X6" = 0.5, "X7" = 0.5)

sim_results <- SuperNOVA(
  w = w,
  a = a,
  y = y,
  delta = deltas,
  estimator = "tmle",
  fluctuation = "standard",
  n_folds = 4,
  family = "continuous",
  quantile_thresh = .25,
  verbose = TRUE,
  parallel = TRUE,
  seed = 7845300)

sim_results$`Indiv Shift Results`
full_indiv_results <- do.call(rbind, sim_results$`Indiv Shift Results`)
full_indiv_results <- full_indiv_results[,-"Type"]
full_indiv_results <- full_indiv_results[,-"Variables"]

full_joint_results <- do.call(bind_rows, sim_results$`Joint Shift Results`)
drops <- c("Type","Variables", "Delta X5", "Delta X7", "Delta X2", "Delta X1")
full_joint_results <- full_joint_results[ , !(names(full_joint_results) %in% drops)]



