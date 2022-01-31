joint_stoch_shift_est_Q <- function(exposures, delta, stack, covars, Av, At) {
  future::plan(future::sequential, gc = TRUE)

  # scale the outcome for logit transform
  At$y_star <- scale_to_unit(vals = At$Y)
  Av$y_star <- scale_to_unit(vals = Av$Y)

  Y_names <- "y_star"

  results <- list()

  for (i in 1:length(exposures)) {
    exposure <- exposures[[i]]

    # need a data set with the exposure stochastically shifted DOWNWARDS A-delta
    Av_downshifted <- data.table::copy(Av)
    data.table::set(Av_downshifted, j = exposure, value = shift_additive(
      A = subset(Av, select = exposure), delta = -delta
    ))

    # need a data set with the exposure stochastically shifted UPWARDS A+delta
    Av_upshifted <- data.table::copy(Av)
    data.table::set(Av_upshifted, j = exposure, value = shift_additive(
      A = subset(Av, select = exposure), delta = delta
    ))

    # need a data set with the exposure stochastically shifted UPWARDS A+2delta
    Av_upupshifted <- data.table::copy(Av)
    data.table::set(Av_upupshifted, j = exposure, value = shift_additive(
      A = subset(Av, select = exposure), delta = 2 * delta
    ))

    # Outcome mechanism
    sl <- Lrnr_sl$new(
      learners = stack
    )

    At_task_noshift <- sl3::sl3_Task$new(
      data = At,
      covariates = covars,
      outcome = Y_names
    )

    Av_task_noshift <- sl3::sl3_Task$new(
      data = Av,
      covariates = covars,
      outcome = Y_names
    )

    Av_task_upshift <- sl3::sl3_Task$new(
      data = Av_upshifted,
      covariates = covars,
      outcome = Y_names
    )

    Av_task_upupshift <- sl3::sl3_Task$new(
      data = Av_upupshifted,
      covariates = covars,
      outcome = Y_names
    )

    Av_task_downshift <- sl3::sl3_Task$new(
      data = Av_downshifted,
      covariates = covars,
      outcome = Y_names
    )

    sl_fit <- sl$train(At_task_noshift)

    # fit new Super Learner to the natural (no shift) data and predict
    Av_pred_star_Qn <- sl_fit$predict(Av_task_noshift)
    Av_pred_star_Qn <- ifelse(Av_pred_star_Qn >= 1.0, 1.0, Av_pred_star_Qn)
    Av_pred_star_Qn <- ifelse(Av_pred_star_Qn <= 0, 0, Av_pred_star_Qn)

    # predict with Super Learner from unshifted data on the shifted data
    Av_pred_star_Qn_upshifted <- sl_fit$predict(Av_task_upshift)
    Av_pred_star_Qn_upshifted <- ifelse(Av_pred_star_Qn_upshifted >= 1.0, 1.0, Av_pred_star_Qn_upshifted)
    Av_pred_star_Qn_upshifted <- ifelse(Av_pred_star_Qn_upshifted <= 0, 0, Av_pred_star_Qn_upshifted)

    # predict with Super Learner from unshifted data on the shifted data
    Av_pred_star_Qn_upupshifted <- sl_fit$predict(Av_task_upupshift)
    Av_pred_star_Qn_upupshifted <- ifelse(Av_pred_star_Qn_upupshifted >= 1.0, 1.0, Av_pred_star_Qn_upupshifted)
    Av_pred_star_Qn_upupshifted <- ifelse(Av_pred_star_Qn_upupshifted <= 0, 0, Av_pred_star_Qn_upupshifted)

    Av_pred_star_Qn_downshifted <- sl_fit$predict(Av_task_downshift)
    Av_pred_star_Qn_downshifted <- ifelse(Av_pred_star_Qn_downshifted >= 1.0, 1.0, Av_pred_star_Qn_downshifted)
    Av_pred_star_Qn_downshifted <- ifelse(Av_pred_star_Qn_downshifted <= 0, 0, Av_pred_star_Qn_downshifted)

    # create output data frame and return result
    out <- data.table::as.data.table(cbind(Av_pred_star_Qn, Av_pred_star_Qn_upshifted, Av_pred_star_Qn_upupshifted, Av_pred_star_Qn_downshifted))
    data.table::setnames(out, c("noshift", "upshift", "upupshift", "downshift"))

    results[[i]] <- out
  }

  return(results)
}
