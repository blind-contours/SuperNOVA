joint_stoch_shift_est_g_exp <- function(exposures, delta, Density_stack, covars, Av, At) {
  future::plan(future::sequential, gc = TRUE)

  results <- list()

  for (i in 1:length(exposures)) {
    if (i == 3) {
      exposure <- exposures[[i]][2]
      covariates <- c(covars, exposures[[1]])
    } else {
      exposure <- exposures[[i]]
      covariates <- covars
    }

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

    sl_task <- sl3::sl3_Task$new(
      data = At,
      outcome = exposure,
      covariates = covariates
    )

    sl_task_noshift <- sl3::sl3_Task$new(
      data = Av,
      outcome = exposure,
      covariates = covariates
    )

    # sl3 task for data with exposure shifted DOWNWARDS A-delta
    sl_task_downshifted <- sl3::sl3_Task$new(
      data = Av_downshifted,
      outcome = exposure,
      covariates = covariates
    )

    # sl3 task for data with exposure shifted UPWARDS A+delta
    sl_task_upshifted <- sl3::sl3_Task$new(
      data = Av_upshifted,
      outcome = exposure,
      covariates = covariates
    )

    # sl3 task for data with exposure shifted UPWARDS A+2delta
    sl_task_upupshifted <- sl3::sl3_Task$new(
      data = Av_upupshifted,
      outcome = exposure,
      covariates = covariates
    )

    sl_fit <- Density_stack$train(sl_task)

    pred_g_exp_noshift <- sl_fit$predict(sl_task_noshift)
    pred_g_exp_downshifted <- sl_fit$predict(sl_task_downshifted)
    pred_g_exp_upshifted <- sl_fit$predict(sl_task_upshifted)
    pred_g_exp_upupshifted <- sl_fit$predict(sl_task_upupshifted)

    # create output data.tables
    out <- data.table::as.data.table(cbind(
      pred_g_exp_downshifted,
      pred_g_exp_noshift,
      pred_g_exp_upshifted,
      pred_g_exp_upupshifted
    ))
    data.table::setnames(out, c("downshift", "noshift", "upshift", "upupshift"))

    results[[i]] <- out
  }


  return(results)
}
