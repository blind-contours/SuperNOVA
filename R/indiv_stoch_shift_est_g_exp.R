indiv_stoch_shift_est_g_exp <- function(exposure, delta, sl_density_lrnr, covars, Av, At) {
  future::plan(future::sequential, gc = TRUE)

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
    covariates = covars
  )

  sl_task_noshift <- sl3::sl3_Task$new(
    data = Av,
    outcome = exposure,
    covariates = covars
  )

  # sl3 task for data with exposure shifted DOWNWARDS A-delta
  sl_task_downshifted <- sl3::sl3_Task$new(
    data = Av_downshifted,
    outcome = exposure,
    covariates = covars
  )

  # sl3 task for data with exposure shifted UPWARDS A+delta
  sl_task_upshifted <- sl3::sl3_Task$new(
    data = Av_upshifted,
    outcome = exposure,
    covariates = covars
  )

  # sl3 task for data with exposure shifted UPWARDS A+2delta
  sl_task_upupshifted <- sl3::sl3_Task$new(
    data = Av_upupshifted,
    outcome = exposure,
    covariates = covars
  )

  sl_fit <- sl_density_lrnr$train(sl_task)

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
  return(out)
}
