indiv_stoch_shift <- function(target, delta, sl_dens_lrnr, A, W, V, Av, Av_shift, At) {
  covars <- c(W, V)

  Av_shift[, target] <- Av_shift[, target] + delta

  task_At <- make_sl3_Task(
    data = At,
    covariates = covars,
    outcome = "M1",
    outcome_type = "density"
  )

  task_Av <- make_sl3_Task(
    data = Av,
    covariates = covars,
    outcome = "M1",
    outcome_type = "density"
  )

  task_Av_shift <- make_sl3_Task(
    data = Av_shift,
    covariates = covars,
    outcome = target,
    outcome_type = "density"
  )

  sl_fit <- sl_dens_lrnr$train(task_At)

  sl_preds_At <- sl_fit$predict()
  sl_preds_Av <- sl_fit$predict(task_Av)
  sl_preds_Av_shift <- sl_fit$predict(task_Av_shift)

  # Outcome mechanism
  sl <- Lrnr_sl$new(
    learners = stack
  )

  At_task_noshift <- sl3::sl3_Task$new(
    data = At[, c(W, V, A, "y_scaled")],
    covariates = c(W, V, A),
    outcome = "y_scaled"
  )

  Av_task_noshift <- sl3::sl3_Task$new(
    data = Av[, c(W, V, A, "y_scaled")],
    covariates = c(W, V, A),
    outcome = "y_scaled"
  )

  Av_task_shift <- sl3::sl3_Task$new(
    data = Av_shift[, c(W, V, A, "y_scaled")],
    covariates = c(W, V, A),
    outcome = "y_scaled"
  )

  sl_fit <- sl$train(At_task_noshift)

  # fit new Super Learner to the natural (no shift) data and predict
  Av_pred_star_Qn <- sl_fit$predict(Av_task_noshift)

  # predict with Super Learner from unshifted data on the shifted data
  Av_pred_star_Qn_shifted <- sl_fit$predict(Av_task_shift)

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
