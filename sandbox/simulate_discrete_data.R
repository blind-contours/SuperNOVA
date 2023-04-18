simulate_discrete_mediation_data <- function(n_obs = 100000, delta = 1) {
  # simulate some baseline covariates
  w_1 <- rnorm(n_obs, mean = 20, sd = 2)
  w_2 <- rbinom(n_obs, 1, 0.5)

  covars <- data.frame(w_1, w_2)

  # probabilities
  b0i <- round(rnorm(4, 0.3, 0.01), 2)
  b1i <- round(rnorm(4, 0.4, 0.01), 2)
  b2i <- round(rnorm(4, 0.5, 0.01), 2)

  gen_probs <- function(index, b0i, b1i, b2i, w_1, w_2) {
    denominator <- 1 + exp(b0i[index] + (b1i[index] * w_1) + (b2i[index] * w_2))
    exp(b0i[index] + (b1i[index] * w_1) + (b2i[index] * w_2)) / denominator
  }

  # Generate a_1_quant and a_1_quant_shift
  a_1_quant <- sapply(seq(n_obs), function(i) {
    w_1_i <- w_1[i]
    w_2_i <- w_2[i]
    probs <- sapply(seq(4), gen_probs, b0i, b1i, b2i, w_1_i, w_2_i)
    sample(1:4, 1, prob = probs)
  })
  a_1_quant_shift <- ifelse(a_1_quant + delta <= max(a_1_quant), a_1_quant + delta, max(a_1_quant))

  covars_exposure <- covars
  covars_exposure$a_1_quant <- a_1_quant
  covars_exposure$a_1_quant_shift <- a_1_quant_shift

  # Generate z_1_quant and z_1_quant_shift
  z_1_quant <- sapply(seq(n_obs), function(i) {
    w_1_i <- w_1[i]
    a_1_quant_i <- a_1_quant[i]
    probs <- sapply(seq(4), gen_probs, b0i, b1i, b2i, w_1_i, a_1_quant_i)
    sample(1:4, 1, prob = probs)
  })
  z_1_quant_shift <- ifelse(z_1_quant + delta <= max(z_1_quant), z_1_quant + delta, max(z_1_quant))

  # Calculate outcome variables
  y_quant <- 10 * z_1_quant + 40 * a_1_quant
  y_shift_a_1_quant <- 10 * z_1_quant + 40 * a_1_quant_shift
  y_shift_a_1_z_1_quant <- 10 * z_1_quant_shift + 40 * a_1_quant_shift

  # Calculate direct, indirect, and total effects
  nde_a1_quant <- mean(y_shift_a_1_quant - y_quant)
  nie_a1_quant <- mean(y_shift_a_1_z_1_quant - y_shift_a_1_quant)
  ate_a1_quant <- nde_a1_quant + nie_a1_quant

  data <- as.data.frame(cbind(
    w_1, w_2,
    a_1_quant, a_1_quant_shift,
    z_1_quant, z_1_quant_shift,
    y_quant, y_shift_a_1_quant, y_shift_a_1_z_1_quant
  ))

  return(list(
    "data" = data,
    "nde_a1_quant" = nde_a1_quant,
    "nie_a1_quant" = nie_a1_quant,
    "ate_a1_quant" = ate_a1_quant
  ))
}
