# data generating mechanism function factory
make_dgp <- function() {
  # treatment mechanism
  g_mech <- function(w) {
    w <- as.data.frame(w)
    p_score <- (rowSums(w) / 4) + 0.1
    return(p_score)
  }

  # mediation mechanism
  z_mech <- function(w, a) {
    w <- as.data.frame(w)
    z1_prob <- 1 - plogis((a + w[[1]]) / (a + w[[1]]^3 + 0.5))
    z2_prob <- plogis((a - 1) + w[[2]] / (w[[3]] + 3))
    z3_prob <- plogis((a - 1) + 2 * w[[1]]^3 - 1 / (2 * w[[1]] + 0.5))
    return(list(z1_prob, z2_prob, z3_prob))
  }

  # outcome mechanism
  m_mech <- function(w, a, z, eps_sd = 0.5) {
    w <- as.data.frame(w)
    z <- as.data.frame(z)
    y <- z[[1]] + z[[2]] - z[[3]] + exp(a + z[[3]] / (1 + rowSums(w)^2))
    return(y)
  }

  # return DGP functions
  return(list(g_mech = g_mech, z_mech = z_mech, m_mech = m_mech))
}
