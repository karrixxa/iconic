# ============================================================
# Internal data-generating process for the SCENIC toy simulation.
# ============================================================

#' Generate one synthetic dataset (internal)
#'
#' @param n          Sample size. Default 500.
#' @param n_features Number of outcome and negative-control features. Default 20.
#' @param beta_Z     Direct effect of Z on Y. Default 0.10.
#' @param alpha_M    Effect of Z on mediator M. Default 0.50.
#' @param beta_M     Effect of M on Y. Default 0.30.
#' @param conf_str   Confounding strength delta. Default 0.80.
#' @param w_signal   Proxy quality omega (0 = noise, 1 = perfect U proxy). Default 0.70.
#' @param seed       Optional integer RNG seed for reproducibility.
#'
#' @return A named list with elements Z, G, Y, W, U1, M, synthetic_data, true_total.
#' @keywords internal
generate_toy_data <- function(n          = 500,
                              n_features = 20,
                              beta_Z     = 0.10,
                              alpha_M    = 0.50,
                              beta_M     = 0.30,
                              conf_str   = 0.80,
                              w_signal   = 0.70,
                              seed       = NULL) {
  if (!is.null(seed)) set.seed(seed)

  U1 <- rnorm(n)
  U2 <- rnorm(n)
  G  <- rnorm(n)

  Z_raw <- 0.6 * G + conf_str * U1 + rnorm(n, 0, 0.5)
  Z <- as.numeric(scale(Z_raw))

  M <- alpha_M * Z + rnorm(n, 0, 0.05)

  W <- matrix(NA_real_, n, n_features)
  for (f in seq_len(n_features))
    W[, f] <- w_signal * U1 + (1 - w_signal) * U2 + rnorm(n, 0, 0.3)
  W <- scale(W)

  Y <- matrix(NA_real_, n, n_features)
  for (f in seq_len(n_features)) {
    gamma_f <- runif(1, 0.4, 0.8) * conf_str
    Y[, f]  <- beta_M * M + beta_Z * Z + gamma_f * U1 + rnorm(n, 0, 0.2)
  }

  list(
    Z = Z,
    G = matrix(rep(G, n_features), n, n_features),
    Y = Y,
    W = W,
    U1 = U1,
    M = M,
    synthetic_data = data.frame(
      fetal_sex = rbinom(n, 1, 0.5),
      gestational_age = rnorm(n)
    ),
    true_total = beta_Z + alpha_M * beta_M
  )
}
