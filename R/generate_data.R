# ============================================================
# Internal data-generating process for the SCENIC simulation.
# ============================================================

#' Generate one synthetic dataset (internal)
#'
#' @param n         Sample size.
#' @param n_features Number of outcome / negative-control features.
#' @param beta_Z    Direct effect of Z -> Y.
#' @param alpha_M   Effect of Z -> M (mediator path).
#' @param beta_M    Effect of M -> Y.
#' @param conf_str  Strength of U -> Z and U -> Y confounding.
#' @param w_signal  Signal fraction of U in W (0 = pure noise, 1 = perfect proxy).
#' @param seed      Optional RNG seed.
#'
#' @return A named list with elements Z, G, Y, W, U, M, synthetic_data, true_total.
#'
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

  U  <- rnorm(n)
  U2 <- rnorm(n)
  G  <- rnorm(n)

  # Exposure: G is a valid instrument (G -> Z only, not G -> Y directly)
  Z_raw <- 0.6 * G + conf_str * U + rnorm(n, 0, 0.5)
  Z     <- as.numeric(scale(Z_raw))

  # Mediator
  M <- alpha_M * Z + rnorm(n, 0, 0.05)

  # Negative-control outcomes W (n x n_features)
  W <- matrix(NA_real_, n, n_features)
  for (f in seq_len(n_features)) {
    W[, f] <- w_signal * U + (1 - w_signal) * U2 + rnorm(n, 0, 0.3)
  }
  W <- scale(W)

  # Primary outcomes Y (n x n_features) — each feature has a random U-loading
  Y <- matrix(NA_real_, n, n_features)
  for (f in seq_len(n_features)) {
    gamma_f <- runif(1, 0.4, 0.8) * conf_str
    Y[, f]  <- beta_M * M + beta_Z * Z + gamma_f * U + rnorm(n, 0, 0.2)
  }

  list(
    Z            = Z,
    G            = matrix(rep(G, n_features), n, n_features),
    Y            = Y,
    W            = W,
    U            = U,
    M            = M,
    synthetic_data = data.frame(
      fetal_sex       = rbinom(n, 1, 0.5),
      gestational_age = rnorm(n)
    ),
    true_total   = beta_Z + alpha_M * beta_M
  )
}




usethis::use_git()       # initializes git inside scenic/
usethis::use_github()    # creates the repo on GitHub and pushes
