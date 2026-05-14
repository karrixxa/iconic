# ============================================================
# Internal loop wrappers that apply all four estimators across
# the n_features columns of a dataset returned by generate_toy_data().
# ============================================================


#' Apply all estimators across features (internal)
#'
#' @param dat        List returned by \code{generate_toy_data()}.
#' @param n_features Number of outcome columns to process.
#'
#' @return Data frame with columns: feature, method, beta, se, pvalue.
#' @keywords internal
run_methods <- function(dat, n_features = ncol(dat$Y)) {

  Z  <- dat$Z
  G  <- dat$G
  Y  <- dat$Y
  W  <- dat$W
  cv <- dat$synthetic_data

  W_avg <- rowMeans(W)

  results <- vector("list", n_features * 5L)
  j <- 1L

  for (f in seq_len(n_features)) {
    y <- Y[, f]; w <- W[, f]; g <- G[, f]
    ok <- complete.cases(cbind(y, w, g, Z, cv))
    if (sum(ok) < 20) { j <- j + 5L; next }

    y_f  <- y[ok]
    w_f  <- w[ok]
    g_f  <- g[ok]
    Z_f  <- Z[ok]
    Wa_f <- W_avg[ok]
    cv_f <- cv[ok, , drop = FALSE]

    # 0. UNADJ
    res <- tryCatch({
      fit <- lm(y_f ~ Z_f)
      sm  <- summary(fit)$coefficients
      list(beta = coef(fit)["Z_f"], se = sm["Z_f", 2], pvalue = sm["Z_f", 4])
    }, error = function(e) list(beta = NA_real_, se = NA_real_, pvalue = NA_real_))
    results[[j]] <- data.frame(feature = f, method = "UNADJ",
                               beta = as.numeric(res$beta),
                               se = as.numeric(res$se),
                               pvalue = as.numeric(res$pvalue),
                               stringsAsFactors = FALSE)
    j <- j + 1L

    # 1. DIRECT
    res <- tryCatch(
      fit_direct(y_f, Z_f, g_f, w_f, cv_f),
      error = function(e) list(beta = NA_real_, se = NA_real_, pvalue = NA_real_)
    )
    results[[j]] <- data.frame(feature = f, method = "DIRECT",
                               beta = as.numeric(res$beta),
                               se = as.numeric(res$se),
                               pvalue = as.numeric(res$pvalue),
                               stringsAsFactors = FALSE)
    j <- j + 1L

    # 2. COCA
    res <- fit_coca(y_f, Z_f, Wa_f, cv_f)
    results[[j]] <- data.frame(feature = f, method = "COCA",
                               beta  = as.numeric(res$beta),
                               se = as.numeric(res$se),
                               pvalue = as.numeric(res$pvalue),
                               stringsAsFactors = FALSE)
    j <- j + 1L

    # 3. IV2SLS
    res <- tryCatch(
      fit_iv2sls(y_f, Z_f, g_f, w_f, cv_f),
      error = function(e) list(beta = NA_real_, se = NA_real_, pvalue = NA_real_)
    )
    results[[j]] <- data.frame(feature = f, method = "IV2SLS",
                               beta = as.numeric(res$beta),
                               se = as.numeric(res$se),
                               pvalue = as.numeric(res$pvalue),
                               stringsAsFactors = FALSE)
    j <- j + 1L

    # 4. PGC
    res <- tryCatch(
      fit_pgc(y_f, Z_f, g_f, Wa_f, cv_f),
      error = function(e) list(beta = NA_real_, se = NA_real_, pvalue = NA_real_)
    )
    results[[j]] <- data.frame(feature = f, method = "PGC",
                               beta = as.numeric(res$beta),
                               se = as.numeric(res$se),
                               pvalue = as.numeric(res$pvalue),
                               stringsAsFactors = FALSE)
    j <- j + 1L
  }

  do.call(rbind, Filter(Negate(is.null), results))
}


#' Summarise simulation results across features (internal)
#'
#' @param combined   Data frame from \code{run_methods()}.
#' @param true_total Scalar true total causal effect.
#'
#' @return Data frame with one row per method: mean, median, sd, bias,
#'         abs_bias, rmse, power, n.
#' @keywords internal
summarise_results <- function(combined, true_total) {
  methods <- c("UNADJ", "DIRECT", "COCA", "IV2SLS", "PGC")

  rows <- lapply(methods, function(m) {
    x <- combined$beta[combined$method == m]
    p <- combined$pvalue[combined$method == m]
    data.frame(
      method   = m,
      mean     = mean(x, na.rm = TRUE),
      median   = median(x, na.rm = TRUE),
      sd       = sd(x, na.rm = TRUE),
      bias     = mean(x, na.rm = TRUE) - true_total,
      abs_bias = abs(mean(x, na.rm = TRUE) - true_total),
      rmse     = sqrt(mean((x - true_total)^2, na.rm = TRUE)),
      power    = mean(p < 0.05, na.rm = TRUE),
      n        = sum(!is.na(x)),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}


