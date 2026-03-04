# ============================================================
# Exported user-facing simulation functions:
#   run_simulation()   – single parameter set, many iterations
#   sweep_param()      – grid sweep over one parameter
#   run_null_sim()     – type-I error experiment
# ============================================================


#' Run repeated simulations for a single parameter configuration
#'
#' Repeats the data-generating process \code{n_iter} times and applies all
#' five estimators (UNADJ, DIRECT, COCA, IV2SLS, PGC) to each replicate.
#' Seeds are set as \code{base_seed + iteration_index} to ensure full
#' reproducibility while varying across iterations.
#'
#' @param n_iter     Number of simulation replicates. Default 100.
#' @param n_samples  Observations per replicate. Default 500.
#' @param n_features Number of outcome / NC features per replicate. Default 20.
#' @param beta_Z     Direct effect of Z -> Y. Default 0.10.
#' @param alpha_M    Effect of Z -> Mediator. Default 0.50.
#' @param beta_M     Effect of Mediator -> Y. Default 0.30.
#' @param conf_str   Confounding strength (U -> Z and U -> Y). Default 0.80.
#' @param w_signal   Proxy quality of W (0 = noise, 1 = perfect U proxy). Default 0.70.
#' @param base_seed  Starting seed; replicate i uses \code{base_seed + i}. Default 100.
#' @param n_cores    Number of parallel workers (uses \pkg{parallel}). Default 1.
#'
#' @return A list with:
#' \describe{
#'   \item{raw}{Data frame of per-feature per-iteration estimates (beta, se, pvalue).}
#'   \item{summary}{Per-method summary: mean, median, sd, bias, rmse, power.}
#'   \item{true_total}{Scalar true total causal effect used.}
#'   \item{params}{Named list of the simulation parameters.}
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' res <- run_simulation(n_iter = 50, beta_Z = 0.1, conf_str = 0.8, n_cores = 4)
#' res$summary
#' }
run_simulation <- function(n_iter     = 100,
                           n_samples  = 500,
                           n_features = 20,
                           beta_Z     = 0.10,
                           alpha_M    = 0.50,
                           beta_M     = 0.30,
                           conf_str   = 0.80,
                           w_signal   = 0.70,
                           base_seed  = 100,
                           n_cores    = 1) {

  true_total <- beta_Z + alpha_M * beta_M

  worker <- function(i) {
    dat <- generate_toy_data(
      n          = n_samples,
      n_features = n_features,
      beta_Z     = beta_Z,
      alpha_M    = alpha_M,
      beta_M     = beta_M,
      conf_str   = conf_str,
      w_signal   = w_signal,
      seed       = base_seed + i
    )
    res       <- run_methods(dat, n_features)
    res$iter  <- i
    res
  }

  iter_results <- .parallel_lapply(seq_len(n_iter), worker, n_cores)
  combined     <- do.call(rbind, iter_results)
  smry         <- summarise_results(combined, true_total)

  list(
    raw        = combined,
    summary    = smry,
    true_total = true_total,
    params     = list(
      n_iter     = n_iter,
      n_samples  = n_samples,
      n_features = n_features,
      beta_Z     = beta_Z,
      alpha_M    = alpha_M,
      beta_M     = beta_M,
      conf_str   = conf_str,
      w_signal   = w_signal
    )
  )
}


#' Sweep a single simulation parameter across a grid
#'
#' For each value in \code{param_grid}, runs \code{n_iter} full simulation
#' replicates and returns a tidy summary data frame.  Iterations are
#' parallelised across \code{n_cores} workers; seeds vary by iteration AND
#' grid point to avoid overlap.
#'
#' @param param       Character name of the parameter to vary. Must be one of
#'   \code{"beta_Z"}, \code{"conf_str"}, \code{"w_signal"},
#'   \code{"alpha_M"}, \code{"beta_M"}, \code{"n_samples"}.
#' @param param_grid  Numeric vector of values to sweep.
#' @param n_iter      Replicates per grid point. Default 100.
#' @param n_samples   Observations per replicate. Default 500.
#' @param n_features  Features per replicate. Default 20.
#' @param beta_Z      Baseline direct effect (overridden when \code{param = "beta_Z"}). Default 0.10.
#' @param alpha_M     Baseline mediator path. Default 0.50.
#' @param beta_M      Baseline mediator effect. Default 0.30.
#' @param conf_str    Baseline confounding strength. Default 0.80.
#' @param w_signal    Baseline proxy quality. Default 0.70.
#' @param base_seed   Seed offset (grid index * 1000 + iter). Default 0.
#' @param n_cores     Parallel workers. Default 1.
#'
#' @return A data frame with columns: \code{param_value}, \code{true_total},
#'   \code{method}, \code{mean}, \code{median}, \code{sd}, \code{bias},
#'   \code{abs_bias}, \code{rmse}, \code{power}, \code{n}.
#' @export
#'
#' @examples
#' \dontrun{
#' df <- sweep_param("conf_str", c(0.2, 0.5, 0.8, 1.0), n_iter = 50, n_cores = 4)
#' }
sweep_param <- function(param,
                        param_grid,
                        n_iter     = 100,
                        n_samples  = 500,
                        n_features = 20,
                        beta_Z     = 0.10,
                        alpha_M    = 0.50,
                        beta_M     = 0.30,
                        conf_str   = 0.80,
                        w_signal   = 0.70,
                        base_seed  = 0,
                        n_cores    = 1) {

  allowed <- c("beta_Z", "conf_str", "w_signal", "alpha_M", "beta_M", "n_samples")
  param   <- match.arg(param, allowed)

  # Build baseline arg list
  base_args <- list(
    n_samples  = n_samples,
    n_features = n_features,
    beta_Z     = beta_Z,
    alpha_M    = alpha_M,
    beta_M     = beta_M,
    conf_str   = conf_str,
    w_signal   = w_signal
  )

  grid_results <- lapply(seq_along(param_grid), function(gi) {
    pval <- param_grid[gi]
    args <- base_args
    args[[param]] <- pval

    true_total <- args$beta_Z + args$alpha_M * args$beta_M

    worker <- function(i) {
      dat <- do.call(generate_toy_data,
                     c(args, list(seed = base_seed + gi * 1000L + i)))
      res      <- run_methods(dat, args$n_features)
      res$iter <- i
      res
    }

    iter_results <- .parallel_lapply(seq_len(n_iter), worker, n_cores)
    combined     <- do.call(rbind, iter_results)
    smry         <- summarise_results(combined, true_total)

    smry$param_value <- pval
    smry$true_total  <- true_total
    smry
  })

  result           <- do.call(rbind, grid_results)
  result$param     <- param
  result[, c("param", "param_value", "true_total",
             setdiff(names(result), c("param", "param_value", "true_total")))]
}


#' Run null simulations to estimate Type I error rates
#'
#' Sets all causal effects to zero (\code{beta_Z = alpha_M = beta_M = 0})
#' and estimates the false-positive rate for each method.
#'
#' @param n_iter     Number of replicates. Default 200.
#' @param n_samples  Observations per replicate. Default 500.
#' @param n_features Features per replicate. Default 20.
#' @param conf_str   Confounding strength (U is present but has no causal path to Y). Default 0.80.
#' @param w_signal   Proxy quality of W. Default 0.70.
#' @param base_seed  Seed offset. Default 300.
#' @param n_cores    Parallel workers. Default 1.
#' @param alpha      Significance threshold. Default 0.05.
#'
#' @return A data frame with columns \code{method} and \code{type1_error}.
#' @export
run_null_sim <- function(n_iter     = 200,
                         n_samples  = 500,
                         n_features = 20,
                         conf_str   = 0.80,
                         w_signal   = 0.70,
                         base_seed  = 300,
                         n_cores    = 1,
                         alpha      = 0.05) {

  worker <- function(i) {
    dat <- generate_toy_data(
      n          = n_samples,
      n_features = n_features,
      beta_Z     = 0,
      alpha_M    = 0,
      beta_M     = 0,
      conf_str   = conf_str,
      w_signal   = w_signal,
      seed       = base_seed + i
    )
    res      <- run_methods(dat, n_features)
    res$iter <- i
    res
  }

  iter_results  <- .parallel_lapply(seq_len(n_iter), worker, n_cores)
  null_combined <- do.call(rbind, iter_results)

  methods <- c("UNADJ", "DIRECT", "COCA", "IV2SLS", "PGC")
  rates   <- sapply(methods, function(m) {
    p <- null_combined$pvalue[null_combined$method == m]
    mean(p < alpha, na.rm = TRUE)
  })

  data.frame(
    method      = methods,
    type1_error = rates,
    flag        = ifelse(rates > 3 * alpha, "INFLATED",
                  ifelse(rates < alpha / 2.5, "CONSERVATIVE", "OK")),
    stringsAsFactors = FALSE
  )
}


#Internal parallelisation helper

#' Parallel lapply dispatcher (internal)
#'
#' Uses \code{parallel::mclapply} on Unix/macOS when \code{n_cores > 1}
#' and falls back to plain \code{lapply} on Windows or when
#' \code{n_cores == 1}.
#'
#' @keywords internal
.parallel_lapply <- function(X, FUN, n_cores = 1) {
  if (n_cores > 1 && .Platform$OS.type != "windows") {
    parallel::mclapply(X, FUN, mc.cores = n_cores)
  } else {
    lapply(X, FUN)
  }
}
