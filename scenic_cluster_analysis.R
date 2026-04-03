################################################################################
# Run on the cluster
################################################################################

#0. Setup

N_CORES    <- 12L
N_ITER     <- 100L
N_SAMPLES  <- 500L
N_FEATURES <- 20L
ALPHA_M    <- 0.5
BETA_M     <- 0.3
BETA_Z     <- 0.10
CONF_STR   <- 0.80
TRUE_TOTAL <- BETA_Z + ALPHA_M * BETA_M

cat("True total effect of Z on Y:", TRUE_TOTAL, "\n")
cat("Using", N_CORES, "cores,", N_ITER, "iterations per condition\n\n")

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
if (!dir.exists("figs")) dir.create("figs")


#1. Experiment 1: Sweep direct effect beta_Z

cat("Experiment 1\n")
sweep1 <- sweep_param(
  param      = "beta_Z",
  param_grid = c(0, 0.05, 0.10, 0.20, 0.40),
  n_iter = N_ITER, n_samples = N_SAMPLES, n_features = N_FEATURES,
  alpha_M = ALPHA_M, beta_M = BETA_M, conf_str = CONF_STR, w_signal = 0.70,
  base_seed = 1000L, n_cores = N_CORES
)


#2. Experiment 2: Sweep confounding strength

cat("Experiment 2\n")
sweep2 <- sweep_param(
  param      = "conf_str",
  param_grid = c(0.2, 0.4, 0.6, 0.8, 1.0),
  n_iter = N_ITER, n_samples = N_SAMPLES, n_features = N_FEATURES,
  beta_Z = BETA_Z, alpha_M = ALPHA_M, beta_M = BETA_M, w_signal = 0.70,
  base_seed = 2000L, n_cores = N_CORES
)


#3. Experiment 3: Null — Type I Error

cat("Experiment 3\n")
null_res <- run_null_sim(
  n_iter = N_ITER * 2L, n_samples = N_SAMPLES, n_features = N_FEATURES,
  conf_str = CONF_STR, w_signal = 0.70,
  base_seed = 3000L, n_cores = N_CORES
)
cat("\nType I error rates:\n")
print(null_res$rates)


#4. Experiment 4: Sweep W proxy quality

cat("Experiment 4\n")
sweep4 <- sweep_param(
  param      = "w_signal",
  param_grid = c(0.2, 0.4, 0.6, 0.7, 0.8, 0.9),
  n_iter = N_ITER, n_samples = N_SAMPLES, n_features = N_FEATURES,
  beta_Z = BETA_Z, alpha_M = ALPHA_M, beta_M = BETA_M, conf_str = CONF_STR,
  base_seed = 4000L, n_cores = N_CORES
)


#5. Experiment 5: Sweep sample size

cat("Experiment 5\n")
sweep5 <- sweep_param(
  param      = "n_samples",
  param_grid = c(100, 200, 500, 1000, 2000),
  n_iter = N_ITER, n_features = N_FEATURES,
  beta_Z = BETA_Z, alpha_M = ALPHA_M, beta_M = BETA_M,
  conf_str = CONF_STR, w_signal = 0.70,
  base_seed = 5000L, n_cores = N_CORES
)


#6. Experiment 6: Sweep mediator strength alpha_M

cat("Experiment 6\n")
sweep6 <- sweep_param(
  param      = "alpha_M",
  param_grid = c(0.0, 0.2, 0.5, 0.8, 1.0),
  n_iter = N_ITER, n_samples = N_SAMPLES, n_features = N_FEATURES,
  beta_Z = BETA_Z, beta_M = BETA_M, conf_str = CONF_STR, w_signal = 0.70,
  base_seed = 6000L, n_cores = N_CORES
)


#7. Experiment 7: Baseline 100-seed variance

cat("Experiment 7\n")
baseline_res <- run_simulation(
  n_iter = N_ITER, n_samples = N_SAMPLES, n_features = N_FEATURES,
  beta_Z = BETA_Z, alpha_M = ALPHA_M, beta_M = BETA_M,
  conf_str = CONF_STR, w_signal = 0.70,
  base_seed = 7000L, n_cores = N_CORES
)
cat("\nBaseline summary:\n")
print(baseline_res$summary[, c("method", "mean", "bias", "rmse", "power")])


#8. Experiment 8: Type I error vs confounding strength

cat("Experiment 8\n")
t1e_sweep <- sweep_null_by_conf(
  conf_grid  = c(0.2, 0.4, 0.6, 0.8, 1.0),
  n_iter = N_ITER, n_samples = N_SAMPLES, n_features = N_FEATURES,
  w_signal = 0.70, base_seed = 8000L, n_cores = N_CORES
)


#9. Save results

rds_path <- paste0("scenic_results_", timestamp, ".rds")
saveRDS(
  list(
    sweep1       = sweep1,
    sweep2       = sweep2,
    null_res     = null_res,
    sweep4       = sweep4,
    sweep5       = sweep5,
    sweep6       = sweep6,
    baseline_res = baseline_res,
    t1e_sweep    = t1e_sweep,
    params = list(N_ITER = N_ITER, N_SAMPLES = N_SAMPLES,
                  N_FEATURES = N_FEATURES, BETA_Z = BETA_Z,
                  ALPHA_M = ALPHA_M, BETA_M = BETA_M,
                  CONF_STR = CONF_STR, TRUE_TOTAL = TRUE_TOTAL)
  ),
  file = rds_path
)
cat("\nResults saved to:", rds_path, "\n")


#10. Save figures as individual PNGs into figs/

cat("\nSaving figures...\n")

save_fig <- function(num, width = 1800, height = 1400, expr) {
  path <- sprintf("figs/fig_page-%d.png", num)
  png(path, width = width, height = height, res = 180)
  expr
  dev.off()
  cat(sprintf("  Saved: %s\n", path))
}

save_fig(1, expr = plot_bias_distribution(baseline_res))

save_fig(2, expr = plot_type1_boxplot(null_res, conf_str = CONF_STR))

save_fig(3, expr =
           plot_bias_boxplot(sweep2$iter_bias, c(0.2, 0.4, 0.6, 0.8, 1.0),
                             xlab = "Confounding Strength (delta)",
                             main = "Bias vs Confounding Strength\ntau = 0.25 | omega = 0.7 | each box = 100 seeds",
                             xfmt = "%.1f"))

save_fig(4, expr =
           plot_bias_boxplot(sweep4$iter_bias, c(0.2, 0.4, 0.6, 0.7, 0.8, 0.9),
                             xlab = "Negative Control Proxy Quality (omega)",
                             main = "Bias vs Negative Control Quality\ntau = 0.25 | delta = 0.8 | each box = 100 seeds",
                             xfmt = "%.1f", legend_pos = "topright"))

save_fig(5, expr =
           plot_bias_boxplot(sweep5$iter_bias, c(100, 200, 500, 1000, 2000),
                             xlab = "Sample Size (n)",
                             main = "Bias vs Sample Size\ntau = 0.25 | delta = 0.8 | SCENIC uses N = 200",
                             xfmt = "%d", legend_pos = "topright"))

save_fig(6, width = 2400, expr = {
  par(mfrow = c(1, 2))
  plot_bias_boxplot(sweep1$iter_bias, c(0, 0.05, 0.10, 0.20, 0.40),
                    xlab = "Direct Effect (beta_Z)",
                    ylab = "Bias  (mean estimate - true)",
                    main = "Bias vs Direct Effect (beta_Z)\ntrue total = beta_Z + 0.15",
                    xfmt = "%.2f", legend_pos = "topleft")
  plot_bias_boxplot(sweep6$iter_bias, c(0.0, 0.2, 0.5, 0.8, 1.0),
                    xlab = "Mediator Strength (alpha_M)",
                    ylab = "",
                    main = "Bias vs Mediator Strength (alpha_M)\ntrue total = 0.10 + alpha_M x 0.30",
                    xfmt = "%.1f", legend_pos = "topright")
})

save_fig(7, expr = plot_type1_vs_conf(t1e_sweep))

cat("\nAll figures saved to figs/\n")


#11. Console summary tables

cat("\nType I error (null):\n")
print(null_res$rates, row.names = FALSE)

cat("\nBaseline bias variance:\n")
ib <- baseline_res$iter_bias
smry_table <- do.call(rbind, lapply(scenic_method_order, function(m) {
  b <- ib$bias[ib$method == m]
  data.frame(method    = m,
             mean_bias = round(mean(b, na.rm = TRUE), 4),
             sd_bias   = round(sd(b,   na.rm = TRUE), 4),
             min_bias  = round(min(b,  na.rm = TRUE), 4),
             max_bias  = round(max(b,  na.rm = TRUE), 4))
}))
print(smry_table, row.names = FALSE)

cat("\nType I error vs confounding strength:\n")
print(t1e_sweep, row.names = FALSE, digits = 3)
