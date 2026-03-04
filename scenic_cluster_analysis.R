################################################################################
# Run on the cluster
################################################################################

#0. Setup

N_CORES    <- 12L     # adjust to your cluster allocation
N_ITER     <- 100L    # simulations per grid point
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


#1. Experiment 1: Sweep direct effect beta_Z

cat("Experiment 1\n")
sweep1 <- sweep_param(
  param      = "beta_Z",
  param_grid = c(0, 0.05, 0.10, 0.20, 0.40),
  n_iter     = N_ITER,
  n_samples  = N_SAMPLES,
  n_features = N_FEATURES,
  alpha_M    = ALPHA_M,
  beta_M     = BETA_M,
  conf_str   = CONF_STR,
  w_signal   = 0.70,
  base_seed  = 1000L,
  n_cores    = N_CORES
)


#2. Experiment 2: Sweep confounding strength

cat("Experiment 2\n")
sweep2 <- sweep_param(
  param      = "conf_str",
  param_grid = c(0.2, 0.4, 0.6, 0.8, 1.0),
  n_iter     = N_ITER,
  n_samples  = N_SAMPLES,
  n_features = N_FEATURES,
  beta_Z     = BETA_Z,
  alpha_M    = ALPHA_M,
  beta_M     = BETA_M,
  w_signal   = 0.70,
  base_seed  = 2000L,
  n_cores    = N_CORES
)


#3. Experiment 3: Null — Type I Error

cat("Experiment 3\n")
null_res <- run_null_sim(
  n_iter     = N_ITER * 2L,   # more reps for stable rate estimate
  n_samples  = N_SAMPLES,
  n_features = N_FEATURES,
  conf_str   = CONF_STR,
  w_signal   = 0.70,
  base_seed  = 3000L,
  n_cores    = N_CORES
)
cat("\nType I error rates:\n")
print(null_res)


#4. Experiment 4: Sweep W proxy quality

cat("Experiment 4: sweep w_signal\n")
sweep4 <- sweep_param(
  param      = "w_signal",
  param_grid = c(0.2, 0.4, 0.6, 0.7, 0.8, 0.9),
  n_iter     = N_ITER,
  n_samples  = N_SAMPLES,
  n_features = N_FEATURES,
  beta_Z     = BETA_Z,
  alpha_M    = ALPHA_M,
  beta_M     = BETA_M,
  conf_str   = CONF_STR,
  base_seed  = 4000L,
  n_cores    = N_CORES
)


#5. Experiment 5: Distribution at baseline

cat("Experiment 5: estimate distributions at baseline\n")
dist_res <- run_simulation(
  n_iter     = N_ITER,
  n_samples  = N_SAMPLES,
  n_features = N_FEATURES,
  beta_Z     = BETA_Z,
  alpha_M    = ALPHA_M,
  beta_M     = BETA_M,
  conf_str   = CONF_STR,
  w_signal   = 0.70,
  base_seed  = 5000L,
  n_cores    = N_CORES
)
cat("\nBaseline summary:\n")
print(dist_res$summary[, c("method", "mean", "bias", "rmse", "power")])


#6. Save results

rds_path <- paste0("scenic_results_", timestamp, ".rds")
saveRDS(
  list(
    sweep1   = sweep1,
    sweep2   = sweep2,
    null_res = null_res,
    sweep4   = sweep4,
    dist_res = dist_res,
    params   = list(N_ITER = N_ITER, N_SAMPLES = N_SAMPLES,
                    N_FEATURES = N_FEATURES, BETA_Z = BETA_Z,
                    ALPHA_M = ALPHA_M, BETA_M = BETA_M,
                    CONF_STR = CONF_STR, TRUE_TOTAL = TRUE_TOTAL)
  ),
  file = rds_path
)
cat("\nResults saved to:", rds_path, "\n")


#7. Plots

pdf_path <- paste0("scenic_plots_", timestamp, ".pdf")
pdf(pdf_path, width = 11, height = 8.5)
par(mar = c(5, 4.5, 4, 2) + 0.1)

# Fig 1a/b — beta_Z sweep
par(mfrow = c(1, 2))
plot_estimated_vs_true(
  sweep1,
  title = "Fig 1a: Estimated vs True Effect\n(varying direct effect beta_Z)"
)
plot_bias(
  sweep1,
  param_label = "True Total Effect",
  title       = "Fig 1b: |Bias| vs True Effect\n(varying beta_Z)"
)

# Fig 2a/b — confounding sweep
par(mfrow = c(1, 2))
plot_bias(
  sweep2,
  param_label = "Confounding Strength (U -> Y)",
  title       = "Fig 2a: Bias vs Confounding Strength\n(true total = 0.25, w_signal = 0.7)"
)
plot_power(
  sweep2,
  param_label = "Confounding Strength (U -> Y)",
  title       = "Fig 2b: Detection Rate vs Confounding Strength",
  legend_pos  = "bottomleft"
)

# Fig 3 — Type I error
par(mfrow = c(1, 1))
plot_type1_error(null_res,
                  title = "Fig 3: Type I Error (true total = 0, conf_str = 0.8)")

# Fig 4a/b — proxy quality sweep
par(mfrow = c(1, 2))
plot_bias(
  sweep4,
  param_label = "W Proxy Quality (0=noise, 1=perfect U proxy)",
  title       = "Fig 4a: Bias vs Proxy Quality"
)
plot_power(
  sweep4,
  param_label = "W Proxy Quality",
  title       = "Fig 4b: Detection Rate vs Proxy Quality"
)

# Fig 5 — estimate distributions
par(mfrow = c(1, 1))
plot_estimate_distribution(dist_res)

dev.off()
cat("Plots saved to:", pdf_path, "\n")


#8. Console summary tables

cat("SUMMARY TABLES\n")

cat("Experiment 1: beta_Z sweep\n")
for (m in scenic_method_order) {
  sub <- sweep1[sweep1$method == m, c("param_value", "true_total", "mean", "bias", "rmse")]
  sub <- sub[order(sub$param_value), ]
  cat(sprintf("\n%s:\n", m))
  print(sub, row.names = FALSE, digits = 4)
}

cat("Experiment 3: Type I error\n")
print(null_res, row.names = FALSE)

cat("Experiment 4: w_signal sweep (COCA / IV2SLS / PGC only)\n")
for (m in c("COCA", "IV2SLS", "PGC")) {
  sub <- sweep4[sweep4$method == m, c("param_value", "bias", "power")]
  sub <- sub[order(sub$param_value), ]
  cat(sprintf("\n%s:\n", m))
  print(sub, row.names = FALSE, digits = 4)
}
