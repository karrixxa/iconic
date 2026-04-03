# ============================================================
# Exported plotting helpers.  All functions accept the data
# frames returned by sweep_param() or run_simulation().
# ============================================================

#' Colour palette for scenic methods
#' @export
scenic_method_colors <- c(
  UNADJ  = "#888888",
  DIRECT = "#E07B00",
  COCA   = "#3A9EC2",
  IV2SLS = "#27A062",
  PGC    = "#C455A8"
)

#' Default method display order
#' @export
scenic_method_order <- c("UNADJ", "DIRECT", "COCA", "IV2SLS", "PGC")


# ── Internal cosmetic helpers ─────────────────────────────────────────────────

.set_theme <- function(mar = c(5, 5, 4, 2) + 0.1) {
  par(mar = mar, bg = "white",
      col.axis = "grey25", col.lab = "grey15", col.main = "grey5",
      font.main = 2, cex.main = 1.05, cex.axis = 0.85, cex.lab = 0.92,
      las = 1, tck = -0.018, mgp = c(2.8, 0.55, 0))
}

.add_grid <- function() {
  grid(nx = NULL, ny = NULL, col = "grey90", lty = 1, lwd = 0.7)
}

.draw_box <- function(vals, xc, col, bw = 0.55, pt_cex = 0.38, pt_alpha = 0.55) {
  vals <- vals[!is.na(vals)]
  if (length(vals) < 3) return(invisible())
  boxplot(vals, at = xc, add = TRUE,
          col = adjustcolor(col, 0.22), border = col,
          outline = FALSE, whisklty = 1, staplelwd = 1.8, medlwd = 2.8,
          boxwex = bw, axes = FALSE)
  set.seed(42 + round(xc * 100))
  jx <- xc + runif(length(vals), -bw * 0.35, bw * 0.35)
  points(jx, vals, pch = 16, cex = pt_cex, col = adjustcolor(col, pt_alpha))
  points(xc, mean(vals), pch = 23, bg = "white", col = col, cex = 1.5, lwd = 1.8)
}

.add_legend <- function(pos = "topleft", ncol = 1) {
  legend(pos,
         legend = scenic_method_order,
         col    = scenic_method_colors[scenic_method_order],
         lwd = 2, pch = 16, bty = "n", cex = 0.85, pt.cex = 1.2,
         ncol = ncol)
}


#' Grouped boxplots of per-seed bias across a parameter sweep
#'
#' @param iter_bias   Data frame with columns iter, method, bias, pval.
#' @param param_grid  Sorted numeric vector of parameter values.
#' @param xlab        X-axis label.
#' @param ylab        Y-axis label. Default "Bias  (mean estimate - true)".
#' @param main        Plot title.
#' @param xfmt        Format string for x-axis labels. Default "%.2f".
#' @param legend_pos  Legend position. Default "topleft".
#' @param methods     Methods to include. Default: all five.
#'
#' @export
plot_bias_boxplot <- function(iter_bias,
                              param_grid,
                              xlab,
                              ylab        = "Bias  (mean estimate - true)",
                              main        = "",
                              xfmt        = "%.2f",
                              legend_pos  = "topleft",
                              methods     = scenic_method_order) {

  pvals  <- sort(unique(param_grid))
  n_m    <- length(methods)
  n_p    <- length(pvals)
  gap    <- 1.4; bw <- 0.62
  xc_grp <- seq(1, by = n_m + gap, length.out = n_p)
  offs   <- seq(-(n_m - 1) / 2, (n_m - 1) / 2) * bw
  all_b  <- iter_bias$bias[!is.na(iter_bias$bias)]
  ylim   <- range(all_b) + c(-0.04, 0.04)

  .set_theme(mar = c(6, 5.5, 4.5, 2) + 0.1)
  plot(NA,
       xlim = c(xc_grp[1] - n_m/2 - 0.8, xc_grp[n_p] + n_m/2 + 0.8),
       ylim = ylim, xaxt = "n", xlab = xlab, ylab = ylab, main = main)
  .add_grid()

  for (pi in seq_along(pvals))
    if (pi %% 2 == 0)
      rect(xc_grp[pi] - n_m/2 - 0.5, ylim[1],
           xc_grp[pi] + n_m/2 + 0.5, ylim[2],
           col = "#F4F4F4", border = NA)
  abline(h = 0, lty = 2, lwd = 1.8, col = "grey40")

  for (pi in seq_along(pvals)) {
    sub <- iter_bias[abs(iter_bias$pval - pvals[pi]) < 1e-9, ]
    for (mi in seq_along(methods))
      .draw_box(sub$bias[sub$method == methods[mi]],
                xc_grp[pi] + offs[mi],
                scenic_method_colors[methods[mi]],
                bw = bw * 0.88)
  }

  axis(1, at = xc_grp, labels = sprintf(xfmt, pvals),
       cex.axis = 0.84, col.axis = "grey25")
  legend(legend_pos,
         legend = methods,
         fill   = adjustcolor(scenic_method_colors[methods], 0.35),
         border = scenic_method_colors[methods],
         bty = "n", cex = 0.78, ncol = 2, inset = 0.01)
}


#' Baseline bias distribution (single-setting hero plot)
#'
#' @param sim_result  Object returned by run_simulation().
#' @param methods     Methods to include. Default: all five.
#' @param title       Plot title. If NULL a default is constructed.
#'
#' @export
plot_bias_distribution <- function(sim_result,
                                   methods = scenic_method_order,
                                   title   = NULL) {
  ibias      <- sim_result$iter_bias
  true_total <- sim_result$true_total
  p          <- sim_result$params

  if (is.null(title))
    title <- sprintf(
      "Bias Distribution across %d Seeds -- Baseline\ndelta = %.1f | omega = %.1f | n = %d | tau = %.2f",
      p$n_iter, p$conf_str, p$w_signal, p$n_samples, true_total)

  n_m   <- length(methods)
  ylim1 <- range(ibias$bias, na.rm = TRUE) + c(-0.05, 0.06)

  .set_theme(mar = c(7, 5.5, 4.5, 2) + 0.1)
  plot(NA, xlim = c(0.3, n_m + 0.7), ylim = ylim1,
       xaxt = "n", xlab = "", ylab = "Bias  (mean estimate - true total)",
       main = title)
  .add_grid()
  abline(h = 0, lty = 2, lwd = 2.0, col = "grey35")
  abline(h = c(-0.05, 0.05), lty = 3, lwd = 0.9, col = "grey70")
  axis(1, at = seq_along(methods), labels = methods,
       cex.axis = 1.05, font = 2, col.axis = "grey15")

  for (i in seq_along(methods)) {
    m    <- methods[i]
    vals <- ibias$bias[ibias$method == m]
    .draw_box(vals, i, scenic_method_colors[m], bw = 0.70,
              pt_cex = 0.55, pt_alpha = 0.60)
    mn  <- mean(vals, na.rm = TRUE)
    sd_ <- sd(vals, na.rm = TRUE)
    mtext(sprintf("mean %+.3f\nSD %.4f", mn, sd_),
          side = 1, at = i, line = 2.9, cex = 0.58,
          col = if (abs(mn) < 0.05) "#1a7a1a" else "#b22222")
  }
}


#' Type I error boxplot per method
#'
#' @param null_result Object returned by run_null_sim().
#' @param methods     Methods to include. Default: all five.
#' @param conf_str    Confounding strength used (for the title). Default 0.80.
#' @param alpha       Nominal significance level. Default 0.05.
#'
#' @export
plot_type1_boxplot <- function(null_result,
                               methods  = scenic_method_order,
                               conf_str = 0.80,
                               alpha    = 0.05) {

  null_combined <- null_result$raw
  n_m           <- length(methods)

  null_iter_rates <- do.call(rbind, lapply(sort(unique(null_combined$iter)), function(i) {
    sub <- null_combined[null_combined$iter == i, ]
    do.call(rbind, lapply(methods, function(m) {
      data.frame(iter = i, method = m,
                 rate = mean(sub$pvalue[sub$method == m] < alpha, na.rm = TRUE),
                 stringsAsFactors = FALSE)
    }))
  }))

  ylim2 <- c(0, max(null_iter_rates$rate, na.rm = TRUE) * 1.18)

  .set_theme(mar = c(7, 5.5, 4.5, 2) + 0.1)
  plot(NA, xlim = c(0.3, n_m + 0.7), ylim = ylim2,
       xaxt = "n", xlab = "",
       ylab = "Type I Error Rate  (prop. p < 0.05 per seed)",
       main = sprintf("Type I Error Rate by Method\ntrue total = 0 | delta = %.1f | each point = one seed",
                      conf_str))
  .add_grid()
  abline(h = alpha, lty = 2, lwd = 2.2, col = "#c0392b")
  axis(1, at = seq_along(methods), labels = methods,
       cex.axis = 1.05, font = 2, col.axis = "grey15")

  for (i in seq_along(methods)) {
    m    <- methods[i]
    vals <- null_iter_rates$rate[null_iter_rates$method == m]
    .draw_box(vals, i, scenic_method_colors[m], bw = 0.70,
              pt_cex = 0.55, pt_alpha = 0.60)
    mtext(sprintf("mean\n%.3f", mean(vals, na.rm = TRUE)),
          side = 1, at = i, line = 2.9, cex = 0.62,
          col = if (mean(vals, na.rm = TRUE) < alpha * 1.6) "#1a7a1a" else "#b22222")
  }
  mtext(sprintf("Red dashed = nominal alpha = %.2f", alpha),
        side = 1, line = 5.8, cex = 0.72, col = "#c0392b")
}


#' Type I error rate vs confounding strength
#'
#' @param t1e_df   Data frame returned by sweep_null_by_conf().
#' @param methods  Methods to plot. Default: all five.
#' @param alpha    Nominal significance level. Default 0.05.
#' @param title    Plot title.
#'
#' @export
plot_type1_vs_conf <- function(t1e_df,
                               methods = scenic_method_order,
                               alpha   = 0.05,
                               title   = "Type I Error Rate vs Confounding Strength") {

  .set_theme(mar = c(5, 5.5, 4.5, 2) + 0.1)
  ylim <- c(0, min(1, max(t1e_df$type1_error, na.rm = TRUE) * 1.15))

  plot(NA,
       xlim = range(t1e_df$conf_str) + c(-0.05, 0.05),
       ylim = ylim,
       xlab = "Confounding Strength (delta)",
       ylab = "Type I Error Rate  (prop. p < 0.05)",
       main = title)
  .add_grid()
  abline(h = alpha, lty = 2, lwd = 2.2, col = "#c0392b")

  for (m in methods) {
    sub <- t1e_df[t1e_df$method == m, ]
    sub <- sub[order(sub$conf_str), ]
    lines(sub$conf_str, sub$type1_error, col = scenic_method_colors[m], lwd = 2.2)
    points(sub$conf_str, sub$type1_error, pch = 16, cex = 1.2,
           col = scenic_method_colors[m])
  }

  legend("topleft", legend = methods,
         col = scenic_method_colors[methods], lwd = 2.2, pch = 16,
         bty = "n", cex = 0.85, inset = 0.02)
  mtext(sprintf("Red dashed = nominal alpha = %.2f", alpha),
        side = 1, line = 3.5, cex = 0.75, col = "#c0392b")
}


#' Plot estimated effect vs true effect
#'
#' @param sweep_df   Data frame from sweep_param()$summary.
#' @param methods    Methods to plot. Default: all five.
#' @param title      Plot title.
#'
#' @export
plot_estimated_vs_true <- function(sweep_df,
                                   methods = scenic_method_order,
                                   title   = "Estimated vs True Effect") {
  tt   <- sort(unique(sweep_df$true_total))
  ylim <- range(sweep_df$mean[sweep_df$method %in% methods], na.rm = TRUE) + c(-0.05, 0.05)

  .set_theme()
  plot(NA, xlim = range(tt) + c(-0.02, 0.02), ylim = ylim,
       xlab = "True Total Effect", ylab = "Estimated Mean Effect",
       main = title, las = 1)
  abline(0, 1, lty = 2, lwd = 1.5)
  abline(h = 0, col = "grey80"); abline(v = 0, col = "grey80")

  for (m in methods) {
    sub <- sweep_df[sweep_df$method == m, ]
    sub <- sub[order(sub$true_total), ]
    lines(sub$true_total, sub$mean, col = scenic_method_colors[m], lwd = 2)
    points(sub$true_total, sub$mean, col = scenic_method_colors[m], pch = 16, cex = 1.1)
  }
  .add_legend()
}


#' Plot absolute bias vs a swept parameter
#'
#' @param sweep_df    Data frame from sweep_param()$summary.
#' @param param_label X-axis label. Default "Parameter Value".
#' @param methods     Methods to include. Default: all five.
#' @param title       Plot title.
#'
#' @export
plot_bias <- function(sweep_df,
                      param_label = "Parameter Value",
                      methods     = scenic_method_order,
                      title       = "Absolute Bias vs Parameter") {
  xvals <- sort(unique(sweep_df$param_value))
  ylim  <- c(0, max(sweep_df$abs_bias[sweep_df$method %in% methods], na.rm = TRUE) * 1.1)

  .set_theme()
  plot(NA, xlim = range(xvals) + c(-0.02, 0.02), ylim = ylim,
       xlab = param_label, ylab = "|Bias|", main = title, las = 1)
  abline(h = 0, lty = 2, lwd = 1.5)

  for (m in methods) {
    sub <- sweep_df[sweep_df$method == m, ]
    sub <- sub[order(sub$param_value), ]
    lines(sub$param_value, sub$abs_bias, col = scenic_method_colors[m], lwd = 2)
    points(sub$param_value, sub$abs_bias, col = scenic_method_colors[m], pch = 16, cex = 1.1)
  }
  .add_legend("topright")
}


#' Plot detection rate (power) vs a swept parameter
#'
#' @param sweep_df    Data frame from sweep_param()$summary.
#' @param param_label X-axis label.
#' @param methods     Methods to include.
#' @param title       Plot title.
#' @param legend_pos  Legend position. Default "bottomleft".
#'
#' @export
plot_power <- function(sweep_df,
                       param_label = "Parameter Value",
                       methods     = scenic_method_order,
                       title       = "Detection Rate vs Parameter",
                       legend_pos  = "bottomleft") {
  xvals <- sort(unique(sweep_df$param_value))

  .set_theme()
  plot(NA, xlim = range(xvals) + c(-0.02, 0.02), ylim = c(0, 1.05),
       xlab = param_label, ylab = "Detection Rate (p < 0.05)",
       main = title, las = 1)
  abline(h = 0.05, lty = 2, col = "red", lwd = 1.5)

  for (m in methods) {
    sub <- sweep_df[sweep_df$method == m, ]
    sub <- sub[order(sub$param_value), ]
    lines(sub$param_value, sub$power, col = scenic_method_colors[m], lwd = 2)
    points(sub$param_value, sub$power, col = scenic_method_colors[m], pch = 16, cex = 1.1)
  }
  .add_legend(legend_pos)
  mtext("Red dashed = nominal alpha = 0.05", side = 1, line = 3.5, cex = 0.75, col = "red")
}


#' Bar chart of Type I error rates
#'
#' @param null_df   Data frame: either run_null_sim()$rates or a data frame
#'                  with columns method and type1_error.
#' @param title     Plot title.
#'
#' @export
plot_type1_error <- function(null_df,
                             title = "Type I Error Rate by Method") {
  rates <- setNames(null_df$type1_error, null_df$method)
  rates <- rates[scenic_method_order]

  .set_theme()
  bp <- barplot(rates, names.arg = scenic_method_order,
                col = scenic_method_colors[scenic_method_order],
                border = NA,
                ylim = c(0, max(rates, na.rm = TRUE) * 1.3),
                ylab = "Type I Error Rate (prop. p < 0.05)",
                main = title, las = 1)
  abline(h = 0.05, lty = 2, col = "red", lwd = 2)
  text(bp, rates + max(rates, na.rm = TRUE) * 0.02, sprintf("%.3f", rates), cex = 0.88)
  mtext("Red line = nominal alpha = 0.05", side = 1, line = 3.5, cex = 0.8, col = "red")
}


#' Boxplot of estimate distributions from run_simulation()
#'
#' @param sim_result  Object returned by run_simulation().
#' @param methods     Methods to include. Default: all five.
#' @param title       Plot title.
#'
#' @export
plot_estimate_distribution <- function(sim_result,
                                       methods = scenic_method_order,
                                       title   = NULL) {
  combined   <- sim_result$raw
  true_total <- sim_result$true_total

  if (is.null(title)) {
    p     <- sim_result$params
    title <- sprintf("Distribution of Estimates\n(true=%.2f, conf=%.1f, w_sig=%.1f)",
                     true_total, p$conf_str, p$w_signal)
  }

  bp_data <- lapply(methods, function(m) combined$beta[combined$method == m])
  names(bp_data) <- methods

  .set_theme()
  boxplot(bp_data, col = scenic_method_colors[methods], border = "grey30",
          notch = FALSE, outline = FALSE, las = 1,
          ylab = "Estimated Effect of Z on Y", main = title)
  abline(h = true_total, lty = 2, col = "black", lwd = 2)
  abline(h = 0, col = "grey60", lty = 3)
  mtext(sprintf("Dashed = true total (%.2f)", true_total), side = 1, line = 4, cex = 0.8)

  for (i in seq_along(methods)) {
    m    <- methods[i]
    bias <- mean(combined$beta[combined$method == m], na.rm = TRUE) - true_total
    mtext(sprintf("bias\n%+.3f", bias), side = 1, at = i, line = 1.5, cex = 0.65,
          col = if (abs(bias) < 0.05) "darkgreen" else "firebrick")
  }
}
