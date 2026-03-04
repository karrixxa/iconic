# ============================================================
# Exported plotting helpers.  All functions accept the data
# frames returned by sweep_param() or run_simulation().
# ============================================================

#' Colour palette and order for scenic methods
#' @export
scenic_method_colors <- c(
  UNADJ  = "#999999",
  DIRECT = "#E69F00",
  COCA   = "#56B4E9",
  IV2SLS = "#009E73",
  PGC    = "#CC79A7"
)

#' Default method display order
#' @export
scenic_method_order <- c("UNADJ", "DIRECT", "COCA", "IV2SLS", "PGC")


#' Plot estimated effect vs true effect (Experiment 1 style)
#'
#' @param sweep_df   Data frame returned by \code{sweep_param()}.
#' @param methods    Methods to plot. Default: all five.
#' @param title      Plot title.
#'
#' @export
plot_estimated_vs_true <- function(sweep_df,
                                   methods = scenic_method_order,
                                   title   = "Estimated vs True Effect") {
  tt  <- sort(unique(sweep_df$true_total))
  ylim <- range(sweep_df$mean[sweep_df$method %in% methods], na.rm = TRUE)
  ylim <- ylim + c(-0.05, 0.05)

  plot(NA, xlim = range(tt) + c(-0.02, 0.02), ylim = ylim,
       xlab = "True Total Effect", ylab = "Estimated Mean Effect",
       main = title, las = 1)
  abline(0, 1, lty = 2, lwd = 1.5)
  abline(h = 0, col = "grey80"); abline(v = 0, col = "grey80")

  for (m in methods) {
    sub <- sweep_df[sweep_df$method == m, ]
    sub <- sub[order(sub$true_total), ]
    lines(sub$true_total, sub$mean,
          col = scenic_method_colors[m], lwd = 2)
    points(sub$true_total, sub$mean,
           col = scenic_method_colors[m], pch = 16, cex = 1.1)
  }
  .add_legend()
}


#' Plot absolute bias vs a swept parameter
#'
#' @param sweep_df    Data frame from \code{sweep_param()}.
#' @param param_label X-axis label. Default \code{"Parameter Value"}.
#' @param methods     Methods to include. Default: all five.
#' @param title       Plot title.
#'
#' @export
plot_bias <- function(sweep_df,
                      param_label = "Parameter Value",
                      methods     = scenic_method_order,
                      title       = "Absolute Bias vs Parameter") {

  xvals <- sort(unique(sweep_df$param_value))
  ylim  <- c(0, max(sweep_df$abs_bias[sweep_df$method %in% methods],
                    na.rm = TRUE) * 1.1)

  plot(NA, xlim = range(xvals) + c(-0.02, 0.02), ylim = ylim,
       xlab = param_label, ylab = "|Bias|",
       main = title, las = 1)
  abline(h = 0, lty = 2, lwd = 1.5)

  for (m in methods) {
    sub <- sweep_df[sweep_df$method == m, ]
    sub <- sub[order(sub$param_value), ]
    lines(sub$param_value, sub$abs_bias,
          col = scenic_method_colors[m], lwd = 2)
    points(sub$param_value, sub$abs_bias,
           col = scenic_method_colors[m], pch = 16, cex = 1.1)
  }
  .add_legend("topright")
}


#' Plot detection rate (power) vs a swept parameter
#'
#' @param sweep_df    Data frame from \code{sweep_param()}.
#' @param param_label X-axis label.
#' @param methods     Methods to include.
#' @param title       Plot title.
#' @param legend_pos  Legend position string. Default \code{"bottomleft"}.
#'
#' @export
plot_power <- function(sweep_df,
                       param_label = "Parameter Value",
                       methods     = scenic_method_order,
                       title       = "Detection Rate vs Parameter",
                       legend_pos  = "bottomleft") {

  xvals <- sort(unique(sweep_df$param_value))

  plot(NA, xlim = range(xvals) + c(-0.02, 0.02), ylim = c(0, 1.05),
       xlab = param_label, ylab = "Detection Rate (p < 0.05)",
       main = title, las = 1)
  abline(h = 0.05, lty = 2, col = "red", lwd = 1.5)

  for (m in methods) {
    sub <- sweep_df[sweep_df$method == m, ]
    sub <- sub[order(sub$param_value), ]
    lines(sub$param_value, sub$power,
          col = scenic_method_colors[m], lwd = 2)
    points(sub$param_value, sub$power,
           col = scenic_method_colors[m], pch = 16, cex = 1.1)
  }
  .add_legend(legend_pos)
  mtext("Red dashed = nominal alpha = 0.05", side = 1,
        line = 3.5, cex = 0.75, col = "red")
}


#' Bar chart of Type I error rates
#'
#' @param null_df   Data frame returned by \code{run_null_sim()}.
#' @param title     Plot title.
#'
#' @export
plot_type1_error <- function(null_df,
                              title = "Type I Error Rate by Method") {
  rates <- setNames(null_df$type1_error, null_df$method)
  rates <- rates[scenic_method_order]

  bp <- barplot(rates,
                names.arg = scenic_method_order,
                col       = scenic_method_colors[scenic_method_order],
                border    = NA,
                ylim      = c(0, max(rates, na.rm = TRUE) * 1.3),
                ylab      = "Type I Error Rate (prop. p < 0.05)",
                main      = title, las = 1)
  abline(h = 0.05, lty = 2, col = "red", lwd = 2)
  text(bp, rates + max(rates, na.rm = TRUE) * 0.02,
       sprintf("%.3f", rates), cex = 0.88)
  mtext("Red line = nominal alpha = 0.05",
        side = 1, line = 3.5, cex = 0.8, col = "red")
}


#' Boxplot of estimate distributions from run_simulation()
#'
#' @param sim_result  Object returned by \code{run_simulation()}.
#' @param methods     Methods to include.
#' @param title       Plot title.
#'
#' @export
plot_estimate_distribution <- function(sim_result,
                                       methods = scenic_method_order,
                                       title   = NULL) {
  combined   <- sim_result$raw
  true_total <- sim_result$true_total

  if (is.null(title)) {
    p <- sim_result$params
    title <- sprintf(
      "Distribution of Estimates\n(true=%.2f, conf=%.1f, w_sig=%.1f)",
      true_total, p$conf_str, p$w_signal
    )
  }

  bp_data <- lapply(methods, function(m)
    combined$beta[combined$method == m])
  names(bp_data) <- methods

  boxplot(bp_data,
          col     = scenic_method_colors[methods],
          border  = "grey30",
          notch   = FALSE,
          outline = FALSE,
          las     = 1,
          ylab    = "Estimated Effect of Z on Y",
          main    = title)
  abline(h = true_total, lty = 2, col = "black", lwd = 2)
  abline(h = 0, col = "grey60", lty = 3)
  mtext(sprintf("Dashed = true total (%.2f)", true_total),
        side = 1, line = 4, cex = 0.8)

  for (i in seq_along(methods)) {
    m    <- methods[i]
    bias <- mean(combined$beta[combined$method == m], na.rm = TRUE) - true_total
    mtext(sprintf("bias\n%+.3f", bias), side = 1, at = i,
          line = 1.5, cex = 0.65,
          col = if (abs(bias) < 0.05) "darkgreen" else "firebrick")
  }
}


#Internal

.add_legend <- function(pos = "topleft") {
  legend(pos,
         legend = scenic_method_order,
         col    = scenic_method_colors[scenic_method_order],
         lwd = 2, pch = 16, bty = "n", cex = 0.85, pt.cex = 1.2)
}
