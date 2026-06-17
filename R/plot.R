#' Diagnostic plots for a fitted `bdpaft` model
#'
#' Produces one of several diagnostic plots from a `bdpaft` fit:
#' MCMC trace plots, posterior inclusion probability summary, or
#' cluster-specific baseline survivor curves.
#'
#' @param x A fitted object of class `"bdpaft"`.
#' @param type Character; one of `"trace"`, `"pip"`, or `"survival"`.
#' @param top_pip Integer; for `type = "pip"`, the number of top features
#'   to display. Default `30`.
#' @param ... Further graphical arguments passed to plotting calls.
#'
#' @return Invisibly returns `x`.
#' @export
plot.bdpaft <- function(x, type = c("trace", "pip", "survival"),
                        top_pip = 30, ...) {
  type <- match.arg(type)
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)

  if (type == "trace") {
    graphics::par(mfrow = c(3, 1), mar = c(3, 4, 2, 1))

    graphics::plot(x$Kplus_draws, type = "l", col = "steelblue",
                   xlab = "iteration (post burn-in)",
                   ylab = "K+",
                   main = "Number of active clusters")
    graphics::abline(h = mean(x$Kplus_draws), col = "red", lty = 2)

    graphics::plot(x$alpha_draws, type = "l", col = "darkorange",
                   xlab = "iteration",
                   ylab = expression(alpha),
                   main = "DP concentration")

    if (!is.null(x$beta_draws) && nrow(x$beta_draws) > 0) {
      q <- nrow(x$beta_draws)
      cov_names <- if (length(x$cov_names) == q) x$cov_names
                   else paste0("beta", seq_len(q))
      cols <- grDevices::hcl.colors(q, palette = "Dark 3")
      ylim <- range(x$beta_draws)
      graphics::plot(x$beta_draws[1, ], type = "l", col = cols[1],
                     ylim = ylim,
                     xlab = "iteration", ylab = expression(beta),
                     main = "Covariate coefficients")
      if (q > 1) for (j in 2:q) {
        graphics::lines(x$beta_draws[j, ], col = cols[j])
      }
      graphics::legend("topright", legend = cov_names,
                       col = cols, lty = 1, bty = "n", cex = 0.8)
    }

  } else if (type == "pip") {
    graphics::par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))

    graphics::hist(x$pip, breaks = 30, col = "steelblue", border = "white",
                   xlab = "PIP",
                   main = sprintf("PIP histogram (n>0.5: %d)",
                                  sum(x$pip > 0.5)))
    graphics::abline(v = c(0.5, 0.95), col = c("red", "darkred"),
                     lty = c(2, 3))

    top_pip <- min(top_pip, length(x$pip))
    ord <- order(x$pip, decreasing = TRUE)[seq_len(top_pip)]
    feat_names <- if (!is.null(names(x$pip))) names(x$pip)[ord]
                  else paste0("f", ord)
    graphics::barplot(rev(x$pip[ord]), horiz = TRUE,
                      names.arg = rev(feat_names),
                      las = 1, col = "steelblue", border = NA,
                      xlim = c(0, 1), cex.names = 0.7,
                      main = sprintf("Top %d features by PIP", top_pip),
                      xlab = "PIP")
    graphics::abline(v = 0.5, col = "red", lty = 2)

  } else if (type == "survival") {
    z_hat <- posterior_z(x)
    Kp <- length(unique(z_hat))
    cluster_sizes <- table(z_hat)
    cols <- grDevices::hcl.colors(Kp, palette = "Dark 3")

    # Empirical survivor by cluster (if no time/status stored, use log times)
    if (!is.null(x$logC) && !is.null(x$Y)) {
      time_obs <- exp(x$logC)
      status_obs <- x$Y
    } else {
      # Try to estimate from posterior log-times
      message("Note: original time/status not stored; using posterior mean log-times.")
      time_obs <- exp(rowMeans(x$t_draws))
      status_obs <- rep(1, length(time_obs))
    }

    graphics::par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))
    tmax <- max(time_obs, na.rm = TRUE)
    graphics::plot(NA, xlim = c(0, tmax), ylim = c(0, 1),
                   xlab = "time", ylab = "survival probability",
                   main = sprintf("Cluster-specific survivor (K+ = %d)", Kp))

    for (k in seq_len(Kp)) {
      idx <- z_hat == sort(unique(z_hat))[k]
      if (sum(idx) < 2) next
      t_k <- time_obs[idx]; s_k <- status_obs[idx]
      ord <- order(t_k); t_k <- t_k[ord]; s_k <- s_k[ord]
      at_risk <- length(t_k)
      surv <- 1; t_step <- 0
      for (i in seq_along(t_k)) {
        if (s_k[i] == 1) surv <- c(surv, tail(surv, 1) *
                                    (1 - 1 / at_risk))
        else surv <- c(surv, tail(surv, 1))
        t_step <- c(t_step, t_k[i])
        at_risk <- at_risk - 1
      }
      graphics::lines(t_step, surv, col = cols[k], lwd = 2, type = "s")
    }

    graphics::legend("topright",
                     legend = sprintf("Cluster %d (n=%d)",
                                      sort(unique(z_hat)),
                                      as.numeric(cluster_sizes)),
                     col = cols, lty = 1, lwd = 2, bty = "n")
  }

  invisible(x)
}
