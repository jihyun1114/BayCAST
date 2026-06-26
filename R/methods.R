#' @export
print.baycast <- function(x, ...) {
  cat("BayCAST fit\n")
  cat(sprintf("  n = %d  p = %d  q = %d  d = %d\n",
              x$dims["n"], x$dims["p"], x$dims["q"], x$dims["d"]))
  cat(sprintf("  iters = %d  burn = %d  thin = %d  (kept = %d)\n",
              x$iters, x$burn, x$thin, x$kept))
  cat(sprintf("  slice = %s   init_A = %s   init_z = %s\n",
              x$slice, x$init_A, x$init_z))
  cat(sprintf("  pi0 = %.3g   sigmaA2 = %.3g\n", x$pi0, x$sigmaA2))
  cat(sprintf("  K+ posterior mode = %d   PIP > 0.5: %d / %d\n",
              .Kplus_mode(x$Kplus_draws),
              sum(x$pip > 0.5), length(x$pip)))
  invisible(x)
}

#' @export
summary.baycast <- function(object, pip_threshold = 0.5, ...) {
  Kpl  <- object$Kplus_draws
  z_last <- object$z_draws[, ncol(object$z_draws)]
  ctab  <- sort(as.integer(table(z_last)), decreasing = TRUE)

  beta_summary <- NULL
  if (!is.null(object$beta_draws)) {
    bm <- rowMeans(object$beta_draws)
    bq <- t(apply(object$beta_draws, 1, stats::quantile,
                  probs = c(0.025, 0.5, 0.975)))
    beta_summary <- data.frame(
      covariate = object$cov_names,
      mean      = round(bm, 4),
      `2.5%`    = round(bq[, 1], 4),
      `50%`     = round(bq[, 2], 4),
      `97.5%`   = round(bq[, 3], 4),
      check.names = FALSE
    )
  }

  out <- list(
    dims          = object$dims,
    kept          = object$kept,
    slice         = object$slice,
    Kplus_mode    = .Kplus_mode(Kpl),
    Kplus_mean    = round(mean(Kpl), 2),
    Kplus_q       = stats::quantile(Kpl, c(0.025, 0.5, 0.975)),
    cluster_sizes = ctab,
    n_PIP_gt      = sum(object$pip > pip_threshold),
    pip_threshold = pip_threshold,
    alpha_mean    = round(mean(object$alpha_draws), 3),
    alpha_q       = stats::quantile(object$alpha_draws, c(0.025, 0.5, 0.975)),
    beta_summary  = beta_summary
  )
  class(out) <- "summary.baycast"
  out
}

#' @export
print.summary.baycast <- function(x, ...) {
  cat("BayCAST posterior summary\n")
  cat(sprintf("  n=%d, p=%d, q=%d, d=%d   kept draws: %d   slice=%s\n",
              x$dims["n"], x$dims["p"], x$dims["q"], x$dims["d"],
              x$kept, x$slice))

  cat("\nCluster structure:\n")
  cat(sprintf("  K+ posterior   mode = %d, mean = %.2f, 95%% CI = [%d, %d]\n",
              x$Kplus_mode, x$Kplus_mean,
              as.integer(x$Kplus_q[1]), as.integer(x$Kplus_q[3])))
  cat(sprintf("  Last-iter cluster sizes (desc): %s\n",
              paste(x$cluster_sizes, collapse = ", ")))

  cat("\nVariable selection:\n")
  cat(sprintf("  PIP > %.2f: %d features (out of %d)\n",
              x$pip_threshold, x$n_PIP_gt, x$dims["p"]))

  cat("\nDP concentration:\n")
  cat(sprintf("  alpha   mean = %.3f, 95%% CI = [%.3f, %.3f]\n",
              x$alpha_mean, x$alpha_q[1], x$alpha_q[3]))

  if (!is.null(x$beta_summary)) {
    cat("\nCovariate coefficients:\n")
    print(x$beta_summary, row.names = FALSE)
  }
  invisible(x)
}

# internal helper
.Kplus_mode <- function(Kpl) {
  tab <- table(Kpl)
  as.integer(names(sort(tab, decreasing = TRUE))[1])
}
