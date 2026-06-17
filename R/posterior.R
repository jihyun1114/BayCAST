#' Posterior mode of cluster assignment
#'
#' For each subject, returns the most frequent cluster label in the
#' posterior draws.
#'
#' @param x A fitted `bdpaft` object.
#'
#' @return Integer vector of length `n` (number of subjects).
#' @export
posterior_z <- function(x) {
  stopifnot(inherits(x, "bdpaft"))
  if (is.null(x$z_draws)) stop("`z_draws` not found in fit.")
  apply(x$z_draws, 1, function(zi) {
    tab <- table(zi)
    as.integer(names(which.max(tab)))
  })
}

#' Features selected at a given PIP threshold
#'
#' @param x A fitted `bdpaft` object.
#' @param threshold PIP threshold. Default `0.5`.
#'
#' @return Integer vector of feature indices with PIP above the threshold,
#'   sorted by PIP (descending). Has attribute `"pip"` with the values.
#' @export
posterior_pip <- function(x, threshold = 0.5) {
  stopifnot(inherits(x, "bdpaft"))
  if (is.null(x$pip)) stop("`pip` not found in fit.")
  sel <- which(x$pip > threshold)
  ord <- sel[order(x$pip[sel], decreasing = TRUE)]
  attr(ord, "pip") <- x$pip[ord]
  ord
}

#' Posterior summary of covariate coefficients
#'
#' @param x A fitted `bdpaft` object.
#' @param cov_names Optional character vector of covariate names.
#' @param prob Width of the credible interval. Default `0.95`.
#'
#' @return Data frame with mean, lower CI, upper CI, and significance flag
#'   (CI excludes zero).
#' @export
posterior_beta <- function(x, cov_names = NULL, prob = 0.95) {
  stopifnot(inherits(x, "bdpaft"))
  if (is.null(x$beta_draws)) stop("`beta_draws` not found in fit.")
  q <- nrow(x$beta_draws)
  if (q == 0) return(data.frame())

  alpha <- (1 - prob) / 2
  if (is.null(cov_names)) {
    cov_names <- if (length(x$cov_names) == q) x$cov_names
                 else paste0("beta", seq_len(q))
  }

  means <- rowMeans(x$beta_draws)
  ci_lo <- apply(x$beta_draws, 1, stats::quantile, probs = alpha)
  ci_hi <- apply(x$beta_draws, 1, stats::quantile, probs = 1 - alpha)
  signif <- (ci_lo > 0) | (ci_hi < 0)

  data.frame(
    covariate = cov_names,
    mean      = round(means, 4),
    ci_lo     = round(ci_lo, 4),
    ci_hi     = round(ci_hi, 4),
    signif    = signif,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}
