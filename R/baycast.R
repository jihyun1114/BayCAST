#' Fit the BayCAST model
#'
#' Bayesian Dirichlet-process accelerated failure-time model with a shared
#' factor-analytic representation of the high-dimensional covariates and
#' spike-and-slab variable selection on the loading matrix.
#'
#' The observed time is modeled on the log scale as
#'
#' \deqn{\log T_i = W_i^\top \beta + \mu_{z_i} + \epsilon_i,
#'       \quad \epsilon_i \sim \mathcal{N}(0, \sigma^2_{z_i}),}
#'
#' with cluster assignment \eqn{z_i} drawn from a Dirichlet process over
#' \eqn{(\mu_k, \sigma^2_k, \nu_k, \Sigma_k)}. The high-dimensional features are
#' modeled as \eqn{X_i = U_i A^\top + e_i} with a spike-and-slab prior on the
#' rows of \eqn{A}. Two variants are supported:
#' \itemize{
#'   \item \strong{Point-mass} (\code{sigma_spike2 = 0}, default): inactive rows
#'     are exactly zero.
#'   \item \strong{Continuous} (\code{sigma_spike2 > 0}): inactive rows are drawn
#'     from a narrow Normal with variance \code{sigma_spike2}.
#' }
#'
#' @param time Numeric vector of length `n`: observed event or censoring time
#'   (`> 0`).
#' @param status Integer vector of length `n`: `1` if event observed, `0` if
#'   right-censored.
#' @param X Numeric matrix `n x p` of high-dimensional features.
#' @param covariates Optional numeric matrix `n x q` of covariates that always
#'   enter the model. Default `NULL`.
#' @param d Latent factor dimension. Default `5`.
#' @param iters Total MCMC iterations. Default `4000`.
#' @param burn Iterations discarded as burn-in. Default `iters %/% 2`.
#' @param thin Thinning interval. Default `10`.
#' @param K_max Truncation level for the DP. Default `50`.
#' @param K_init Number of mixture components at initialization. Default `10`.
#' @param pi0 Prior inclusion probability in SSVS. Default `0.02`.
#' @param sigma_slab2 Slab variance for the SSVS prior on rows of `A`. Default
#'   `1.0`.
#' @param sigma_spike2 Spike variance for continuous spike-slab. Default `0`
#'   (point-mass spike). Set `> 0` (e.g. `0.01`) to use the continuous variant.
#' @param slice Slice sampler for the Dirichlet-process stick-breaking weights.
#'   `"kgw"` (default) is the Kalli-Griffin-Walker slice sampler
#'   (Kalli, Griffin and Walker, 2011, \emph{Statistics and Computing} 21(1),
#'   93--105, \doi{10.1007/s11222-009-9150-y}), which uses an adaptive set of
#'   deterministic decreasing weights; `"walker"` is the original slice sampler
#'   of Walker (2007, \emph{Communications in Statistics} 36(1), 45--54).
#' @param init_A Loading initialization. `"pca"` (default) or `"zero"`.
#' @param init_z Cluster initialization. `"kmeans"` (default) or `"random"`.
#' @param control Advanced control parameters from [baycast_control()].
#' @param seed Optional integer seed.
#' @param verbose Logical; print MCMC progress.
#'
#' @return An object of class `"baycast"`.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 60; p <- 80
#' X <- matrix(rnorm(n * p), n, p)
#' time <- rexp(n, rate = 0.1)
#' status <- rbinom(n, 1, 0.7)
#'
#' # Point-mass spike (default)
#' fit_pm <- baycast(time, status, X, d = 3, iters = 400, burn = 200, thin = 2)
#'
#' # Continuous spike-slab
#' fit_co <- baycast(time, status, X, d = 3, iters = 400, burn = 200, thin = 2,
#'                  sigma_spike2 = 0.01)
#' }
#' @seealso [baycast_control()]
#' @export
baycast <- function(time,
                   status,
                   X,
                   covariates = NULL,
                   d = 5L,
                   iters = 4000L,
                   burn = NULL,
                   thin = 10L,
                   K_max = 50L,
                   K_init = 10L,
                   pi0 = 0.02,
                   sigma_slab2 = 1.0,
                   sigma_spike2 = 0.0,
                   slice = c("kgw", "walker"),
                   init_A = c("pca", "zero"),
                   init_z = c("kmeans", "random"),
                   control = baycast_control(),
                   seed = NULL,
                   verbose = TRUE) {

  slice  <- match.arg(slice)
  init_A <- match.arg(init_A)
  init_z <- match.arg(init_z)

  time   <- as.numeric(time)
  status <- as.integer(status)
  X      <- as.matrix(X)
  n <- length(time)
  p <- ncol(X)
  stopifnot(
    length(status) == n,
    nrow(X) == n,
    all(time > 0),
    all(status %in% c(0L, 1L)),
    d >= 1, K_max >= 2, K_init >= 2, K_init <= K_max,
    iters > 0, thin > 0,
    pi0 > 0, pi0 < 1, sigma_slab2 > 0,
    sigma_spike2 >= 0,
    sigma_spike2 < sigma_slab2 || sigma_spike2 == 0
  )

  if (is.null(burn)) burn <- iters %/% 2L
  burn <- as.integer(burn)
  stopifnot(iters > burn)

  if (!inherits(control, "baycast_control"))
    stop("`control` must be the output of baycast_control().")

  if (is.null(covariates)) {
    W <- matrix(0.0, nrow = n, ncol = 0L)
    q <- 0L
    cov_names <- character(0)
  } else {
    W <- as.matrix(covariates)
    stopifnot(nrow(W) == n)
    q <- ncol(W)
    cov_names <- if (!is.null(colnames(W))) colnames(W)
                 else paste0("X", seq_len(q))
  }

  logC <- log(time)
  mu0       <- if (is.null(control$mu0)) mean(logC) else control$mu0
  b0_eff    <- if (is.null(control$b0))  0.5 * stats::var(logC) else control$b0
  nu_w_pri  <- if (is.null(control$nu_w_prior)) d + 2 else control$nu_w_prior

  b0_beta <- rep(0.0, max(q, 1L))
  B0_beta <- diag(control$beta_prior_var, max(q, 1L))
  if (q == 0L) { b0_beta <- numeric(0); B0_beta <- matrix(0.0, 0, 0) }

  if (is.null(control$diag_feat_idx)) {
    nkeep <- min(50L, p)
    diag_idx <- order(apply(X, 2, stats::var), decreasing = TRUE)[seq_len(nkeep)]
  } else {
    diag_idx <- as.integer(control$diag_feat_idx)
  }

  if (!is.null(seed)) set.seed(seed)

  if (!verbose) {
    sink_con <- textConnection("baycast_silent_log", "w", local = TRUE)
    sink(sink_con, type = "output")
    on.exit({ sink(NULL, type = "output"); close(sink_con) }, add = TRUE)
  }

  res <- baycast_cpp(
    Y              = as.numeric(status),
    logC           = logC,
    X              = X,
    W              = W,
    b0_beta        = b0_beta,
    B0_beta        = B0_beta,
    mu0            = mu0,
    kappa_mu0      = control$kappa_mu0,
    a0             = control$a0,
    b0             = b0_eff,
    d              = as.integer(d),
    K_max          = as.integer(K_max),
    K_init         = as.integer(K_init),
    iters          = as.integer(iters),
    burn           = as.integer(burn),
    thin           = as.integer(thin),
    sigma_slab2        = sigma_slab2,
    sigma_spike2     = sigma_spike2,
    pi0            = pi0,
    alpha_fixed    = control$alpha_fixed,
    alpha_init     = control$alpha_init,
    a_alpha        = control$a_alpha,
    b_alpha        = control$b_alpha,
    sample_ab      = control$sample_ab,
    alpha_a_prior  = control$alpha_a_prior,
    beta_a_prior   = control$beta_a_prior,
    alpha_b_prior  = control$alpha_b_prior,
    beta_b_prior   = control$beta_b_prior,
    mh_a_sd        = control$mh_a_sd,
    init_A_pca     = identical(init_A, "pca"),
    init_z_kmeans  = identical(init_z, "kmeans"),
    ind_slice      = identical(slice, "kgw"),
    rho            = control$rho,
    nu_w_prior     = nu_w_pri,
    diag_feat_idx  = diag_idx,
    diag_max_keep  = control$diag_max_keep
  )

  res$call         <- match.call()
  res$dims         <- c(n = n, p = p, q = q, d = as.integer(d))
  res$cov_names    <- cov_names
  res$slice        <- slice
  res$init_A       <- init_A
  res$init_z       <- init_z
  res$pi0          <- pi0
  res$sigma_slab2      <- sigma_slab2
  res$sigma_spike2   <- sigma_spike2
  res$iters        <- as.integer(iters)
  res$burn         <- as.integer(burn)
  res$thin         <- as.integer(thin)

  class(res) <- c("baycast", "list")
  res
}
