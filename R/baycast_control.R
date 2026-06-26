#' Advanced control parameters for [baycast()]
#'
#' Returns a list of hyperparameters that govern the priors and MCMC update
#' details. Most users will not need to change these; pass the result directly
#' to `baycast(..., control = baycast_control(...))`.
#'
#' @param mu0 Prior mean for the cluster-specific intercept. If `NULL`
#'   (default), `baycast()` sets `mu0 = mean(log(time))`.
#' @param kappa_mu0 Prior precision multiplier for `mu_k`. Default `1.0`.
#' @param a0,b0 Shape and rate of the inverse-gamma prior on cluster residual
#'   variance. If `b0` is `NULL`, `baycast()` sets `b0 = 0.5 * var(log(time))`.
#' @param beta_prior_var Variance of the (independent) Normal prior on each
#'   covariate coefficient. Default `10`.
#' @param alpha_init Initial DP concentration. Default `1.0`.
#' @param alpha_fixed If `> 0`, alpha is held fixed at this value; if `< 0`
#'   (default), alpha is sampled.
#' @param a_alpha,b_alpha Shape and rate of the Gamma prior on the DP
#'   concentration. Defaults `2.0`, `2.0` -> prior mean 1.
#' @param sample_ab Logical; if `TRUE` (default), put a hyperprior on
#'   `a_alpha` and `b_alpha`.
#' @param alpha_a_prior,beta_a_prior Gamma hyperprior on `a_alpha`.
#' @param alpha_b_prior,beta_b_prior Gamma hyperprior on `b_alpha`.
#' @param mh_a_sd Proposal SD for the log-`a_alpha` Metropolis-Hastings step.
#'   Default `0.3`.
#' @param nu_w_prior Degrees of freedom for the inverse-Wishart prior on the
#'   cluster covariance. If `NULL` (default), `baycast()` uses `d + 2`.
#' @param rho Geometric tail parameter for the independent slice sampler.
#'   Default `0.7`. Only used when `slice = "walker"`.
#' @param diag_feat_idx Integer vector (1-based) of feature indices to retain
#'   full posterior draws for (loading and idiosyncratic variance). If `NULL`
#'   (default), `baycast()` keeps the top 50 by marginal variance.
#' @param diag_max_keep Maximum number of MCMC draws to store for the
#'   diagnostic feature subset. Default `20`.
#'
#' @return A list of control parameters.
#'
#' @examples
#' ctrl <- baycast_control(beta_prior_var = 100)
#' @export
baycast_control <- function(mu0 = NULL,
                           kappa_mu0 = 1.0,
                           a0 = 4.0,
                           b0 = NULL,
                           beta_prior_var = 10.0,
                           alpha_init = 1.0,
                           alpha_fixed = -1.0,
                           a_alpha = 2.0,
                           b_alpha = 2.0,
                           sample_ab = TRUE,
                           alpha_a_prior = 2.0,
                           beta_a_prior = 1.0,
                           alpha_b_prior = 4.0,
                           beta_b_prior = 2.0,
                           mh_a_sd = 0.3,
                           nu_w_prior = NULL,
                           rho = 0.7,
                           diag_feat_idx = NULL,
                           diag_max_keep = 20L) {

  stopifnot(kappa_mu0 > 0, a0 > 0, beta_prior_var > 0,
            alpha_init > 0, a_alpha > 0, b_alpha > 0,
            is.logical(sample_ab), length(sample_ab) == 1,
            mh_a_sd > 0, rho > 0, rho < 1, diag_max_keep >= 1)

  structure(
    list(
      mu0 = mu0,
      kappa_mu0 = kappa_mu0,
      a0 = a0,
      b0 = b0,
      beta_prior_var = beta_prior_var,
      alpha_init = alpha_init,
      alpha_fixed = alpha_fixed,
      a_alpha = a_alpha,
      b_alpha = b_alpha,
      sample_ab = sample_ab,
      alpha_a_prior = alpha_a_prior,
      beta_a_prior = beta_a_prior,
      alpha_b_prior = alpha_b_prior,
      beta_b_prior = beta_b_prior,
      mh_a_sd = mh_a_sd,
      nu_w_prior = nu_w_prior,
      rho = rho,
      diag_feat_idx = diag_feat_idx,
      diag_max_keep = as.integer(diag_max_keep)
    ),
    class = "baycast_control"
  )
}
