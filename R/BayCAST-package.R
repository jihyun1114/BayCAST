#' baycast: Bayesian Dirichlet-Process AFT with Factor-Analytic Spike-and-Slab Loadings
#'
#' A Bayesian semiparametric model for survival data with high-dimensional
#' covariates. The latent factor model is shared across a Dirichlet-process
#' mixture of accelerated-failure-time components; spike-and-slab variable
#' selection (SSVS) is applied to the loading matrix.
#'
#' @section Main function:
#' [baycast()] is the user-facing entry point. Use [baycast_control()] to
#' adjust advanced hyperparameters.
#'
#' @docType package
#' @name baycast-package
#' @aliases baycast-package
#' @useDynLib BayCAST, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom utils tail
"_PACKAGE"
