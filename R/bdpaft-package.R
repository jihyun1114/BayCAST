#' bdpaft: Bayesian Dirichlet-Process AFT with Factor-Analytic Spike-and-Slab Loadings
#'
#' A Bayesian semiparametric model for survival data with high-dimensional
#' covariates. The latent factor model is shared across a Dirichlet-process
#' mixture of accelerated-failure-time components; spike-and-slab variable
#' selection (SSVS) is applied to the loading matrix.
#'
#' @section Main function:
#' [bdpaft()] is the user-facing entry point. Use [bdpaft_control()] to
#' adjust advanced hyperparameters.
#'
#' @docType package
#' @name bdpaft-package
#' @aliases bdpaft-package
#' @useDynLib bdpaft, .registration = TRUE
#' @importFrom Rcpp evalCpp
"_PACKAGE"
