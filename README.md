# baycast

<!-- badges: start -->
<!-- badges: end -->

**Bayesian Dirichlet-Process Accelerated Failure-Time model with factor-analytic spike-and-slab loadings.**

`baycast` fits a survival model that simultaneously (i) clusters patients via a Dirichlet-process mixture on a shared latent factor `U`, and (ii) selects high-dimensional features via spike-and-slab on the loading matrix `A`. Both point-mass and continuous spike variants are supported through a single argument `tau_spike2`.

## Model

The log event time is modeled as

```
log T_i = W_i' β + μ_{z_i} + ε_i,    ε_i ~ N(0, σ²_{z_i})
```

with cluster assignment `z_i` drawn from a Dirichlet process. High-dimensional features follow

```
X_i = U_i A' + e_i
```

with spike-and-slab on the rows of `A`.

## Installation

```r
# install.packages("devtools")
devtools::install_github("jihyun1114/baycast")
```

Requires `Rcpp`, `RcppArmadillo`. Tested on R ≥ 4.2.

## Quick start

```r
library(BayCAST)

set.seed(1)
n <- 100; p <- 200
X <- matrix(rnorm(n * p), n, p)
time   <- rexp(n, rate = 0.1)
status <- rbinom(n, 1, 0.7)

# Point-mass spike (default)
fit <- baycast(time, status, X,
              d = 5,
              iters = 2000, burn = 1000, thin = 5,
              seed = 1)

# Continuous spike-slab
fit_co <- baycast(time, status, X,
                 d = 5,
                 iters = 2000, burn = 1000, thin = 5,
                 tau_spike2 = 0.01,
                 seed = 1)
```

## What the fit returns

| field          | description                                   |
|----------------|-----------------------------------------------|
| `pip`          | posterior inclusion probability per feature   |
| `beta_draws`   | posterior draws for covariate coefficients    |
| `z_draws`      | posterior draws of cluster assignments        |
| `Kplus_draws`  | number of active clusters per iteration       |
| `alpha_draws`  | DP concentration parameter draws              |
| `ssvs_mode`    | `"point_mass"` or `"continuous"`              |
| `tau_spike2`   | spike variance used                           |

## Variable selection and clustering

```r
# Selected features at PIP > 0.5
selected <- which(fit$pip > 0.5)

# Posterior mode of cluster assignment
z_hat <- posterior_z(fit)

# Posterior summary of covariate effects
posterior_beta(fit, cov_names = c("HXMI", "SEX", "RACE_B"))
```

## Diagnostic plots

```r
plot(fit, type = "trace")    # Kplus, alpha, beta trace
plot(fit, type = "pip")      # PIP histogram
plot(fit, type = "survival") # cluster-specific baseline survivor
```

## Reference

Kim, J., et al. (2026). *Bayesian Dirichlet-Process Accelerated Failure-Time Model for Clustering Patient Survival with Variable Selection.* (in preparation).
