test_that("baycast runs end-to-end on tiny synthetic data", {
  skip_on_cran()

  set.seed(42)
  n <- 50; p <- 60
  X <- matrix(rnorm(n * p), n, p)
  time   <- rexp(n, rate = 0.1)
  status <- rbinom(n, 1, 0.7)
  W      <- matrix(rnorm(n * 1), n, 1)
  colnames(W) <- "covA"

  fit <- baycast(time, status, X, covariates = W,
                d = 3, iters = 200L, burn = 100L, thin = 2L,
                K_init = 3L, K_max = 10L,
                verbose = FALSE, seed = 1L)

  expect_s3_class(fit, "baycast")
  expect_equal(fit$dims["n"], c(n = n))
  expect_equal(fit$dims["p"], c(p = p))
  expect_equal(fit$dims["q"], c(q = 1L))
  expect_true(fit$kept > 0)
  expect_true(all(fit$pip >= 0 & fit$pip <= 1))
  expect_true(all(fit$Kplus_draws >= 1))
  expect_equal(length(fit$beta_mean), 1L)
})

test_that("summary and print work", {
  skip_on_cran()

  set.seed(7)
  n <- 40; p <- 30
  X <- matrix(rnorm(n * p), n, p)
  time   <- rexp(n)
  status <- rbinom(n, 1, 0.8)

  fit <- baycast(time, status, X, d = 2,
                iters = 100L, burn = 50L, thin = 2L,
                K_init = 3L, K_max = 6L,
                verbose = FALSE, seed = 2L)

  expect_output(print(fit), "BayCAST fit")
  s <- summary(fit)
  expect_s3_class(s, "summary.baycast")
  expect_output(print(s), "Cluster structure")
})

test_that("baycast_control validates inputs", {
  expect_s3_class(baycast_control(), "baycast_control")
  expect_error(baycast_control(kappa_mu0 = -1))
  expect_error(baycast_control(rho = 1.5))
})

test_that("baycast rejects bad inputs", {
  X <- matrix(rnorm(40), 10, 4)
  expect_error(baycast(time = c(-1, rep(1, 9)), status = rep(1, 10), X = X),
               "time > 0")
  expect_error(baycast(time = rep(1, 10), status = rep(2, 10), X = X),
               "status %in%")
})
