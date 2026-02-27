library(Rnetrics)

test_that("ols_fit matches lm() coefficients", {
  set.seed(42)
  n <- 1000L
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  beta_true <- c(1.0, 2.0, -1.5)
  Y <- beta_true[1] + beta_true[2] * X$x1 + beta_true[3] * X$x2 + rnorm(n, sd = 0.5)

  fit_ours <- ols_fit(Y, X, silent = TRUE)
  fit_lm   <- lm(Y ~ x1 + x2, data = X)

  expect_equal(as.numeric(fit_ours$beta), as.numeric(coef(fit_lm)), tolerance = 1e-8)
})

test_that("ols_fit HC vcov is close to sandwich::vcovHC", {
  skip_if_not_installed("sandwich")
  set.seed(42)
  n <- 1000L
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y <- 1 + 2 * X$x1 - X$x2 + rnorm(n)

  fit_ours <- ols_fit(Y, X, silent = TRUE)

  # sandwich HC0 for comparison (our formula uses HC1 via fsc = n/(n-K))
  fit_lm   <- lm(Y ~ x1 + x2, data = X)
  vcov_sw  <- sandwich::vcovHC(fit_lm, type = "HC1")

  se_ours <- unname(sqrt(diag(fit_ours$vcov)))
  se_sw   <- unname(sqrt(diag(vcov_sw)))

  # Should agree to 3 significant figures
  expect_equal(se_ours, se_sw, tolerance = 1e-3)
})

test_that("ols_fit cluster-robust SE works without error", {
  set.seed(42)
  n  <- 500L
  cl <- rep(1:50, each = 10L)
  X  <- data.frame(x1 = rnorm(n))
  Y  <- 1 + X$x1 + rnorm(n)

  fit <- ols_fit(Y, X, c_id = cl, silent = TRUE)
  expect_length(as.numeric(fit$beta), 2L)   # constant + x1
  expect_true(all(diag(fit$vcov) >= 0))
})

test_that("ols_fit returns score_i with correct dimensions", {
  set.seed(1)
  n <- 200L
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y <- rnorm(n)

  fit <- ols_fit(Y, X, silent = TRUE)
  expect_equal(dim(fit$score_i), c(n, 3L))   # 3 = constant + 2 regressors
})
