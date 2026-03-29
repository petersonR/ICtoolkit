# ---------------------------------------------------------------------------
# mBIC tests for lm / glm
# ---------------------------------------------------------------------------

test_that("compute_mbic.lm computes correctly", {
  fit    <- lm(mpg ~ wt + hp, data = mtcars)
  n      <- nobs(fit)
  k      <- length(coef(fit))
  k_pred <- k - 1L
  P      <- 10
  kappa  <- 4
  expected <- as.numeric(stats::BIC(fit)) + 2 * k_pred * log(P / kappa - 1)

  result <- compute_mbic(fit, P = P, kappa = kappa)
  expect_equal(as.numeric(result), expected)
  expect_equal(attr(result, "fit_class"), "lm")
  expect_equal(attr(result, "k"), k)
  expect_equal(attr(result, "P"), P)
  expect_equal(attr(result, "kappa"), kappa)
  expect_equal(attr(result, "criterion"), "mBIC")
})

test_that("compute_mbic.lm with P inferred from fit", {
  fit    <- lm(mpg ~ wt + hp, data = mtcars)
  result <- compute_mbic(fit, kappa = 1)
  # inferred P = 2, kappa = 1 => P > kappa is satisfied
  expect_equal(attr(result, "P"), 2L)
})

test_that("compute_mbic.lm >= compute_bic.lm for P >> kappa", {
  fit <- lm(mpg ~ wt + hp, data = mtcars)
  mbic <- as.numeric(compute_mbic(fit, P = 20, kappa = 4))
  bic  <- as.numeric(compute_bic(fit))
  expect_gt(mbic, bic)
})

test_that("compute_mbic.lm errors when P <= kappa", {
  fit <- lm(mpg ~ wt + hp, data = mtcars)
  expect_error(compute_mbic(fit, P = 3, kappa = 5), "'P' must be greater than 'kappa'")
  expect_error(compute_mbic(fit, P = 5, kappa = 5), "'P' must be greater than 'kappa'")
})

test_that("compute_mbic.lm errors for invalid kappa", {
  fit <- lm(mpg ~ wt + hp, data = mtcars)
  expect_error(compute_mbic(fit, P = 10, kappa = -1), "'kappa' must be a positive")
  expect_error(compute_mbic(fit, P = 10, kappa = 0), "'kappa' must be a positive")
})

test_that("compute_mbic.glm works", {
  fit <- glm(am ~ wt + hp, data = mtcars, family = binomial)
  result <- compute_mbic(fit, P = 10)
  expect_equal(attr(result, "fit_class"), "glm")
  expect_equal(attr(result, "criterion"), "mBIC")
})

# ---------------------------------------------------------------------------
# mBIC2 tests for lm / glm
# ---------------------------------------------------------------------------

test_that("compute_mbic2.lm computes correctly", {
  fit    <- lm(mpg ~ wt + hp, data = mtcars)
  k      <- length(coef(fit))
  k_pred <- k - 1L
  P      <- 10
  kappa  <- 4
  expected <- as.numeric(stats::BIC(fit)) +
    2 * k_pred * log(P / kappa - 1) - 2 * lchoose(P, k_pred)

  result <- compute_mbic2(fit, P = P, kappa = kappa)
  expect_equal(as.numeric(result), expected)
  expect_equal(attr(result, "criterion"), "mBIC2")
  expect_equal(attr(result, "P"), P)
  expect_equal(attr(result, "kappa"), kappa)
})

test_that("compute_mbic2.lm <= compute_mbic.lm (combinatorial offset)", {
  fit  <- lm(mpg ~ wt + hp, data = mtcars)
  mbic  <- as.numeric(compute_mbic(fit, P = 20, kappa = 4))
  mbic2 <- as.numeric(compute_mbic2(fit, P = 20, kappa = 4))
  # mBIC2 subtracts 2*log(C(P,k)), so mBIC2 <= mBIC
  expect_lte(mbic2, mbic)
})

test_that("compute_mbic2.lm equals compute_mbic.lm for intercept-only model", {
  fit   <- lm(mpg ~ 1, data = mtcars)
  # k_pred = 0 => log(C(P, 0)) = 0 => mBIC2 = mBIC
  mbic  <- as.numeric(compute_mbic(fit, P = 10, kappa = 4))
  mbic2 <- as.numeric(compute_mbic2(fit, P = 10, kappa = 4))
  expect_equal(mbic2, mbic)
})

test_that("compute_mbic2.lm errors when P <= kappa", {
  fit <- lm(mpg ~ wt + hp, data = mtcars)
  expect_error(compute_mbic2(fit, P = 3, kappa = 5), "'P' must be greater than 'kappa'")
})

# ---------------------------------------------------------------------------
# mBIC / mBIC2 tests for glmnet
# ---------------------------------------------------------------------------

test_that("compute_mbic.glmnet returns vector of correct length", {
  skip_if_not_installed("glmnet")
  library(glmnet)
  set.seed(42)
  X   <- matrix(rnorm(100 * 10), 100, 10)
  y   <- X[, 1] * 2 + rnorm(100)
  fit <- glmnet(X, y, family = "gaussian")

  result <- compute_mbic(fit, P = 10, kappa = 4)
  expect_length(result, length(fit$lambda))
  expect_equal(attr(result, "fit_class"), "glmnet")
  expect_equal(attr(result, "criterion"), "mBIC")
  expect_equal(attr(result, "P"), 10)
  expect_equal(attr(result, "kappa"), 4)
})

test_that("compute_mbic2.glmnet returns vector of correct length", {
  skip_if_not_installed("glmnet")
  library(glmnet)
  set.seed(42)
  X   <- matrix(rnorm(100 * 10), 100, 10)
  y   <- X[, 1] * 2 + rnorm(100)
  fit <- glmnet(X, y, family = "gaussian")

  result <- compute_mbic2(fit, P = 10, kappa = 4)
  expect_length(result, length(fit$lambda))
  expect_equal(attr(result, "criterion"), "mBIC2")
})

test_that("compute_mbic.glmnet >= compute_bic.glmnet for P >> kappa", {
  skip_if_not_installed("glmnet")
  library(glmnet)
  set.seed(42)
  X   <- matrix(rnorm(100 * 10), 100, 10)
  y   <- X[, 1] + rnorm(100)
  fit <- glmnet(X, y)

  mbic <- as.numeric(compute_mbic(fit, P = 10, kappa = 2))
  bic  <- as.numeric(compute_bic(fit))
  expect_true(all(mbic >= bic - 1e-10))
})

test_that("compute_mbic2.glmnet <= compute_mbic.glmnet", {
  skip_if_not_installed("glmnet")
  library(glmnet)
  set.seed(42)
  X   <- matrix(rnorm(100 * 10), 100, 10)
  y   <- X[, 1] + rnorm(100)
  fit <- glmnet(X, y)

  mbic  <- as.numeric(compute_mbic(fit, P = 10, kappa = 4))
  mbic2 <- as.numeric(compute_mbic2(fit, P = 10, kappa = 4))
  expect_true(all(mbic2 <= mbic + 1e-10))
})

# ---------------------------------------------------------------------------
# mBIC / mBIC2 tests for ncvreg
# ---------------------------------------------------------------------------

test_that("compute_mbic.ncvreg returns vector of correct length", {
  skip_if_not_installed("ncvreg")
  library(ncvreg)
  set.seed(42)
  X   <- matrix(rnorm(100 * 10), 100, 10)
  y   <- X[, 1] * 2 + rnorm(100)
  fit <- ncvreg(X, y, family = "gaussian")

  result <- compute_mbic(fit, P = 10, kappa = 4)
  expect_length(result, length(fit$lambda))
  expect_equal(attr(result, "fit_class"), "ncvreg")
  expect_equal(attr(result, "criterion"), "mBIC")
})

test_that("compute_mbic2.ncvreg returns vector of correct length", {
  skip_if_not_installed("ncvreg")
  library(ncvreg)
  set.seed(42)
  X   <- matrix(rnorm(100 * 10), 100, 10)
  y   <- X[, 1] * 2 + rnorm(100)
  fit <- ncvreg(X, y, family = "gaussian")

  result <- compute_mbic2(fit, P = 10, kappa = 4)
  expect_length(result, length(fit$lambda))
  expect_equal(attr(result, "criterion"), "mBIC2")
})
