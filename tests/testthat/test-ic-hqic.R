# ---------------------------------------------------------------------------
# HQIC tests for lm / glm
# ---------------------------------------------------------------------------

test_that("compute_hqic.lm computes correctly", {
  fit    <- lm(mpg ~ wt + hp, data = mtcars)
  n      <- nobs(fit)
  k      <- length(coef(fit))
  ll     <- logLik(fit)
  df     <- attr(ll, "df")   # includes sigma for lm
  expected <- -2 * as.numeric(ll) + 2 * df * log(log(n))

  result <- compute_hqic(fit)
  expect_equal(as.numeric(result), expected)
  expect_equal(attr(result, "fit_class"), "lm")
  expect_equal(attr(result, "k"), k)
  expect_equal(attr(result, "criterion"), "HQIC")
})

test_that("compute_hqic.glm computes correctly", {
  fit <- glm(am ~ wt + hp, data = mtcars, family = binomial)
  result <- compute_hqic(fit)
  expect_equal(attr(result, "fit_class"), "glm")
  expect_equal(attr(result, "criterion"), "HQIC")
})

test_that("HQIC lies between AIC and BIC for large enough n", {
  fit  <- lm(mpg ~ wt + hp, data = mtcars)
  aic  <- as.numeric(compute_aic(fit))
  hqic <- as.numeric(compute_hqic(fit))
  bic  <- as.numeric(compute_bic(fit))
  # For n = 32: 2 < 2*log(log(32)) ~ 2.48 < log(32) ~ 3.47
  expect_gt(hqic, aic)
  expect_lt(hqic, bic)
})

# ---------------------------------------------------------------------------
# HQIC tests for glmnet
# ---------------------------------------------------------------------------

test_that("compute_hqic.glmnet returns vector of correct length", {
  skip_if_not_installed("glmnet")
  library(glmnet)
  set.seed(42)
  X   <- matrix(rnorm(100 * 10), 100, 10)
  y   <- X[, 1] * 2 + rnorm(100)
  fit <- glmnet(X, y, family = "gaussian")

  result <- compute_hqic(fit)
  expect_length(result, length(fit$lambda))
  expect_equal(attr(result, "fit_class"), "glmnet")
  expect_equal(attr(result, "criterion"), "HQIC")
  expect_length(attr(result, "lambda"), length(fit$lambda))
})

# ---------------------------------------------------------------------------
# HQIC tests for ncvreg
# ---------------------------------------------------------------------------

test_that("compute_hqic.ncvreg returns vector of correct length", {
  skip_if_not_installed("ncvreg")
  library(ncvreg)
  set.seed(42)
  X   <- matrix(rnorm(100 * 10), 100, 10)
  y   <- X[, 1] * 2 + rnorm(100)
  fit <- ncvreg(X, y, family = "gaussian")

  result <- compute_hqic(fit)
  expect_length(result, length(fit$lambda))
  expect_equal(attr(result, "fit_class"), "ncvreg")
  expect_equal(attr(result, "criterion"), "HQIC")
})
